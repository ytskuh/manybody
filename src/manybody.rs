use rand::rngs::ThreadRng;
use rand::prelude::SliceRandom;
use rand::Rng;
use kiddo::KdTree;
use kiddo::distance::squared_euclidean;
use std::ops::{Add, Sub, AddAssign, SubAssign, Mul, Div, MulAssign, DivAssign};
use std::iter::Sum;

pub trait Particle<const N:usize> {
    type Point: AsRef<[f64; N]> + PartialEq
        + Add<Output = Self::Point>
        + Sub<Output = Self::Point>
        + AddAssign + SubAssign 
        + Mul<f64, Output = Self::Point>
        + Div<f64, Output = Self::Point>
        + MulAssign<f64> + DivAssign<f64>
        + Sum;

    const R_SPLIT: f64;

    fn norm(&self) -> f64;

    fn zero_point() -> Self::Point;
    fn standard_normal<R: Rng>(rng: &mut R) -> Self::Point;
    fn standard_uni(rng: &mut ThreadRng) -> Self::Point;

    fn id (&self) -> usize;
    fn point (&self) -> Self::Point;
    fn new_position (&self, point: &Self::Point) -> Self;

    fn available (&self) -> bool;
    fn v    (&self)                 -> f64;
    fn w    (&self, other: &Self)   -> f64;
    fn w1   (&self, other: &Self)   -> f64;
    fn w2   (&self, other: &Self)   -> f64;

    fn dv   (&self)                 -> Self::Point;
    fn dw   (&self, other: &Self)   -> Self::Point;
    fn dw1  (&self, other: &Self)   -> Self::Point;
    fn dw2  (&self, other: &Self)   -> Self::Point;
}

struct ManybodyParam {
    dt: f64,
    p: usize,
    m: usize,
    c1: f64,
    c2: f64,
    std: f64,
    c3:f64
}

pub struct Manybody<T: Particle<DIM>, const DIM:usize, R: Rng> {
    num: usize,
    particles: Vec<T>,
    kdtree: KdTree<f64, DIM>,
    rng: R,
    param: ManybodyParam,

    pub count: usize
}

impl<T: Particle<N> + Clone, const N:usize, R: Rng> Manybody<T, N, R> {
    pub fn new(particles: Vec<T>, rng: R, dt: f64, p: usize, m: usize, c1: f64, c2: f64, std: f64, c3: f64) -> Self {
        let mut kdtree = KdTree::new();
        for particle in particles.iter() {
            kdtree.add(particle.point().as_ref(), particle.id());
        }
        let param = ManybodyParam {dt, p, m, c1, c2, std, c3};

        Manybody { 
            num: particles.len(),
            particles,
            kdtree,
            rng,
            param,
            count: 0
        }
    }

    pub fn particles(&self) -> &Vec<T> {
        &self.particles
    }

    fn rbmc_step1(&mut self, xi: T) -> T {
        let range: Vec<usize> = (0..self.num-1).collect();
        let sum: <T as Particle<N>>::Point = range
            .choose_multiple(&mut self.rng, self.param.p-1)
            .map(|&xeta_rawindex| {
                if xeta_rawindex < xi.id() {
                    xi.dw1(&self.particles[xeta_rawindex])
                } else {
                    xi.dw1(&self.particles[xeta_rawindex+1])
                }
            }).sum();
        T::new_position(&xi, 
            &(
                xi.point() 
                - (
                    xi.dv()*self.param.c1
                    +sum/(self.param.p-1) as f64*self.param.c2
                )*self.param.dt
                + T::standard_normal(&mut self.rng) 
                * self.param.std * (2.0*self.param.dt).sqrt()
            )
        )
    }

    pub fn rbmc (&mut self) -> &T
    {
        let i = self.rng.gen_range(0..self.num);

        let mut xstar = self.particles[i].clone();
        for _ in 0..self.param.m {
            xstar = self.rbmc_step1(xstar);
        }

        let mut sum = 0f64;
        let found = self.kdtree.within(
            xstar.point().as_ref(),
            T::R_SPLIT.powi(2),
            &squared_euclidean
        );
        sum += found.iter().map(|neighbour| neighbour.item).filter(|&index| index != i).map(|index| xstar.w2(&self.particles[index])).sum::<f64>();

        let found = self.kdtree.within(
            self.particles[i].point().as_ref(),
            T::R_SPLIT.powi(2),
            &squared_euclidean
        );
        // for (_, &index) in found.into() {
        //     if index != i {
        //         sum -= self.particles[i].w2(&self.particles[index]);
        //     }
        // }
        sum -= found.iter().map(|neighbour| neighbour.item).filter(|&index| index != i).map(|index| self.particles[i].w2(&self.particles[index])).sum::<f64>();

        let alpha = (-self.param.c3*sum).exp();
        let zeta = self.rng.gen_range(0.0..1.0);
        if zeta < alpha && xstar.available() {
            self.kdtree.remove(self.particles[i].point().as_ref(), i);
            self.particles[i]=xstar;
            self.kdtree.add(self.particles[i].point().as_ref(), i);
            self.count += 1;
        }
        &self.particles[i]
    }

    pub fn mh(&mut self) -> &T {
        let i = self.rng.gen_range(0..self.num);
        let xi = &self.particles[i];
        let xstar = T::new_position(xi, &(xi.point()+T::standard_normal(&mut self.rng)*self.param.std));
        let mut sum = 0f64;
        for j in 0..self.num {
            if j != i {
                sum += xstar.w(&self.particles[j]) - xi.w(&self.particles[j]);
            }
        }
        let alpha = (-self.param.c1*(xstar.v()-xi.v())-self.param.c2*sum).exp();
        let zeta = self.rng.gen_range(0.0..1.0);
//        println!("{}", zeta < alpha);
        if zeta <= alpha && xstar.available() {
            self.particles[i]=xstar;
            self. count +=1;
        }
        &self.particles[i]
    }
}

