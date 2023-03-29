use rand::rngs::ThreadRng;
use rand::prelude::SliceRandom;
use rand::Rng;
use kiddo::KdTree;
use kiddo::distance::squared_euclidean;
use std::ops::{Add, Sub, AddAssign, SubAssign, Mul, Div, MulAssign, DivAssign};

pub trait AsArrRef<const N:usize> {
    fn as_aref(&self) -> &[f64; N];
}

pub trait Particle<const N:usize> {
    type Point: AsArrRef<N> + PartialEq
        + Add<Output = Self::Point>
        + Sub<Output = Self::Point>
        + AddAssign + SubAssign 
        + Mul<f64, Output = Self::Point>
        + Div<f64, Output = Self::Point>
        + MulAssign<f64> + DivAssign<f64>;

    const R_SPLIT: f64;

    fn norm(&self) -> f64;

    fn zero_point() -> Self::Point;
    fn standard_normal(rng: &mut ThreadRng) -> Self::Point;
    fn standard_uni(rng: &mut ThreadRng) -> Self::Point;

    fn id (&self) -> usize;
    fn point (&self) -> Self::Point;
    fn new_position (&self, point: &Self::Point) -> Self;

    fn reflection (point: &Self::Point) -> Self::Point;
    fn v    (&self)                 -> f64;
    fn w    (&self, other: &Self)   -> f64;
    fn w1   (&self, other: &Self)   -> f64;
    fn w2   (&self, other: &Self)   -> f64;

    fn dv   (&self)                 -> Self::Point;
    fn dw   (&self, other: &Self)   -> Self::Point;
    fn dw1  (&self, other: &Self)   -> Self::Point;
    fn dw2  (&self, other: &Self)   -> Self::Point;
}

pub struct Manybody<T: Particle<DIM>, const DIM:usize> {
    num: usize,
    particles: Vec<T>,
    kdtree: KdTree<f64, usize, DIM>,
    rng: ThreadRng,
    pub count: usize
}

impl<T: Particle<N> + Clone, const N:usize> Manybody<T, N> {
    pub fn new(particles: Vec<T>, rng: ThreadRng) -> Self {
        let mut kdtree = KdTree::new();
        for particle in particles.iter() {
            kdtree.add(particle.point().as_aref(), particle.id()).unwrap();
        }
        Manybody { 
            num: particles.len(),
            particles,
            kdtree,
            rng,
            count: 0
        }
    }

    pub fn particles(&self) -> &Vec<T> {
        &self.particles
    }

    fn rbmc_step1(&mut self, xi: T, dt: f64, p:usize, std: f64, a: f64) -> T {
        let mut sum = T::zero_point();
        let range: Vec<usize> = (0..self.num-1).collect();
        for xeta_rawindex in range.choose_multiple(&mut self.rng, p-1) {
            if *xeta_rawindex < xi.id() {
                sum += xi.dw1(&self.particles[*xeta_rawindex]);
            } else {
                sum += xi.dw1(&self.particles[*xeta_rawindex+1]);
            }
        }
        T::new_position(&xi, &T::reflection(&(xi.point() - (xi.dv()/a+sum/(p-1) as f64)*dt + T::standard_normal(&mut self.rng) * std)))
    }

    pub fn rbmc (
        &mut self, dt: f64, beta: f64, omega: f64, p: usize, m: usize
    ) -> &T
    {
        let i = self.rng.gen_range(0..self.num);
        let a = omega*(self.num as f64 - 1.0);
        let std = (2.0*dt / ((self.num as f64 - 1.0)*omega*omega*beta)).sqrt();
        let mut xstar = self.particles[i].clone();
        for _ in 0..m {
            xstar = self.rbmc_step1(xstar, dt, p, std, a);
        }

        let mut sum = 0f64;
        let found = self.kdtree.within(
            xstar.point().as_aref(),
            T::R_SPLIT.powi(2),
            &squared_euclidean
        ).unwrap();
        for (_, &index) in found {
            if index != i {
                sum += xstar.w2(&self.particles[index]);
            }
        }

        let found = self.kdtree.within(
            self.particles[i].point().as_aref(),
            T::R_SPLIT.powi(2),
            &squared_euclidean
        ).unwrap();
        for (_, &index) in found {
            if index != i {
                sum -= self.particles[i].w2(&self.particles[index]);
            }
        }

        let alpha = (-beta*omega*omega*sum).exp();
        let zeta = self.rng.gen_range(0.0..1.0);
        if zeta < alpha {
            self.kdtree.remove(&self.particles[i].point().as_aref(), &i).unwrap();
            self.particles[i]=xstar;
            self.kdtree.add(self.particles[i].point().as_aref(), i).unwrap();
            self.count += 1;
        }
        &self.particles[i]
    }

    pub fn mh(&mut self, beta: f64, omega: f64) -> &T {
        let i = self.rng.gen_range(0..self.num);
        let xi = &self.particles[i];
        let xstar = T::new_position(xi, &(xi.point()+T::standard_normal(&mut self.rng)));
        let mut sum = 0f64;
        for j in 0..self.num {
            if j != i {
                sum += xstar.w(xi) - self.particles[j].w(xi);
            }
        }
        let alpha = (-beta*omega*(xstar.v()-xi.v())-beta*omega*omega*sum).exp();
        let zeta = self.rng.gen_range(0.0..1.0);
//        println!("{}", zeta < alpha);
        if zeta <= alpha {
            self.particles[i]=xstar;
            self. count +=1;
        }
        &self.particles[i]
    }
}

