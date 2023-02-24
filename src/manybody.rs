pub mod manybody {

use rand::Rng;
use rand::prelude::SliceRandom;
use kiddo::KdTree;
use kiddo::distance::squared_euclidean;
use std::ops::{Add, Sub, AddAssign, SubAssign, Mul, Div, MulAssign, DivAssign};

use std::collections::HashMap;

fn convert_vec(vec: Vec<(f64, &usize)>) -> Vec<(f64, usize)> {
    vec.iter().map(|(f, &u)| (*f, u)).collect()
}

fn union(vec1: Vec<(f64, usize)>, vec2: Vec<(f64, usize)>) -> Vec<(f64, usize)> {
    let mut map: HashMap<usize, f64> = HashMap::new();

    for elem in vec1.iter().chain(vec2.iter()) {
        let (f, u) = *elem;
        if let Some(max_f) = map.get_mut(&u) {
            if *max_f < f {
                *max_f = f;
            }
        } else {
            map.insert(u, f);
        }
    }

    map.into_iter().map(|(u, f)| (f, u)).collect()
}

pub trait Particle<const N:usize> {
    type Point: AsRef<[f64; N]> + PartialEq
        + Add + Sub + AddAssign + SubAssign 
        + Mul<f64, Output = Self::Point>
        + Div<f64, Output = Self::Point>
        + MulAssign<f64> + DivAssign<f64> + std::fmt::Debug;
    fn zero_point() -> Self::Point;
    fn standard_normal(rng: &mut rand::rngs::ThreadRng) -> Self::Point;

    fn id (&self) -> usize;
    fn point (&self) -> Self::Point;
    fn new (id:usize, point: &Self::Point) -> Self;

    fn r_split () -> f64;

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
    kdtree: KdTree<f64, usize, DIM>
}

impl<T: Particle<N>, const N:usize> Manybody<T, N> {
    pub fn new(particles: Vec<T>) -> Self {
        let mut kdtree = KdTree::new();
        for particle in particles.iter() {
            kdtree.add(particle.point().as_ref(), particle.id()).unwrap();
        }
        Manybody { 
            num: particles.len(),
            particles,
            kdtree
        }
    }

    pub fn particles(&self) -> &Vec<T> {
        &self.particles
    }

    pub fn rbmc(
        &mut self, dt: f64, beta: f64, p: usize,
        rng: &mut rand::rngs::ThreadRng
    ) {
        let i = rng.gen_range(0..self.num);
        let mut xstar_point = self.particles[i].point();
        let mut sum = T::zero_point();

        let range: Vec<usize> = (0..self.num-1).collect();
        for xeta_rawindex in range.choose_multiple(rng, p) {
            let mut xeta_index = *xeta_rawindex;
            if *xeta_rawindex >= i { xeta_index += 1; }
            sum += self.particles[i].dw1(&self.particles[xeta_index]);
        }
        sum /= (p-1) as f64;
        xstar_point -= self.particles[i].dv()*dt;
//        println!("{:?}", self.particles[i].dv()*dt);
        xstar_point -= sum*dt;
//        println!("{:?}", xstar_point);
        xstar_point += T::standard_normal(rng) * (2.0 / ((self.num-1) as f64 * beta)).sqrt();

        let found1 = self.kdtree.within(
            xstar_point.as_ref(),
            T::r_split(),
            &squared_euclidean
        ).unwrap();
        let found2 = self.kdtree.within(
            self.particles[i].point().as_ref(),
            T::r_split(),
            &squared_euclidean
        ).unwrap();
        let found = union(convert_vec(found1), convert_vec(found2));
        let mut sum = 0f64;
        let xstar = T::new(0, &xstar_point);
        for (_, index) in found {
            sum += xstar.w2(&self.particles[index]) - self.particles[i].w2(&self.particles[index]);
        }

        let alpha = (-beta*sum).exp();
        let zeta = rng.gen_range(0.0..1.0);
        if zeta < alpha {
            self.kdtree.remove(&self.particles[i].point().as_ref(), &i).unwrap();
            self.particles[i]=T::new(xstar.id(), &xstar_point);
            self.kdtree.add(self.particles[i].point().as_ref(), i).unwrap();
        }
    }
}
}
