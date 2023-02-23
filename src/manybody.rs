pub mod manybody {

use rand::Rng;
use rand::prelude::SliceRandom;
use kdtree::KdTree;
use kdtree::distance::squared_euclidean;
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

pub trait Particle {
    type Point: AsRef<[f64]> + PartialEq
        + Add + Sub + AddAssign + SubAssign 
        + Mul<f64, Output = Self::Point>
        + Div<f64, Output = Self::Point>
        + MulAssign<f64> + DivAssign<f64>;
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

pub struct Manybody<T: Particle> {
    num: usize,
    particles: Vec<T>,
    kdtree: KdTree<f64, usize, T::Point>
}

impl<T: Particle> Manybody<T> {
    pub fn new(dimension: usize, particles: Vec<T>) -> Self {
        let mut kdtree = KdTree::new(dimension);
        for particle in particles.iter() {
            kdtree.add(particle.point(), particle.id()).unwrap();
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

        for xeta in self.particles.choose_multiple(rng, p) {
            sum += self.particles[i].dw1(xeta);
        }
        sum /= (p-1) as f64;
        xstar_point -= (self.particles[i].dv() / (self.num-1) as f64)*dt;
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
            self.kdtree.remove(&self.particles[i].point(), &i).unwrap();
            self.particles[i]=T::new(xstar.id(), &xstar_point);
            self.kdtree.add(self.particles[i].point(), i).unwrap();
        }
    }
}
}

// impl<T0: Particle + Interaction, A: Float ,T: PartialEq, U: AsRef<[A]> + PartialEq> Manybody<T0, A, T, U> {
//     fn rbmc(&mut self, dt: f64, beta: f64, p: usize,
//         rng: &mut rand::rngs::ThreadRng) {
//         let num = self.1.len();
//         let i = rng.gen_range(0..num);
//         let mut xstar_point = VectorD::zeros();
//         let mut normal = VectorD::zeros();

//         xstar_point.copy_from(&self.1[i].point());
//         let mut sum = VectorD::zeros();
//         for xeta in self.1.choose_multiple(rng, p) {
//             sum += self.1[i].dw1(xeta);
//         }
//         sum/=(p-1) as f64;
//         xstar_point-=dt*(self.1[i].dv()/((num-1) as f64));
//         for i in 0..DIM {
//             normal[i]=rng.sample(StandardNormal);
//         }
//         xstar_point += (2.0/((num-1) as f64)*beta).sqrt()*normal;
//         let xstar = T0::new(self.1[i].id(), &xstar_point);

//         let found1 = self.0.within(&xstar, T0::r_update());
//         let found2 = self.0.within(&self.1[i], T0::r_update());
//         let found = found1.union(found2);
//         let mut sum = 0.0;
//         for xj in found {
//             sum += xstar.w2(xj) - self.1[i].w2(xj);
//         }
//         let alpha = (-beta*sum).exp();
//         let zeta = rng.gen_range(0.0..1.0);
//         if zeta < alpha {
//             self.1.
//             self.1[i]=T0::new(xstar.id(), xstar.point());
//         }
//     }
// }