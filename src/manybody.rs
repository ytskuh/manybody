use rand::Rng;
use rand::prelude::SliceRandom;
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
    fn zero_point() -> Self::Point;
    fn standard_normal(rng: &mut rand::rngs::ThreadRng) -> Self::Point;
    fn standard_uni(rng: &mut rand::rngs::ThreadRng) -> Self::Point;

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
            kdtree.add(particle.point().as_aref(), particle.id()).unwrap();
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

    pub fn rbmc (
        &mut self, dt: f64, beta: f64, omega: f64, p: usize,
        rng: &mut rand::rngs::ThreadRng
    ) -> T::Point
    {
        let i = rng.gen_range(0..self.num);
        let mut xstar_point = self.particles[i].point();
        let mut sum = T::zero_point();

        let range: Vec<usize> = (0..self.num-1).collect();
        for xeta_rawindex in range.choose_multiple(rng, p-1) {
            let mut xeta_index = *xeta_rawindex;
            if *xeta_rawindex >= i { xeta_index += 1; }
            sum += self.particles[i].dw1(&self.particles[xeta_index]);
        }
        sum /= (p-1) as f64;
        xstar_point -= self.particles[i].dv()*dt/(omega*(self.num as f64 - 1.0));
        xstar_point -= sum*dt;
        xstar_point += T::standard_normal(rng) * (2.0*dt / ((self.num as f64-1.0)*omega*omega*beta)).sqrt();

        let xstar = T::new(0, &xstar_point);
        let mut sum = 0f64;
        let found = self.kdtree.within(
            xstar_point.as_aref(),
            T::r_split().powf(2.0),
            &squared_euclidean
        ).unwrap();
        for (_, &index) in found {
            if index != i {
            sum += xstar.w2(&self.particles[index])
            }
        }

        let found = self.kdtree.within(
            self.particles[i].point().as_aref(),
            T::r_split().powf(2.0),
            &squared_euclidean
        ).unwrap();
        for (_, &index) in found {
            if index != i {
            sum -= self.particles[i].w2(&self.particles[index])
            }
        }

        let alpha = (-beta*omega*omega*sum).exp();
        let zeta = rng.gen_range(0.0..1.0);
        if zeta < alpha {
            self.kdtree.remove(&self.particles[i].point().as_aref(), &i).unwrap();
            self.particles[i]=T::new(xstar.id(), &xstar_point);
            self.kdtree.add(self.particles[i].point().as_aref(), i).unwrap();
        }
        self.particles[i].point()
    }

    pub fn mh(&mut self, beta: f64, omega: f64, rng: &mut rand::rngs::ThreadRng) -> T::Point {
        let i = rng.gen_range(0..self.num);
        let xi = &self.particles[i];
        let xstar_point = xi.point() 
        + T::standard_normal(rng);
        let xstar = T::new(0, &xstar_point);
        let mut sum = 0f64;
        for j in 0..self.num {
            if j != i {
                sum += xstar.w(xi) - self.particles[j].w(xi);
            }
        }
        let alpha = (-beta*omega*(xstar.v()-xi.v())-beta*omega*omega*sum).exp();
        let zeta = rng.gen_range(0.0..1.0);
//        println!("{}", zeta < alpha);
        if zeta <= alpha {
            self.particles[i]=T::new(xstar.id(), &xstar_point);
        }
        self.particles[i].point()
    }
}

