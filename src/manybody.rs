pub mod manybody {

use kdtree::KdTree;

pub trait Particle {
    type Point: AsRef<[f64]> + PartialEq;

    fn id (&self) -> usize;
    fn point (&self) -> Self::Point;
    fn new (id:usize, point: &Self::Point) -> Self;

    fn r_split () -> f64;
    fn r_update () -> f64;

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
    dimension: usize,
    kdtree: KdTree<f64, usize, T::Point>
}

impl<T: Particle> Manybody<T> {
    fn new(dimension: usize, particles: Vec<T>) -> Self {
        let mut item = Manybody { 
            num: particles.len(),
            dimension,
            kdtree: KdTree::new(dimension) 
        };
        for particle in particles.iter() {
            item.kdtree.add(particle.point(), particle.id()).unwrap();
        }
        item
    }

    fn rbmc(
        &mut self, dt: f64, beta: f64, p: usize,
        rng: &mut rand::rngs::ThreadRng
    ) {

    }
}
}
