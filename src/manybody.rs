pub mod manybody {

use kdtree::KdTree;

pub trait Interaction {
    type Point: AsRef<[f64]> + PartialEq;

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

struct Manybody<'a, T: Interaction> {
    num: usize,
    dimension: usize,
    particles: Vec<T>,
    kdtree: KdTree<f64, usize, &'a T::Point>
}

impl<'a, T: Interaction> Manybody<'a, T> {
    fn new(dimension: usize, particles: Vec<T>) -> Self {
        Manybody { num: particles.len(), dimension, particles, kdtree: KdTree::new(dimension) }
    }

    fn rbmc(
        &mut self, dt: f64, beta: f64, p: usize,
        rng: &mut rand::rngs::ThreadRng
    ) {

    }
}
}
