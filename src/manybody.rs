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
            dimension,
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