use std::fmt;
use typenum::Unsigned;
extern crate nalgebra as na;
use kd_tree::{KdPoint, KdTree};

type DimType = typenum::U2;
const DIM: usize = DimType::U8 as usize;

type VectorD = na::SVector<f64, DIM>;
struct Particle(VectorD);
struct Manybody<T: KdPoint>(KdTree<T>);
pub trait Interaction {
    fn v_field(&self) -> f64;
    fn w_field(&self, other: &Self) -> f64;
    fn dv(&self) -> VectorD;
    fn dw(&self, other: &Self) -> VectorD;
}
fn main() {
    // let particle_num = 128;
    // let time_length = 1.0;
    // let group_step_time = 0.1;
    // let single_step_time = group_step_time/particle_num as f64;
    // let step_num = (time_length/single_step_time) as u32;
    // println!("{}", step_num);

    // for step in 0..step_num {
        
    // }
    let x = VectorD::new(1.0, 2.0);
    let y = Particle::new(x);
    println!("{}", y);
}

impl<T: KdPoint> Manybody<T> {
    fn rbmc(&self) {
        
    }
}

impl Particle {
    fn new(some: VectorD) -> Self {
        Self(some)
    }
}

impl Interaction for Particle {
    fn v_field(&self) -> f64 {
        self.0.norm_squared()/2.0
    }

    fn w_field(&self, other: &Particle) -> f64 {
        (self.0-other.0).norm().ln()
    }

    fn dv(&self) -> VectorD {
        self.0
    }

    fn dw(&self, other: &Particle) -> VectorD {
        let delta = self.0 - other.0;
        delta/delta.norm_squared()
    }
}

impl std::fmt::Display for Particle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl KdPoint for Particle {
    type Scalar = f64;
    type Dim = DimType;
    fn at (&self, k: usize) -> f64 {
        self.0[k]
    }
}
