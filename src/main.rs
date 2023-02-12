use kd_tree::KdTree;
use rand;

extern crate nalgebra as na;

const DIM: usize = 1;

type Particle = na::SVector<f64, DIM>;
pub trait Interaction {
    fn v_field(&self) -> f64;
    fn w_field(&self, other: &Self) -> f64;
    fn dv(&self) -> Self;
    fn dw(&self, other: &Self) -> Self;
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
    let x = Particle::new(1.0);
    let mut y = x.dv();
    y = Particle::new(2.0);
    println!("{}", y);
}

impl Interaction for Particle {
    fn v_field(&self) -> f64 {
        self.norm_squared()/2.0
    }

    fn w_field(&self, other: &Particle) -> f64 {
        (self-other).norm().ln()
    }

    fn dv(&self) -> Particle {
        *self
    }

    fn dw(&self, other: &Particle) -> Particle {
        let delta = self - other;
        delta/delta.norm_squared()
    }
}
