use std::fmt;
use typenum::Unsigned;
extern crate nalgebra as na;
use kdtree::KdTree;
use rand::Rng;
use rand::seq::SliceRandom;
use rand_distr::{StandardNormal};
use array_tool::vec::Union;
use num_traits::float::Float;

type DimType = typenum::U2;
const DIM: usize = DimType::U8 as usize;

type VectorD = na::SVector<f64, DIM>;
trait Particle {
    fn id (&self) -> usize;
    fn point (&self) -> &VectorD;
    fn new (id:usize, point: &VectorD) -> Self;
}
trait Interaction {
    fn r_split () ->f64;
    fn r_update () -> f64;

    fn v    (&self)                 -> f64;
    fn w    (&self, other: &Self)   -> f64;
    fn w1   (&self, other: &Self)   -> f64;
    fn w2   (&self, other: &Self)   -> f64;

    fn dv   (&self)                 -> VectorD;
    fn dw   (&self, other: &Self)   -> VectorD;
    fn dw1  (&self, other: &Self)   -> VectorD;
    fn dw2  (&self, other: &Self)   -> VectorD;
}


struct DBParticle {
    id: usize,
    point: VectorD
}
struct Manybody<T0: Particle + Interaction, A, T: PartialEq, U: AsRef<[A]> + PartialEq>(KdTree<A, T, U>, Vec<T0>);

fn main() {
    // let particle_num = 128;
    // let time_length = 1.0;
    // let group_step_time = 0.1;
    // let single_step_time = group_step_time/particle_num as f64;
    // let step_num = (time_length/single_step_time) as u32;
    // println!("{}", step_num);

    // for step in 0..step_num {
        
    // }
    let y = VectorD::zeros();
    println!("{}", y);
}

impl<T0: Particle + Interaction, A: Float ,T: PartialEq, U: AsRef<[A]> + PartialEq> Manybody<T0, A, T, U> {
    fn rbmc(&mut self, dt: f64, beta: f64, p: usize,
        rng: &mut rand::rngs::ThreadRng) {
        let num = self.1.len();
        let i = rng.gen_range(0..num);
        let mut xstar_point = VectorD::zeros();
        let mut normal = VectorD::zeros();

        xstar_point.copy_from(&self.1[i].point());
        let mut sum = VectorD::zeros();
        for xeta in self.1.choose_multiple(rng, p) {
            sum += self.1[i].dw1(xeta);
        }
        sum/=(p-1) as f64;
        xstar_point-=dt*(self.1[i].dv()/((num-1) as f64));
        for i in 0..DIM {
            normal[i]=rng.sample(StandardNormal);
        }
        xstar_point += (2.0/((num-1) as f64)*beta).sqrt()*normal;
        let xstar = T0::new(self.1[i].id(), &xstar_point);

        let found1 = self.0.within(&xstar, T0::r_update());
        let found2 = self.0.within(&self.1[i], T0::r_update());
        let found = found1.union(found2);
        let mut sum = 0.0;
        for xj in found {
            sum += xstar.w2(xj) - self.1[i].w2(xj);
        }
        let alpha = (-beta*sum).exp();
        let zeta = rng.gen_range(0.0..1.0);
        if zeta < alpha {
            self.1.
            self.1[i]=T0::new(xstar.id(), xstar.point());
        }
    }
}

impl Particle for DBParticle{
    fn id (&self) -> usize {
        self.id
    }

    fn point (&self) -> &VectorD {
        &self.point
    }

    fn new (id:usize, point: &VectorD) -> Self {
        let mut p = VectorD::zeros();
        p.copy_from(point);
        DBParticle{id, point: p}
    }
}

impl PartialEq for DBParticle {
    fn eq(&self, other: &Self) -> bool { self.id == other.id }
    fn ne(&self, other: &Self) -> bool { self.id != other.id }
}

impl Interaction for DBParticle {
    fn r_split () -> f64 {
        0.01
    }

    fn r_update () -> f64 {
        0.01
    }

    fn v(&self) -> f64 {
        self.point.norm_squared()/2.0
    }

    fn w(&self, other: &DBParticle) -> f64 {
        -(self.point-other.point).norm().ln()
    }

    fn w1(&self, other: &DBParticle) -> f64 {
        let r = (self.point-other.point).norm();
        if r>0.01 { return -r.ln(); } 
        100.0_f64.ln()-100.0*r+1.0
    }

    fn w2(&self, other: &DBParticle) -> f64 {
        let r = (self.point-other.point).norm();
        if r>0.01 { return 0.0; }
        -100.0_f64.ln()+100.0*r-1.0   
    }

    fn dv(&self) -> VectorD {
        self.point
    }

    fn dw(&self, other: &DBParticle) -> VectorD {
        let delta = self.point - other.point;
        delta/delta.norm_squared()
    }

    fn dw1(&self, other: &DBParticle) -> VectorD {
        let delta = self.point - other.point;
        let r2 = delta.norm_squared();
        if r2>0.0001 { return -delta/r2; }
        -100.0*delta/r2.sqrt()
    }

    fn dw2(&self, other: &DBParticle) -> VectorD {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r>0.01 {return VectorD::zeros();}
        100.0*delta/r
    }
}

impl std::fmt::Display for DBParticle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.point.fmt(f)
    }
}

