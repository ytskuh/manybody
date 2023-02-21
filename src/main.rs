mod manybody;

extern crate nalgebra as na;
use rand::Rng;
use rand::seq::SliceRandom;
use rand_distr::{StandardNormal};
use array_tool::vec::Union;
use num_traits::float::Float;
use crate::manybody::manybody::Interaction;

const DIM: usize = 2;

type VectorD = na::SVector<f64, DIM>;

#[derive(Clone)]
struct PointD(VectorD);
trait Particle {
    fn id (&self) -> usize;
    fn point (&self) -> &VectorD;
    fn new (id:usize, point: &VectorD) -> Self;
}

struct DBParticle {
    id: usize,
    point: PointD
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
    let y = VectorD::zeros();
    println!("{}", y);
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

// impl Particle for DBParticle{
//     fn id (&self) -> usize {
//         self.id
//     }

//     fn point (&self) -> &VectorD {
//         &self.point
//     }

//     fn new (id:usize, point: &VectorD) -> Self {
//         let mut p = VectorD::zeros();
//         p.copy_from(point);
//         DBParticle{id, point: p}
//     }
// }

impl PartialEq for DBParticle {
    fn eq(&self, other: &Self) -> bool { self.id == other.id }
    fn ne(&self, other: &Self) -> bool { self.id != other.id }
}

impl AsRef<[f64]> for PointD {
    fn as_ref(&self) -> &[f64] {
        &self.0.data.0[0]
    }
}

impl PartialEq for PointD {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl Interaction for DBParticle {
    type Point = PointD;

    fn r_split () -> f64 {
        0.01
    }

    fn r_update () -> f64 {
        0.01
    }

    fn v(&self) -> f64 {
        self.point.0.norm_squared()/2.0
    }

    fn w(&self, other: &DBParticle) -> f64 {
        -(self.point.0-other.point.0).norm().ln()
    }

    fn w1(&self, other: &DBParticle) -> f64 {
        let r = (self.point.0-other.point.0).norm();
        if r>0.01 { return -r.ln(); } 
        100.0_f64.ln()-100.0*r+1.0
    }

    fn w2(&self, other: &DBParticle) -> f64 {
        let r = (self.point.0-other.point.0).norm();
        if r>0.01 { return 0.0; }
        -100.0_f64.ln()+100.0*r-1.0   
    }

    fn dv(&self) -> PointD {
        self.point.clone()
    }

    fn dw(&self, other: &DBParticle) -> PointD {
        let delta = self.point.0 - other.point.0;
        PointD(delta/delta.norm_squared())
    }

    fn dw1(&self, other: &DBParticle) -> PointD {
        let delta = self.point.0 - other.point.0;
        let r2 = delta.norm_squared();
        if r2>0.0001 { return PointD(-delta/r2); }
        PointD(-100.0*delta/r2.sqrt())
    }

    fn dw2(&self, other: &DBParticle) -> PointD {
        let delta = self.point.0 - other.point.0;
        let r = delta.norm();
        if r>0.01 {return PointD(VectorD::zeros());}
        PointD(100.0*delta/r)
    }
}



