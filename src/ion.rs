use nalgebra::SVector;
use std::fmt::Display;
use rand_distr::{StandardNormal, Uniform};

use crate::manybody::{AsArrRef, Particle};

const DIM: usize = 3;

type PointD = SVector<f64, DIM>;

pub struct Ion {
    id: usize,
    z: f64,
    point: PointD
}

impl PartialEq for Ion {
    fn eq(&self, other: &Self) -> bool { self.id == other.id }
    fn ne(&self, other: &Self) -> bool { self.id != other.id }
}

impl AsArrRef<3> for PointD {
    fn as_aref(&self) -> &[f64; DIM] {
        &self.data.0[0]
    }
}

// impl Particle for Ion {
//     type Point = PointD;

//     fn id (&self) -> usize { self.id }
//     fn point (&self) -> PointD { self.point.clone() }
//     fn new (id:usize, point: &Self::Point) -> Self {
        
//     }

// }
