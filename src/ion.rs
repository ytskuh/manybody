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

impl AsArrRef<DIM> for PointD {
    fn as_aref(&self) -> &[f64; DIM] {
        &self.data.0[0]
    }
}

const Q_PLUS: f64 = 10.0;
const R_SPILT: f64 = 1.0;
const Q: f64 = 0.1;
const N: usize = 100;


impl Particle<DIM> for Ion {
    type Point = PointD;

    fn id (&self) -> usize { self.id }
    fn point (&self) -> PointD { self.point.clone() }
    fn new_position (&self, point: &Self::Point) -> Self {
        Ion {id: self.id, z: self.z, point: point.clone()}
    }

    fn r_split () -> f64 {
        R_SPILT
    }

    fn v(&self) -> f64 {
        
    }
}
