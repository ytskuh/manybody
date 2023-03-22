use nalgebra::SVector;
use std::{fmt::Display, f64::consts::PI};
use rand_distr::{StandardNormal, Uniform};

use crate::manybody::{AsArrRef, Particle};

const DIM: usize = 3;

type PointD = SVector<f64, DIM>;

pub struct Ion {
    pub id: usize,
    pub z: f64,
    pub point: PointD
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
const R: f64 = 1.0;
pub const Q: f64 = 0.1;
pub const N: usize = 100;
const R_C: f64 = 0.1;
const SIGMA: f64 = 0.215;

impl Particle<DIM> for Ion {
    type Point = PointD;

    const R_SPLIT: f64 = R_C;

    fn id (&self) -> usize { self.id }
    fn point (&self) -> PointD { self.point.clone() }
    fn new_position (&self, point: &Self::Point) -> Self {
        Ion {id: self.id, z: self.z, point: point.clone()}
    }

    fn zero_point() -> Self::Point { PointD::zeros() }

    fn standard_normal(rng: &mut rand::rngs::ThreadRng) -> Self::Point {
        PointD::from_distribution(&StandardNormal, rng)
    }

    fn standard_uni(rng: &mut rand::rngs::ThreadRng) -> Self::Point {
        PointD::from_distribution(&Uniform::new(-1.0, 1.0), rng)
    }

    fn reflection(point: &Self::Point) -> Self::Point {
        point*R.powi(2)/point.norm_squared()
    }

    fn v(&self) -> f64 {
        self.z*Q_PLUS/(4.0*PI*self.point.norm())
    }

    fn w(&self, other: &Self) -> f64 {
        let r = (self.point-other.point).norm();
        if r > SIGMA {
            self.z*other.z/(4.0*PI*r) 
        } else {
            self.z*other.z/(4.0*PI*SIGMA)
        }
        
    }

    fn w1(&self, other: &Self) -> f64 {
        let r = (self.point-other.point).norm();
        if r > Self::R_SPLIT {
            if r > SIGMA {
                self.z*other.z/(4.0*PI*r) 
            } else {
                self.z*other.z/(4.0*PI*SIGMA)
            }
        } else {
            if Self::R_SPLIT > SIGMA {
                r/Self::R_SPLIT*self.z*other.z/(4.0*PI*Self::R_SPLIT)
            } else {
                r/Self::R_SPLIT*self.z*other.z/(4.0*PI*SIGMA)
            }
        }
    }

    fn w2(&self, other: &Self) -> f64 {
        let r = (self.point-other.point).norm();
        if r > Self::R_SPLIT { return 0.0; }
        if Self::R_SPLIT < SIGMA {
            return self.z*other.z/(4.0*PI*SIGMA) - r/Self::R_SPLIT*self.z*other.z/(4.0*PI*SIGMA);
        }
        if r > SIGMA {
            self.z*other.z/(4.0*PI*r) - r/Self::R_SPLIT*self.z*other.z/(4.0*PI*Self::R_SPLIT)
        } else {
            self.z*other.z/(4.0*PI*SIGMA) - r/Self::R_SPLIT*self.z*other.z/(4.0*PI*Self::R_SPLIT)
        }
    }

    fn dv(&self) -> Self::Point {
        -Q_PLUS*self.z*self.point/(4.0*PI*self.point.norm_squared().powf(1.5))
    }

    fn dw (&self, other: &Self) -> Self::Point {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r > SIGMA {
            -self.z*other.z*delta/(4.0*PI*r.powi(3))
        } else {
            Self::zero_point()
        }
    }

    fn dw1  (&self, other: &Self)   -> Self::Point {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r > Self::R_SPLIT {
            if r > SIGMA {
                -self.z*other.z*delta/(4.0*PI*r.powi(3))
            } else {
                Self::zero_point()
            }
        } else {
            if Self::R_SPLIT > SIGMA {
                self.z*other.z*delta/(4.0*PI*r*Self::R_SPLIT.powi(2))
            } else {
                self.z*other.z*delta/(4.0*PI*r*Self::R_SPLIT*SIGMA)
            }
        }
    }

    fn dw2  (&self, other: &Self)   -> Self::Point {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r > Self::R_SPLIT { return Self::zero_point() }
        if Self::R_SPLIT < SIGMA {
            return -self.z*other.z*delta/(4.0*PI*r*Self::R_SPLIT*SIGMA)
        }
        if r > SIGMA {
            -self.z*other.z*delta/(4.0*PI*r.powi(3)) - self.z*other.z*delta/(4.0*PI*r*Self::R_SPLIT.powi(2))
        } else {
            -self.z*other.z*delta/(4.0*PI*r*Self::R_SPLIT*SIGMA)
        }
    }
}
