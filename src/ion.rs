use nalgebra::SVector;
use std::{fmt::Display, f64::consts::PI};
use rand_distr::{StandardNormal, Uniform};

use crate::manybody::{AsArrRef, Particle};

pub const DIM: usize = 3;

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


pub static QF: f64 = 10.0;
const R: f64 = 1.0;
const V: f64 = 0.1;

const R_C: f64 = 0.3;
static SIGMA: f64 = 0.014;

const R2: f64 = R*R;


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

    fn norm(&self) -> f64 {
        self.point.norm()
    }

    fn available (&self) -> bool {
        let r2 = self.point.norm_squared();
        if r2 > R2 { true }
        else { false }
    }

    fn v(&self) -> f64 {
        self.z*QF/(4.0*PI*self.point.norm())/V
    }

    fn w(&self, other: &Self) -> f64 {
        let r = (self.point-other.point).norm();
        if r > SIGMA {
            self.z*other.z/(4.0*PI*r)/V
        } else {
            self.z*other.z/(4.0*PI*SIGMA)/V
        }
    }

    fn w1(&self, other: &Self) -> f64 {
        let r = (self.point-other.point).norm();
        if r > Self::R_SPLIT {
            if r > SIGMA {
                self.z*other.z/(4.0*PI*r)/V
            } else {
                self.z*other.z/(4.0*PI*SIGMA)/V
            }
        } else {
            if Self::R_SPLIT > SIGMA {
                r/Self::R_SPLIT*self.z*other.z/(4.0*PI*Self::R_SPLIT)/V
            } else {
                r/Self::R_SPLIT*self.z*other.z/(4.0*PI*SIGMA)/V
            }
        }
    }

    fn w2(&self, other: &Self) -> f64 {
        let r = (self.point-other.point).norm();
        if r > Self::R_SPLIT { return 0.0; }
        if Self::R_SPLIT < SIGMA {
            return self.z*other.z/(4.0*PI*SIGMA) * (1.0 - r/Self::R_SPLIT)/V;
        }
        if r > SIGMA {
            (self.z*other.z/(4.0*PI)*(1.0/r - r/Self::R_SPLIT.powi(2)))/V
        } else {
            (self.z*other.z/(4.0*PI)*(1.0/SIGMA - r/Self::R_SPLIT.powi(2)))/V
        }
    }

    fn dv(&self) -> Self::Point {
        self.point*(-QF*self.z/(4.0*PI*self.point.norm_squared().powf(1.5))/V)
    }

    fn dw (&self, other: &Self) -> Self::Point {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r > SIGMA {
            delta*(-self.z*other.z/(4.0*PI*r.powi(3))/V)
        } else {
            Self::zero_point()
        }
    }

    fn dw1  (&self, other: &Self)   -> Self::Point {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r > Self::R_SPLIT {
            if r > SIGMA {
                delta*(-self.z*other.z/(4.0*PI*r.powi(3))/V)
            } else {
                Self::zero_point()
            }
        } else {
            if Self::R_SPLIT > SIGMA {
                delta*(self.z*other.z/(4.0*PI*r*Self::R_SPLIT.powi(2))/V)
            } else {
                delta*(self.z*other.z/(4.0*PI*r*Self::R_SPLIT*SIGMA)/V)
            }
        }
    }

    fn dw2  (&self, other: &Self)   -> Self::Point {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r > Self::R_SPLIT { return Self::zero_point() }
        if Self::R_SPLIT < SIGMA {
            return delta*(-self.z*other.z/(4.0*PI*r*Self::R_SPLIT*SIGMA)/V)
        }
        if r > SIGMA {
            delta*(-(self.z*other.z/(4.0*PI)*(1.0/r.powi(3) + 1.0/r*Self::R_SPLIT.powi(2)))/V)
        } else {
            delta*(-self.z*other.z/(4.0*PI*r*Self::R_SPLIT.powi(2))/V)
        }
    }
}

impl Clone for Ion {
    fn clone(&self) -> Self {
        Ion {id: self.id, z: self.z, point: self.point.clone()}
    }
}

impl Display for Ion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, num) in <PointD as AsRef<[f64; DIM]>>::as_ref(&self.point).iter().enumerate() {
            if i>0 { write!(f, ",")?; }
            write!(f, "{}", num)?;
        }
        Ok(())
    }
}
