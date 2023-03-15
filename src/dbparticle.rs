use nalgebra::SVector;
use std::fmt::Display;
use rand_distr::{StandardNormal, Uniform};

use crate::manybody::{AsArrRef, Particle};

const DIM: usize = 1;

type VectorD = SVector<f64, DIM>;

type PointD = VectorD;

pub struct DBParticle {
    pub id: usize,
    pub point: PointD
}

impl PartialEq for DBParticle {
    fn eq(&self, other: &Self) -> bool { self.id == other.id }
    fn ne(&self, other: &Self) -> bool { self.id != other.id }
}

impl AsArrRef<1> for PointD {
    fn as_aref(&self) -> &[f64; DIM] {
        &self.data.0[0]
    }
}

impl Display for DBParticle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, num) in <PointD as AsRef<[f64; DIM]>>::as_ref(&self.point).iter().enumerate() {
            if i>0 { write!(f, ",")?; }
            write!(f, "{}", num)?;
        }
        Ok(())
    }
}

impl Particle<DIM> for DBParticle {
    type Point = PointD;

    fn zero_point() -> Self::Point {
        VectorD::zeros()
    }

    fn standard_normal(rng: &mut rand::rngs::ThreadRng) -> Self::Point {
        VectorD::from_distribution(&StandardNormal, rng)
    }

    fn standard_uni(rng: &mut rand::rngs::ThreadRng) -> Self::Point {
        VectorD::from_distribution(&Uniform::new(-1.0, 1.0), rng)
    }

    fn id(&self) -> usize {
        self.id
    }

    fn point(&self) -> PointD {
        self.point.clone()
    }

    fn new_position(&self, point: &PointD) -> Self {
        DBParticle {id: self.id, point: point.clone()}
    }

    fn r_split () -> f64 { 0.01 }  

    fn v(&self) -> f64 {
        self.point[0].powi(2)/2.0
    }

    fn w(&self, other: &DBParticle) -> f64 {
        -(self.point[0]-other.point[0]).abs().ln()
    }

    fn w1(&self, other: &DBParticle) -> f64 {
        let r = (self.point[0]-other.point[0]).abs();
        if r>0.01 { return -r.ln(); } 
        100.0_f64.ln()-100.0*r+1.0
    }

    fn w2(&self, other: &DBParticle) -> f64 {
        let r = (self.point[0]-other.point[0]).abs();
        if r>0.01 { 
            panic!(); }
        -r.ln()-100.0_f64.ln()+100.0*r-1.0   
    }

    fn dv(&self) -> PointD {
        self.point.clone()
    }

    fn dw(&self, other: &DBParticle) -> PointD {
        let delta = self.point - other.point;
        -delta/delta.norm_squared()
    }

    fn dw1(&self, other: &DBParticle) -> PointD {
        let delta = self.point - other.point;
        let r2 = delta.norm_squared();
        if r2>0.0001 { return -delta/r2; }
        -100.0*delta/r2.sqrt()
    }

    fn dw2(&self, other: &DBParticle) -> PointD {
        let delta = self.point - other.point;
        let r = delta[0].abs();
        if r>0.01 { return VectorD::zeros(); }
        delta*(100.0/r - 1.0/(r*r))
    }
}

impl Clone for DBParticle {
    fn clone(&self) -> Self {
        DBParticle { id: self.id, point: self.point.clone() }
    }
}
