use nalgebra::SVector;
use std::{fmt::Display, f64::consts::PI};
use rand_distr::{StandardNormal, Uniform};

use crate::manybody::{AsArrRef, Particle};

pub const WIDTH: usize = 3;

struct Theta {
    alpha: SVector<f64, WIDTH>,
    w: SVector<f64, WIDTH>,
    b: SVector<f64, WIDTH>
}