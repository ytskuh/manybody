use crate::data::Histogram;
use crate::manybody::{Manybody, Particle};

pub struct AutoRBMC<const N:usize, T: Particle<N>> {
    system: Manybody<T, N>,
    hist: Histogram<N>
}