use std::collections::VecDeque;
use ndarray::{ArrayD, ArrayView, Axis, concatenate};

use crate::data::Histogram;
use crate::manybody::{Manybody, Particle, AsArrRef};

pub struct AutoRBMC<const N:usize, T: Particle<N>> {
    system: Manybody<T, N>,
    hist: Histogram<N>,
    p: usize,
    dt: f64,
    beta: f64,
    omega: f64,

    pack_step: usize,
    lr: f64
}

impl<const N:usize, T:Particle<N>> AutoRBMC<N, T> {
    pub fn new(
        system: Manybody<T, N>, hist:Histogram<N>,
        p: usize, dt: f64, beta: f64, omega: f64,
        pack_step: usize, lr: f64
    ) -> Self {
        AutoRBMC {system, hist, p, dt, beta, omega, pack_step, lr}
    }

    pub fn sample(&mut self, packs: usize) -> ArrayD<f64> {
        let max_len = 20;
        let mut hist_list = VecDeque::new();
        let mut h = Histogram::new_from_shape(&self.hist);
        let mut variance_list = VecDeque::new();
        let mut min_variance = f64::MAX;
        for _ in 0..packs {
            h = Histogram::new_from_shape(&self.hist);
            for _ in 0..self.pack_step {
                let x = self.system.rbmc(self.dt, self.beta, self.omega, self.p);
                h.add(x.as_aref());
            }
            hist_list.push_back(h.hist_density());
            if hist_list.len() > max_len {
                hist_list.pop_front();

                let hist_list_slice: &[ArrayView<f64, _>] = &hist_list.iter().map(|a| a.view()).collect::<Vec<_>>();
                let stacked_array = concatenate(Axis(0), hist_list_slice).unwrap();
                let mean = stacked_array.mean_axis(Axis(0)).unwrap();
                let variance = stacked_array.mapv(|x| x.powi(2)) - mean.mapv(|x| x.powi(2));

                let sum_of_variance = variance.sum();
                variance_list.push_back(sum_of_variance);
                if sum_of_variance < min_variance { min_variance = sum_of_variance }
            }

            if variance_list.len() > max_len {
                variance_list.pop_front();
                let mut temp = true;
                for v in variance_list.iter() {
                    if *v == min_variance {
                        temp = false;
                    }
                }
                if temp {
                    self.dt*=self.lr;
                }
            }
        }
        h.hist_density()
    }

}