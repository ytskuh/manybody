use std::fs::{File, OpenOptions}; 
use std::io::{BufWriter, Write};
use std::fmt::Display;
use ndarray::{ArrayD, ArrayBase, IxDyn};

use crate::manybody::AsArrRef;

pub fn append_vec_to_file<T: Display> (data: &Vec<T>, file_name: &str) 
-> std::io::Result<()> {
    let file = OpenOptions::new()
        .create(false).append(true)
        .open(file_name)?;
    let mut writer = BufWriter::new(file);

    for x in data {
        writeln!(writer, "{}", x)?;
    }
    Ok(())
}

pub fn write_str_to_file (data: &str, file_name: &str)
-> std::io::Result<()> {
    let file = File::create(file_name)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{}", data)?;
    Ok(())
}

pub struct Histogram<const N:usize> {
    hist: ArrayD<usize>,
    alpha: [f64; N],
    beta: [f64; N],
    shape: [usize; N]
}

impl<const N:usize> Histogram<N> {
    fn new(a: &[f64; N], b: &[f64; N], s: &[usize; N]) -> Self {
        let hist = ArrayBase::zeros(IxDyn(s));
        let mut alpha = [0.0; N];
        let mut beta = [0.0; N];
        let mut shape = [0; N];

        alpha.copy_from_slice(a);
        beta.copy_from_slice(b);
        shape.copy_from_slice(s);
        Histogram {hist, alpha, beta, shape}
    }

    fn add(mut self, p: &[f64; N]) {
        let mut index = [0; N];
        for i in 0..N {
            index[i] = ((p[i]-self.alpha[i])/(self.beta[i]-self.alpha[i]) * self.shape[i] as f64).floor() as usize;
        }
        let elem = self.hist.get_mut(IxDyn(&index)).unwrap();
        *elem+=1;
    }
}
