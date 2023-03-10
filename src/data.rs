use std::fs::{File, OpenOptions}; 
use std::io::{BufWriter, Write};
use std::fmt::Display;
use ndarray::{ArrayD, ArrayBase, IxDyn};

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

pub fn write_to_file<T: Display> (data: T, file_name: &str)
-> std::io::Result<()> {
    let file = File::create(file_name)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{}", data)?;
    Ok(())
}

pub fn append_to_file<T: Display> (data: T, file_name: &str)
-> std::io::Result<()> {
    let file = OpenOptions::new()
    .create(false).append(true)
    .open(file_name)?;
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
    pub fn new(a: &[f64; N], b: &[f64; N], s: &[usize; N]) -> Self {
        let hist = ArrayBase::zeros(IxDyn(s));
        let mut alpha = [0.0; N];
        let mut beta = [0.0; N];
        let mut shape = [0; N];

        alpha.copy_from_slice(a);
        beta.copy_from_slice(b);
        shape.copy_from_slice(s);
        Histogram {hist, alpha, beta, shape}
    }

    pub fn new_from_shape(other: &Self) -> Self {
        let a = &other.alpha;
        let b = &other.beta;
        let s = &other.shape;
        Self::new(a, b, s)
    }

    pub fn add(&mut self, p: &[f64; N]) {
        let mut index = [0; N];
        for i in 0..N {
            index[i] = ((p[i]-self.alpha[i])/(self.beta[i]-self.alpha[i]) * self.shape[i] as f64).floor() as usize;
        }
        match self.hist.get_mut(IxDyn(&index)) {
            Some(elem) => *elem+=1,
            None => ()
        }
    }

    pub fn hist(&self) -> ArrayD<usize> {
        self.hist.clone()
    }

    pub fn hist_density(&self) -> ArrayD<f64> {
        let mut output = self.hist.clone().mapv(|x| x as f64);
        let mut volume = 1.0;
        for i in 0..N {
            volume *= (self.beta[i]-self.alpha[i])/self.shape[i] as f64;
        }
        output/=output.sum();
        output/=volume;
        output
    }
}
