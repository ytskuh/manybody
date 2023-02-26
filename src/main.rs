mod manybody;

extern crate nalgebra as na;
use std::fs::OpenOptions; 
use std::io::{BufWriter, Write};
use std::fmt::Display;
use std::ops::*;
extern crate derive_more;
use derive_more::{Add, Sub, AddAssign, SubAssign};


use rand_distr::StandardNormal;
use manybody::{AsArrRef, Particle, Manybody};

const DIM: usize = 1;

type VectorD = na::SVector<f64, DIM>;

type PointD = VectorD;

struct DBParticle {
    id: usize,
    point: PointD
}

fn main() {
    let particle_num = 500;
    let time_length = 1f64;
    let step_time = 0.0001;
//    let step_num = (particle_num as f64*time_length/step_time) as u32;
    let step_num = 1e7 as u32;
    let skip_step_num = 3000000;
    let beta = 1f64;
    let p = 2;
    let filename = "result5.csv";
    println!("{}", step_num);

    let mut rng = rand::thread_rng();

    let mut particle_initial = Vec::new();
    for i in 0..particle_num {
        particle_initial.push(DBParticle {
            id: i,
            point: DBParticle::standard_normal(&mut rng)
        });
    }
    let mut particle_system = Manybody::new(particle_initial);

    write_str_to_file("x", filename, true).unwrap();

    for i in 0..step_num {
        particle_system.rbmc(step_time, beta/ (particle_num-1) as f64, p, &mut rng);
        if (i%(particle_num*3) as u32) == 0 {
            write_vec_to_file(particle_system.particles(), filename, true).unwrap();
        }
    }
}

fn write_vec_to_file<T: Display> (data: &Vec<T>, file_name: &str, append: bool) 
-> std::io::Result<()> {
    let file = OpenOptions::new()
        .create(true).append(append)
        .open(file_name)?;
    let mut writer = BufWriter::new(file);

    for x in data {
        writeln!(writer, "{}", x)?;
    }
    Ok(())
}

fn write_str_to_file (data: &str, file_name: &str, append: bool)
-> std::io::Result<()> {
    let file = OpenOptions::new()
    .create(true).append(append)
    .open(file_name)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{}", data)?;
    Ok(())
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

    fn id(&self) -> usize {
        self.id
    }

    fn point(&self) -> PointD {
        self.point.clone()
    }

    fn new(id: usize, point: &PointD) -> Self {
        DBParticle {id, point: point.clone()}
    }

    fn r_split () -> f64 { 0.01 }  

    fn v(&self) -> f64 {
        self.point.norm_squared()/2.0
    }

    fn w(&self, other: &DBParticle) -> f64 {
        -(self.point-other.point).norm().ln()
    }

    fn w1(&self, other: &DBParticle) -> f64 {
        let r = (self.point-other.point).norm();
        if r>0.01 { return -r.ln(); } 
        100.0_f64.ln()-100.0*r+1.0
    }

    fn w2(&self, other: &DBParticle) -> f64 {
        let r = (self.point-other.point).norm();
        if r>0.01 { return 0.0; }
        -100.0_f64.ln()+100.0*r-1.0   
    }

    fn dv(&self) -> PointD {
        self.point.clone()
    }

    fn dw(&self, other: &DBParticle) -> PointD {
        let delta = self.point - other.point;
        delta/delta.norm_squared()
    }

    fn dw1(&self, other: &DBParticle) -> PointD {
        let delta = self.point - other.point;
        let r2 = delta.norm_squared();
        if r2>0.0001 { return -delta/r2; }
        -100.0*delta/r2.sqrt()
    }

    fn dw2(&self, other: &DBParticle) -> PointD {
        let delta = self.point - other.point;
        let r = delta.norm();
        if r>0.01 { return VectorD::zeros(); }
        100.0*delta/r
    }
}