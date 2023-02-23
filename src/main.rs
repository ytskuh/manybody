mod manybody;

extern crate nalgebra as na;
use std::fs::OpenOptions; 
use std::io::{BufWriter, Write};
use std::fmt::Display;
use std::ops::*;
extern crate derive_more;
use derive_more::{Add, Sub, AddAssign, SubAssign};


use rand_distr::StandardNormal;
use crate::manybody::manybody::Particle;
use crate::manybody::manybody::Manybody;

const DIM: usize = 1;

type VectorD = na::SVector<f64, DIM>;

#[derive(Clone, PartialEq, Add, Sub, AddAssign, SubAssign)]
struct PointD(VectorD);

struct DBParticle {
    id: usize,
    point: PointD
}

fn main() {
    let particle_num = 128;
    let time_length = 100.0;
    let group_step_time = 0.1;
    let single_step_time = group_step_time/particle_num as f64;
    let step_num = (time_length/single_step_time) as u32;
    let beta = 1f64;
    let p = 2;
    println!("{}", step_num);

    let mut rng = rand::thread_rng();

    let mut particle_initial = Vec::new();
    for i in 0..particle_num {
        particle_initial.push(DBParticle {
            id: i,
            point: PointD(VectorD::from_distribution(&StandardNormal, &mut rng))
        });
    }
    let mut particle_system = Manybody::new(DIM, particle_initial);

    write_str_to_file("x", "result.csv", true).unwrap();

    for i in 0..step_num {
        particle_system.rbmc(single_step_time, beta, p, &mut rng);
        if i>1000 && (i%particle_num as u32) == 0 {
            write_vec_to_file(particle_system.particles(), "result.csv", true).unwrap();
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

impl AsRef<[f64]> for PointD {
    fn as_ref(&self) -> &[f64] {
        &self.0.data.0[0]
    }
}

impl Display for DBParticle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, num) in self.point.as_ref().iter().enumerate() {
            if i>0 { write!(f, ",")?; }
            write!(f, "{}", num)?;
        }
        Ok(())
    }
}

impl Particle for DBParticle {
    type Point = PointD;

    fn zero_point() -> Self::Point {
        PointD(VectorD::zeros())
    }

    fn standard_normal(rng: &mut rand::rngs::ThreadRng) -> Self::Point {
        PointD(VectorD::from_distribution(&StandardNormal, rng))
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
        self.point.0.norm_squared()/2.0
    }

    fn w(&self, other: &DBParticle) -> f64 {
        -(self.point.0-other.point.0).norm().ln()
    }

    fn w1(&self, other: &DBParticle) -> f64 {
        let r = (self.point.0-other.point.0).norm();
        if r>0.01 { return -r.ln(); } 
        100.0_f64.ln()-100.0*r+1.0
    }

    fn w2(&self, other: &DBParticle) -> f64 {
        let r = (self.point.0-other.point.0).norm();
        if r>0.01 { return 0.0; }
        -100.0_f64.ln()+100.0*r-1.0   
    }

    fn dv(&self) -> PointD {
        self.point.clone()
    }

    fn dw(&self, other: &DBParticle) -> PointD {
        let delta = self.point.0 - other.point.0;
        PointD(delta/delta.norm_squared())
    }

    fn dw1(&self, other: &DBParticle) -> PointD {
        let delta = self.point.0 - other.point.0;
        let r2 = delta.norm_squared();
        if r2>0.0001 { return PointD(-delta/r2); }
        PointD(-100.0*delta/r2.sqrt())
    }

    fn dw2(&self, other: &DBParticle) -> PointD {
        let delta = self.point.0 - other.point.0;
        let r = delta.norm();
        if r>0.01 { return PointD(VectorD::zeros()); }
        PointD(100.0*delta/r)
    }
}

impl Mul<f64> for PointD {
    type Output = PointD;

    fn mul(self, rhs: f64) -> PointD {
        PointD(self.0*rhs)
    }
}

impl Div<f64> for PointD {
    type Output = PointD;
    fn div(self, rhs: f64) -> Self::Output {
        PointD(self.0/rhs)
    }
}

impl MulAssign<f64> for PointD {
    fn mul_assign(&mut self, rhs: f64) {
        self.0 *= rhs;
    }
}

impl DivAssign<f64> for PointD {
    fn div_assign(&mut self, rhs: f64) {
        self.0 /= rhs;
    }
}