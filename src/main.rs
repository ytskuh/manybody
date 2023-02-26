mod manybody;
mod particle;
mod data;

extern crate nalgebra as na;

use crate::manybody::{Particle, Manybody};
use particle::{DBParticle};
use data::{write_str_to_file, write_vec_to_file};

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


