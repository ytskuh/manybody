use clap::Parser;

use manybody::manybody::{Particle, Manybody};
use manybody::dbparticle::DBParticle;
use manybody::data::{write_str_to_file, append_vec_to_file};

#[derive(Parser, Debug)]
#[command(version)]
struct Args {
    #[arg(long, short)]
    output: String,

    #[arg(long)]
    particle_num: usize,

    #[arg(long)]
    iterations: u32,
}

fn main() {
    let args = Args::parse();
    let particle_num = args.particle_num;
    let step_num = args.iterations;
    let beta = 1.0;
    let omega = 1.0/(particle_num as f64 - 1.0);
    let filename = &args.output;

    let mut rng = rand::thread_rng();

    let mut particle_initial = Vec::new();
    for i in 0..particle_num {
        particle_initial.push(DBParticle {
            id: i,
            point: DBParticle::standard_uni(&mut rng)
        });
    }
    let mut particle_system = Manybody::new(particle_initial);

    write_str_to_file("x", filename).unwrap();
    append_vec_to_file(particle_system.particles(), filename).unwrap();

    for i in 0..step_num {
        particle_system.mh(beta, omega, &mut rng);
        if (i%particle_num as u32) == 0 {
            append_vec_to_file(particle_system.particles(), filename).unwrap();
        }
    }
}