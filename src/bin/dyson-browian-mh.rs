use clap::Parser;

use manybody::manybody::{Particle, Manybody};
use manybody::dbparticle::DBParticle;
use manybody::data::*;

#[derive(Parser, Debug)]
#[command(version)]
struct Args {
    #[arg(long, short)]
    output: String,

    #[arg(long)]
    particle_num: usize,

    #[arg(long)]
    iteration: u32,
}

fn main() {
    let args = Args::parse();
    let particle_num = args.particle_num;
    let step_num = args.iteration;

    let filename = &args.output;

    let mut rng = rand::thread_rng();

    let mut particle_initial = Vec::new();
    for i in 0..particle_num {
        particle_initial.push(DBParticle {
            id: i,
            point: DBParticle::standard_uni(&mut rng)
        });
    }
    let mut particle_system = Manybody::new(particle_initial, rng, 0.0, 0, 0, particle_num as f64 - 1.0, 1.0, 1.0, 0.0);

    write_to_file("x", filename).unwrap();
    append_vec_to_file(particle_system.particles(), filename).unwrap();

    for i in 0..step_num {
        particle_system.mh();
        if (i%particle_num as u32) == 0 {
            append_vec_to_file(particle_system.particles(), filename).unwrap();
        }
    }
}
