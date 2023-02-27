use clap::Parser;

use manybody::manybody::{Particle, Manybody};
use manybody::dbparticle::{DBParticle};
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

    #[arg(long)]
    step_time: f64
}

fn main() {
    let args = Args::parse();

    let particle_num = args.particle_num;
    let time_length = args.iterations as f64 * args.step_time;
    let step_time = args.step_time;
    let step_num = (particle_num as f64*time_length/step_time) as u32;
    let beta = 2f64;
    let p = 2;
    let filename = &args.output;
    println!("{}", args.output);

    let mut rng = rand::thread_rng();

    let mut particle_initial = Vec::new();
    for i in 0..particle_num {
        particle_initial.push(DBParticle {
            id: i,
            point: DBParticle::standard_normal(&mut rng)
        });
    }
    let mut particle_system = Manybody::new(particle_initial);

    write_str_to_file("x", filename).unwrap();

    for i in 0..step_num {
        particle_system.rbmc(step_time, beta* (particle_num-1) as f64, p, &mut rng);
        if (i%(particle_num*3) as u32) == 0 {
            append_vec_to_file(particle_system.particles(), filename).unwrap();
        }
    }
}


