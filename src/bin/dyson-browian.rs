use clap::Parser;

use manybody::manybody::{Particle, Manybody};
use manybody::dbparticle::DBParticle;
use manybody::data::*;
use manybody::rbmc::AutoRBMC;

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
    step_time: f64,

    #[arg(long)]
    distribution: bool,
    #[arg(long)]
    low: Option<f64>,
    #[arg(long)]
    high: Option<f64>,
    #[arg(long)]
    pack_step: Option<usize>
}

fn main() {
    let args = Args::parse();

    let particle_num = args.particle_num;
    let time_length = (args.iterations/particle_num as u32) as f64 * args.step_time;
    let step_time = args.step_time;
    let step_num = (particle_num as f64*time_length/step_time) as usize;
    let p = 2;
    let filename = &args.output;

    let low;
    let high;
    let interval_num;

    let mut rng = rand::thread_rng();

    let mut particle_initial = Vec::new();
    for i in 0..particle_num {
        particle_initial.push(DBParticle {
            id: i,
            point: DBParticle::standard_normal(&mut rng)
        });
    }
    let mut particle_system = Manybody::new(particle_initial, rng);

    let dt = step_time;
    let beta = (particle_num as f64 - 1.0).powi(2);
    let omega = 1.0/(particle_num as f64 - 1.0);

    if args.distribution {
        low = args.low.unwrap();
        high = args.high.unwrap();
        let pack_step = args.pack_step.unwrap();
        interval_num = step_num / pack_step;
        let distribution = Histogram::new(&[low], &[high], &[interval_num]);

        let mut lab = AutoRBMC::new(particle_system, distribution, p, dt, beta, omega, pack_step, 0.9);

        let hist = lab.sample(interval_num);
        write_to_file(hist, filename).unwrap();
    } else {
        write_to_file("x", filename).unwrap();
        append_vec_to_file(particle_system.particles(), filename).unwrap();
        for i in 0..step_num {
            particle_system.rbmc(dt, beta, omega, p);
            if (i%particle_num) == 0 {
                append_vec_to_file(particle_system.particles(), filename).unwrap();
            }
        }
    }
}

