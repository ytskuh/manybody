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
    step_time: f64,

    #[arg(long)]
    distribution: bool,
    #[arg(long)]
    low: Option<f64>,
    #[arg(long)]
    high: Option<f64>,
    #[arg(long)]
    interval_num: Option<usize>,
    #[arg(long)]
    burnin_step: Option<usize>
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
    let mut distribution;

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

    let dt = step_time;
    let beta = (particle_num as f64 - 1.0).powf(2.0);
    let omega = 1.0/(particle_num as f64 - 1.0);

    if args.distribution {
        low = args.low.unwrap();
        high = args.high.unwrap();
        interval_num = args.interval_num.unwrap();
        let burnin_step = args.burnin_step.unwrap();
        distribution = vec![0_usize; interval_num];
        for _ in 0..burnin_step {
            particle_system.rbmc(dt, beta, omega, p, &mut rng)[0];
        }
        for _ in burnin_step..step_num {
            let x = particle_system.rbmc(dt, beta, omega, p, &mut rng)[0];
            let index = (x-low)/(high-low);
            if index >= 0.0 && index < 1.0 {
                distribution[(index * interval_num as f64) as usize]+=1;
            }
        }
        append_vec_to_file(&distribution, filename).unwrap();
    } else {
        append_vec_to_file(particle_system.particles(), filename).unwrap();
        for i in 0..step_num {
            particle_system.rbmc(dt, beta, omega, p, &mut rng);
            if (i%particle_num) == 0 {
                append_vec_to_file(particle_system.particles(), filename).unwrap();
            }
        }
    }
}
