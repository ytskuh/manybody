use clap::Parser;

use manybody::manybody::{Particle, Manybody, AsArrRef};
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
    #[arg(long)]
    step_time: f64,

    #[arg(long)]
    burn_in: Option<usize>,
    #[arg(long)]
    p: Option<usize>,
    #[arg(long)]
    m: Option<usize>,

    #[arg(long)]
    raw: bool,
    #[arg(long)]
    low: Option<f64>,
    #[arg(long)]
    high: Option<f64>,
    #[arg(long)]
    interval_num: Option<usize>
}

fn main() {
    let args = Args::parse();

    let particle_num = args.particle_num;
    let time_length = (args.iteration/particle_num as u32) as f64 * args.step_time;
    let step_time = args.step_time;
    let iteration = (particle_num as f64*time_length/step_time) as usize;
    let p = match args.p { Some(p) => p, None => 2 };
    let m = match args.m { Some(m) => m, None => 9 };
    let burn_in = match args.burn_in {
        Some(b) => b,
        None => iteration/2
    };
    let filename = &args.output;

    let mut rng = rand::thread_rng();

    let mut particle_initial = Vec::new();
    for i in 0..particle_num {
        particle_initial.push(DBParticle {
            id: i,
            point: DBParticle::standard_normal(&mut rng)
        });
    }
    let mut particle_system = Manybody::new(particle_initial, rng, step_time, p, m, 1.0, 1.0, 1.0/(particle_num as f64 - 1.0).sqrt(), 1.0);

    for _ in 0..burn_in {
        particle_system.rbmc();
    }
    if args.raw {
        write_to_file("x", filename).unwrap();
        for _ in burn_in..iteration {
            let x = particle_system.rbmc();
            append_to_file(x.point(), filename).unwrap();
        }
    } else {
        let l = args.low.unwrap();
        let h = args.high.unwrap();
        let s = args.interval_num.unwrap();
        let mut hist = Histogram::new(&[l], &[h], &[s]);
        for _ in burn_in..iteration {
            let x = particle_system.rbmc();
            hist.add(x.point().as_aref());
        }
        write_to_file(hist.hist_density(), filename).unwrap();
    }
}
