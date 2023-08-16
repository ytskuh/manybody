use manybody::manybody::{Particle, Manybody};
use manybody::ion::{Ion, DIM, QF};
use manybody::data::*;

use rand::{SeedableRng, Rng};

use clap::Parser;

fn point_initialize<R: Rng>(rng: &mut R)
-> <Ion as Particle<DIM>>::Point {
    let x = Ion::standard_normal(rng);
    x*(1.0+1.0/x.norm())
}

#[derive(Parser, Debug)]
#[command(version)]
struct Args {
    #[arg(long)]
    low: f64,
    #[arg(long)]
    high: f64,
    #[arg(long)]
    interval_num: usize
}

fn main() {
    let mut ions = Vec::new();
    let seed: [u8; 32] = [0; 32];
    let mut rng = rand::rngs::StdRng::from_seed(seed);
    static N_PLUS: usize = 16;
    static N_NEGA: usize = 17;
    static Q_PLUS: f64 = 10.0;

    let q_plus = Q_PLUS/N_PLUS as f64;
    let q_nega = -(QF+Q_PLUS)/N_NEGA as f64;
    let q = q_plus;

    let particle_num = N_PLUS + N_NEGA;
    let dt = 0.0001;
    let p = 2;
    let m = 9;
    for i in 0..N_PLUS {
        ions.push(Ion {id: i, z: q_plus, point: point_initialize(&mut rng)});
    }
    for i in N_PLUS..N_PLUS+N_NEGA {
        ions.push(Ion {id: i, z: q_nega, point: point_initialize(&mut rng)});
    }

    let mut particle_system = Manybody::new(ions, rng, dt, p, m, 1.0/q, (particle_num as f64 - 1.0)/q, 1.0, 1.0/q);
    let nb: usize = 160000000;
    let ne: usize = 200000000;
    for _ in 0..nb {
        particle_system.rbmc();
 //       if i%10000 == 0 {println!("{}",i)}
    }

    let args = Args::parse();

    let histrange = (args.low, args.high);
    let bins_num = args.interval_num;

    let mut hist1 = Histogram::new(&[histrange.0], &[histrange.1], &[bins_num]);
    let mut hist2 = Histogram::new(&[histrange.0], &[histrange.1], &[bins_num]);
    for _ in nb..ne {
        let x = particle_system.rbmc();
        if x.z > 0.0 {
            hist1.add(&[x.point.norm()]);
        } else {
            hist2.add(&[x.point.norm()]);
        }
    }
    println!("{}", hist1.hist_density());
    println!("{}", hist2.hist_density());
    println!("{}", particle_system.count);
}
