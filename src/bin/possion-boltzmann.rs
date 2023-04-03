use manybody::manybody::{Particle, Manybody};
use manybody::ion::{Ion, DIM, N_NEGA, N_PLUS, Q_PLUS, QF};
use manybody::data::*;

fn point_initialize(rng: &mut rand::rngs::ThreadRng)
-> <Ion as Particle<DIM>>::Point {
    let x = Ion::standard_normal(rng);
    x*(1.0+1.0/x.norm())
}

fn main() {
    let mut ions = Vec::new();
    let mut rng = rand::thread_rng();
    let q_plus = Q_PLUS/N_PLUS as f64;
    let q_nega = -(QF+Q_PLUS)/N_NEGA as f64;

    let particle_num = N_PLUS + N_NEGA;
    let dt = 0.01;
    let p = 2;
    let m = 9;

    for i in 0..N_PLUS {
        ions.push(Ion {id: i, z: q_plus, point: point_initialize(&mut rng)});
    }
    for i in N_PLUS..N_PLUS+N_NEGA {
        ions.push(Ion {id: i, z: q_nega, point: point_initialize(&mut rng)});
    }

    let mut particle_system = Manybody::new(ions, rng, dt, p, m, 1.0, particle_num as f64 - 1.0, 1.0, 1.0);
    let nb: usize = 6000000;
    let ne: usize = 10000000;
    for _ in 0..nb {
        particle_system.rbmc();
    }

    let histrange = (0.0, 10.0);
    let bins_num = 30;

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
