use manybody::manybody::{Particle, Manybody, AsArrRef};
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

    for i in 0..N_PLUS {
        ions.push(Ion {id: i, z: q_plus, point: point_initialize(&mut rng)});
    }
    for i in N_PLUS..N_PLUS+N_NEGA {
        ions.push(Ion {id: i, z: q_nega, point: point_initialize(&mut rng)});
    }

    let mut particle_system = Manybody::new(ions, rng);
    let nb: usize = 6000000;
    let ne: usize = 10000000;
    for _ in 0..nb {
        particle_system.rbmc(0.01, 99.0, 1.0/99.0, 2, 9);
    }

    let mut hist1 = Histogram::new(&[0.0], &[100.0], &[100]);
    let mut hist2 = Histogram::new(&[0.0], &[100.0], &[100]);
    for _ in nb..ne {
        let x = particle_system.rbmc(0.01, 99.0, 1.0/99.0, 2, 9);
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
