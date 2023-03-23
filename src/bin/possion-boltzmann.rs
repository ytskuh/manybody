use manybody::manybody::{Particle, Manybody, AsArrRef};
use manybody::ion::{Ion, DIM, N, Q};
use manybody::data::*;

fn point_initialize(rng: &mut rand::rngs::ThreadRng)
-> <Ion as Particle<DIM>>::Point {
    let x = Ion::standard_normal(rng);
    x*(1.0+1.0/x.norm())
}

fn main() {
    let mut ions = Vec::new();
    let mut rng = rand::thread_rng();
    for i in 0..N {
        ions.push(Ion {id: 3*i, z: Q, point: point_initialize(&mut rng) });
        ions.push(Ion {id: 3*i+1, z: -Q, point: point_initialize(&mut rng)});
        ions.push(Ion {id: 3*i+2, z: -Q, point: point_initialize(&mut rng)});
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
