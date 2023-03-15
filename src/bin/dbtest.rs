use manybody::manybody::{Particle, Manybody, AsArrRef};
use manybody::dbparticle::DBParticle;
use manybody::data::*;

fn main() {
    let particle_num = 640;
    let step_time = 0.00002;
    let p = 2;

    let low = -2.0;
    let high = 2.0;
    let interval_num = 100_usize;

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

    for n in 0_u32..13 {
        let mut hist = Histogram::new(&[low], &[high], &[interval_num]);
        for _ in 0..100000*2_u32.pow(n) {
            hist.add(particle_system.rbmc(dt, beta, omega, p, 1).as_aref());
        }
        write_to_file(hist.hist_density(), &n.to_string()).unwrap();
    }
}
