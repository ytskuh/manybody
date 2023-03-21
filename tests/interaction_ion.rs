use approx::assert_relative_eq;
use manybody::{ion::Ion, manybody::Particle};

#[test]
fn wtest() {
    let mut rng = rand::thread_rng();
    let error_margin = 1e-10;
    for _ in 0..100000 {
        let x = Ion { id: 0, z: 1.0, point: Ion::standard_normal(&mut rng) };
        let y = Ion { id: 0, z: 2.0, point: Ion::standard_normal(&mut rng) };
        assert_relative_eq!(x.w1(&y) + x.w2(&y), x.w(&y), max_relative = error_margin);
        assert_relative_eq!(x.dw1(&y) + x.dw2(&y), x.dw(&y), max_relative = error_margin);
    }
}
