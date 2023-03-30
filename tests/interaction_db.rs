use approx::assert_relative_eq;
use manybody::{dbparticle::DBParticle, manybody::Particle};

#[test]
fn wtest() {
    let mut rng = rand::thread_rng();
    let error_margin = 1e-10;
    for _ in 0..100000 {
        let x = DBParticle { id: 0, point: DBParticle::standard_normal(&mut rng) };
        let y = DBParticle { id: 0, point: DBParticle::standard_normal(&mut rng) };
 //       println!("{}", (x.point-y.point).norm());
        assert_relative_eq!(x.w1(&y) + x.w2(&y), x.w(&y), max_relative = error_margin);
        assert_relative_eq!(x.dw1(&y) + x.dw2(&y), x.dw(&y), max_relative = error_margin);
    }
}
