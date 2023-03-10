pub mod manybody;
pub mod dbparticle;
pub mod data;
pub mod rbmc;



//     fn pack_rbmc<const N:usize, T>(system: &mut Manybody<T, N>, iter: usize, a: &[f64; N], b: &[f64; N], s:&[f64; N])
//     -> Histogram<N>
//     where T: Particle<N>
//     {
//         let hist = Histogram::new(a, b, s);
//         for _ in 0..iter {
//             hist.add(system.rbmc())
//         }
//     }
// }