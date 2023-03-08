use std::fs;
use manybody::data::*;

#[test]
fn write_str() {
    write_to_file("Hello World!", "test.txt").unwrap();
    fs::remove_file("test.txt").unwrap();
}

#[test]
fn hist_test() {
    let mut dist = Histogram::new(&[0.0, 0.0], &[1.0, 1.0], &[2, 2]);
    dist.add(&[0.1, 1.2]);
    dist.add(&[0.2, 0.4]);
}
