#[cfg(test)]
mod test {
use kiddo::KdTree;

#[test]
fn test_kdtree() {
    let a = ([1f64, 0f64], 0);
    let b = ([2f64, 0f64], 1);

    let mut kdtree = KdTree::new();

    kdtree.add(&a.0, a.1).unwrap();
    kdtree.add(&b.0, a.1).unwrap();

    kdtree.remove(&[1f64, 0f64], &0).unwrap();
    assert_eq!(kdtree.size(), 1);
}
}

