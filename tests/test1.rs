#[cfg(test)]
mod test {
use kdtree::KdTree;

#[test]
fn test_kdtree() {
    let a = ([0f64, 0f64], 0);
    let b = ([1f64, 0f64], 1);

    let dimensions = 2;
    let mut kdtree = KdTree::new(dimensions);

    kdtree.add(a.0, a.1).unwrap();
    kdtree.add(b.0, b.1).unwrap();

    kdtree.remove(&&a.0, &0).unwrap();
    assert_eq!(kdtree.size(), 1);
}
}
