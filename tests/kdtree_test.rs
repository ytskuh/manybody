use kiddo::KdTree;
use kiddo::distance::squared_euclidean;

#[test]
fn kdtree_structure() {
    let a = ([1f64, 0f64], 0);
    let b = ([2f64, 0f64], 1);

    let mut kdtree = KdTree::new();

    kdtree.add(&a.0, a.1).unwrap();
    kdtree.add(&b.0, a.1).unwrap();

    kdtree.remove(&a.0, &a.1).unwrap();
    assert_eq!(kdtree.size(), 1);
}

#[test]
fn kdtree_within() {
    let a = ([0.0, 0.0], 0);
    let b = ([1.0, 1.0], 1);
    let c = ([2.0, 0.1], 2);

    let mut kdtree = KdTree::new();
    kdtree.add(&a.0, a.1).unwrap();
    kdtree.add(&b.0, b.1).unwrap();
    kdtree.add(&c.0, c.1).unwrap();

    let result = kdtree.within(&[1.0, 0.0], 1.0, &squared_euclidean).unwrap();

    assert_eq!(result.len(), 2);
    assert!(result.contains(&(1.0, &1)));
    assert!(result.contains(&(1.0, &0)));
}
