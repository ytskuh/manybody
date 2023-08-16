use kiddo::KdTree;
use kiddo::distance::squared_euclidean;

#[test]
fn kdtree_structure() {
    let a = ([1f64, 0f64], 0);
    let b = ([2f64, 0f64], 1);

    let mut kdtree = KdTree::new();

    kdtree.add(&a.0, a.1);
    kdtree.add(&b.0, a.1);

    kdtree.remove(&a.0, a.1);
    assert_eq!(kdtree.size(), 1);
}

#[test]
fn kdtree_within() {
    let a = ([0.0, 0.0], 0);
    let b = ([1.0, 1.0], 1);
    let c = ([2.0, 0.1], 2);

    let mut kdtree = KdTree::new();
    kdtree.add(&a.0, a.1);
    kdtree.add(&b.0, b.1);
    kdtree.add(&c.0, c.1);

    let result = kdtree.within(&[1.0, 0.0], 1.01, &squared_euclidean).iter().map(|neighbour| neighbour.item).collect::<Vec<_>>();

    assert_eq!(result.len(), 2);
    assert!(result.contains(&1));
    assert!(result.contains(&0));

    let result = kdtree.within(&[1.0, 0.0], 0.1, &squared_euclidean).iter().map(|neighbour| neighbour.item).collect::<Vec<_>>();
    assert_eq!(result.len(), 0);
}
