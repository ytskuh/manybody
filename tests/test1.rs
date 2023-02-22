#[cfg(test)]
mod test {
use kdtree::KdTree;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};

#[test]
fn test_kdtree() {
    let a = ([1f64, 0f64], 0);
    let b = ([2f64, 0f64], 1);

    let dimensions = 2;
    let mut kdtree = KdTree::new(dimensions);

    kdtree.add(a.0, a.1).unwrap();
    kdtree.add(b.0, a.1).unwrap();

    kdtree.remove(&[1f64, 0f64], &0).unwrap();
    assert_eq!(kdtree.size(), 1);
}

#[test]
fn lab() {
    let a = Some(3);
    let b = a;
    assert_eq!(b.unwrap(), 3);
}

#[test]
fn file_write() {
    let file = OpenOptions::new()
    .create(true)
    .append(true)
    .open("my_file.txt")
    .unwrap();
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{}", "123 123").unwrap();
    
}

}
