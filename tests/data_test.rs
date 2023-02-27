use std::fs;
use manybody::data::write_str_to_file;

#[test]
fn write_str() {
    write_str_to_file("Hello World!", "test.txt").unwrap();
    fs::remove_file("test.txt").unwrap();
}
