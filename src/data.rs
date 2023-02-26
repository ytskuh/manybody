use std::fs::OpenOptions; 
use std::io::{BufWriter, Write};
use std::fmt::Display;

pub fn write_vec_to_file<T: Display> (data: &Vec<T>, file_name: &str, append: bool) 
-> std::io::Result<()> {
    let file = OpenOptions::new()
        .create(true).append(append)
        .open(file_name)?;
    let mut writer = BufWriter::new(file);

    for x in data {
        writeln!(writer, "{}", x)?;
    }
    Ok(())
}

pub fn write_str_to_file (data: &str, file_name: &str, append: bool)
-> std::io::Result<()> {
    let file = OpenOptions::new()
    .create(true).append(append)
    .open(file_name)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{}", data)?;
    Ok(())
}