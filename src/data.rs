use std::fs::{File, OpenOptions}; 
use std::io::{BufWriter, Write};
use std::fmt::Display;

pub fn append_vec_to_file<T: Display> (data: &Vec<T>, file_name: &str) 
-> std::io::Result<()> {
    let file = OpenOptions::new()
        .create(false).append(true)
        .open(file_name)?;
    let mut writer = BufWriter::new(file);

    for x in data {
        writeln!(writer, "{}", x)?;
    }
    Ok(())
}

pub fn write_str_to_file (data: &str, file_name: &str)
-> std::io::Result<()> {
    let file = File::create(file_name)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{}", data)?;
    Ok(())
}
