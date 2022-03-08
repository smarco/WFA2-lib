use cbindgen;

use std::env;

fn main() {
    if env::var("BLOCK_ALIGNER_C").is_ok() {
        let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
        cbindgen::generate(&crate_dir)
            .unwrap()
            .write_to_file(format!("{}/c/block_aligner.h", crate_dir));
    }
}
