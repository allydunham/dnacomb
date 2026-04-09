//! Integration tests for the LibSpec module

use dnacomb::lib_spec;

#[test]
fn load_lib_spec() {
    let path = "config/grna_sensor.json";

    let ls = lib_spec::LibrarySpec::from_file(path, None, None, None, None);

    assert!(!ls.is_err());
}
