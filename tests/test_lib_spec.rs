//! Integration tests for the LibSpec module

use dnacomb::lib_spec;

#[test]
fn load_lib_spec() {
    let paths = vec![
        "config/grna_no_id.json",
        "config/grna_sensor.json",
        "config/grna.json",
        "config/pegrna.json",
    ];

    for path in paths {
        let ls = lib_spec::LibrarySpec::from_file(path, None, None, None, None);

        assert!(!ls.is_err(), "Error in {path}: {:?}", ls);
    }
}
