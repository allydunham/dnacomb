//! Integration tests for the LibSpec module
use dnacomb::{
    Compression, ReadPairParser, SeqFormat, SeqPath, lib_spec, library::Library,
    parsing::ReadPairProducer,
};

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

#[test]
fn load_fasta() {
    let paths = vec![
        ("tests/file_loading/forward.fa", None),
        (
            "tests/file_loading/forward.fa",
            Some("tests/file_loading/reverse.fa"),
        ),
        ("tests/file_loading/forward.fa.gz", None),
        (
            "tests/file_loading/forward.fa.gz",
            Some("tests/file_loading/reverse.fa.gz"),
        ),
        ("tests/file_loading/forward.fq", None),
        (
            "tests/file_loading/forward.fq",
            Some("tests/file_loading/reverse.fq"),
        ),
        ("tests/file_loading/forward.fq.gz", None),
        (
            "tests/file_loading/forward.fq.gz",
            Some("tests/file_loading/reverse.fq.gz"),
        ),
    ];

    for (fwd, rev) in paths {
        let forward = SeqPath::new(fwd.to_string(), SeqFormat::Auto, Compression::Auto);
        let reverse = match rev {
            Some(x) => Some(SeqPath::new(
                x.to_string(),
                SeqFormat::Auto,
                Compression::Auto,
            )),
            None => None,
        };

        let reader = ReadPairParser::from_paths(forward, reverse, None, 0, b'I');

        if let Ok(reader) = reader {
            assert_eq!(
                reader.has_reverse(),
                rev.is_some(),
                "Reader reverse status doesn't match ({}, {:?})",
                fwd,
                rev
            );

            for read in reader {
                assert!(read.is_ok(), "Read failed ({}, {:?})", fwd, rev);
            }
        } else {
            panic!("Error establising reader ({}, {:?})", fwd, rev)
        }
    }
}

#[test]
fn load_library() {
    let paths = vec![
        (
            "config/grna_no_id.json",
            vec!["config/grna_no_id.tsv".to_string()],
        ),
        (
            "config/grna_sensor.json",
            vec!["config/grna_sensor.tsv".to_string()],
        ),
        ("config/grna.json", vec!["config/grna.tsv".to_string()]),
        ("config/pegrna.json", vec!["config/pegrna.tsv".to_string()]),
        (
            "config/pegrna.json",
            vec![
                "config/pegrna_spacer.tsv".to_string(),
                "config/pegrna_target.tsv".to_string(),
            ],
        ),
    ];

    for (spec, libs) in paths {
        let ls = lib_spec::LibrarySpec::from_file(spec, None, None, None, None)
            .expect("LibSpec load failed");

        let lib = Library::from_files(&libs, &ls, 1);

        assert!(!lib.is_err(), "Error in {spec}: {:?}", libs);
    }
}
