mod common;

use primer3_py_rs::{ThalArgs, ThalError, calc_hairpin};

#[test]
fn hairpin_matches_primer3_py() {
    let cases = common::load("hairpin");
    common::check_all("hairpin", &cases, |c| {
        calc_hairpin(&c.seq, &c.args.thal(), c.args.output_structure)
    });
}

#[test]
fn hairpin_rejects_bad_input() {
    let args = ThalArgs::default();
    assert_eq!(calc_hairpin("", &args, false), Err(ThalError::EmptySequence));
    assert_eq!(
        calc_hairpin(&"A".repeat(61), &args, false),
        Err(ThalError::TooLong { len: 61, max: 60 })
    );
    for max_loop in [-1, 31] {
        let a = ThalArgs { max_loop, ..args };
        assert_eq!(
            calc_hairpin("GGGGGAAAACCCCC", &a, false),
            Err(ThalError::InvalidMaxLoop(max_loop))
        );
    }
}

#[test]
fn hairpin_is_case_insensitive() {
    let args = ThalArgs::default();
    for seq in ["CCCCCATCCGATCAGGGGG", "GGGGGAAAACCCCC", "ACGTACGTTTTTACGTACGT"] {
        let upper = calc_hairpin(seq, &args, true).unwrap();
        let lower = calc_hairpin(&seq.to_lowercase(), &args, true).unwrap();
        assert_eq!(upper.tm, lower.tm);
        assert_eq!(upper.dg, lower.dg);
        assert_eq!(upper.ascii_structure, lower.ascii_structure);
    }
}

#[cfg(feature = "parallel")]
#[test]
fn hairpin_batch_matches_single() {
    let cases = common::load("hairpin");
    let args = ThalArgs::default();
    let seqs: Vec<&str> = cases.iter().map(|c| c.seq.as_str()).collect();
    let batch = primer3_py_rs::calc_hairpin_batch(&seqs, &args, false);
    for (seq, b) in seqs.iter().zip(batch) {
        assert_eq!(calc_hairpin(seq, &args, false), b, "{seq}");
    }
}

proptest::proptest! {
    // Any input must give Ok or Err, never a panic (e.g. loop-table indexing).
    #[test]
    fn hairpin_never_panics(
        seq in "[ACGTNacgtn\\-xé ]{0,80}",
        max_loop in -5i32..40,
        temp_c in -50.0f64..150.0,
        mv in 0.0f64..1000.0,
        dv in 0.0f64..100.0,
        dntp in 0.0f64..100.0,
    ) {
        let args = ThalArgs { max_loop, mv, dv, dntp, ..ThalArgs::default() }.with_temp_c(temp_c);
        let _ = calc_hairpin(&seq, &args, true);
    }
}
