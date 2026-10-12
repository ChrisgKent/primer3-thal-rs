mod common;

use primer3_thal_rs::{ThalArgs, ThalError, calc_heterodimer, calc_homodimer};

#[test]
fn homodimer_matches_primer3_py() {
    let cases = common::load("homodimer");
    common::check_all("homodimer", &cases, |c| {
        calc_homodimer(&c.seq, &c.args.thal(), c.args.output_structure)
    });
}

#[test]
fn heterodimer_matches_primer3_py() {
    let cases = common::load("heterodimer");
    common::check_all("heterodimer", &cases, |c| {
        let seq2 = c.seq2.as_deref().expect("heterodimer fixture without seq2");
        calc_heterodimer(&c.seq, seq2, &c.args.thal(), c.args.output_structure)
    });
}

#[test]
fn dimer_length_rules_follow_thal() {
    let args = ThalArgs::default();
    let short = "ACGTACGTACGTACGTACGT";
    let long = "ACGT".repeat(20);
    // One sequence over 60 bp is allowed, both is not.
    assert!(calc_heterodimer(&long, short, &args, false).is_ok());
    assert!(calc_heterodimer(short, &long, &args, false).is_ok());
    assert_eq!(
        calc_heterodimer(&long, &long, &args, false),
        Err(ThalError::TooLong { len: 80, max: 60 })
    );
    assert_eq!(
        calc_homodimer(&long, &args, false),
        Err(ThalError::TooLong { len: 80, max: 60 })
    );
    assert_eq!(
        calc_heterodimer("A".repeat(10_001), short, &args, false),
        Err(ThalError::TooLong {
            len: 10_001,
            max: 10_000
        })
    );
    assert_eq!(
        calc_heterodimer("", short, &args, false),
        Err(ThalError::EmptySequence)
    );
    assert_eq!(
        calc_homodimer("", &args, false),
        Err(ThalError::EmptySequence)
    );
}

#[cfg(feature = "parallel")]
#[test]
fn dimer_batch_matches_single() {
    let args = ThalArgs::default();
    let cases = common::load("heterodimer");
    let pairs: Vec<(String, String)> = cases
        .iter()
        .take(500)
        .map(|c| (c.seq.clone(), c.seq2.clone().unwrap()))
        .collect();
    let batch = primer3_thal_rs::calc_heterodimer_batch(&pairs, &args, false);
    for ((a, b), r) in pairs.iter().zip(batch) {
        assert_eq!(calc_heterodimer(a, b, &args, false), r, "{a}/{b}");
    }
    let seqs: Vec<&str> = pairs.iter().map(|(a, _)| a.as_str()).collect();
    let batch = primer3_thal_rs::calc_homodimer_batch(&seqs, &args, false);
    for (s, r) in seqs.iter().zip(batch) {
        assert_eq!(calc_homodimer(s, &args, false), r, "{s}");
    }
}

proptest::proptest! {
    #[test]
    fn dimers_never_panic(
        a in "[ACGTNacgtn\\-xé ]{0,80}",
        b in "[ACGTNacgtn\\-xé ]{0,80}",
        max_loop in -5i32..40,
        temp_c in -50.0f64..150.0,
        mv in 0.0f64..1000.0,
        dv in 0.0f64..100.0,
        dntp in 0.0f64..100.0,
    ) {
        let args = ThalArgs { max_loop, mv, dv, dntp, ..ThalArgs::default() }.with_temp_c(temp_c);
        let _ = calc_heterodimer(&a, &b, &args, true);
        let _ = calc_homodimer(&a, &args, true);
    }
}
