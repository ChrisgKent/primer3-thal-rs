//! The exported nearest-neighbour model ([`primer3_thal_rs::nn`]) must give
//! exactly thal's numbers, and the optional DMSO/formamide terms must shift
//! only Tm.

use primer3_thal_rs::nn::{self, Duplex};
use primer3_thal_rs::{ThalArgs, ThalError, calc_hairpin, calc_heterodimer, calc_homodimer};

fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        })
        .collect()
}

/// Deterministic random ACGT sequences of 15 to 34 nt.
fn random_seqs(n: usize, mut state: u64) -> Vec<Vec<u8>> {
    let mut next = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    (0..n)
        .map(|_| (0..15 + next() % 20).map(|_| b"ACGT"[next() % 4]).collect())
        .collect()
}

fn conditions() -> Vec<ThalArgs> {
    let d = ThalArgs::default();
    vec![
        d,
        ThalArgs {
            mv: 100.0,
            dv: 2.0,
            dntp: 0.8,
            dna_conc: 15.0,
            ..d
        },
        ThalArgs {
            mv: 40.0,
            dv: 0.0,
            dntp: 0.8,
            dna_conc: 250.0,
            ..d
        },
        ThalArgs {
            mv: 100.0,
            dv: 0.5,
            dntp: 0.8,
            dna_conc: 15.0,
            ..d
        },
        ThalArgs {
            formamide_conc: 0.8,
            ..d
        },
        ThalArgs {
            mv: 100.0,
            dv: 2.0,
            dntp: 0.8,
            dna_conc: 15.0,
            formamide_conc: 0.8,
            ..d
        },
        ThalArgs {
            dmso_conc: 5.0,
            dmso_fact: 0.6,
            formamide_conc: 1.5,
            ..d
        },
    ]
}

fn palindromes() -> Vec<Vec<u8>> {
    [
        "ACGTACGTACGTACGTACGT",
        "GGGGCCCCGGGGCCCCGGGGCCCC",
        "AATTAATTAATTAATTAATT",
        "GAATTCGAATTCGAATTC",
    ]
    .iter()
    .map(|s| s.as_bytes().to_vec())
    .collect()
}

#[test]
fn duplex_tm_equals_thal_heterodimer_with_complement() {
    let mut seqs = random_seqs(400, 11);
    seqs.extend(palindromes());
    for args in conditions() {
        for seq in &seqs {
            let thal = calc_heterodimer(seq, rc(seq), &args, false).unwrap();
            let tm = nn::duplex_tm(seq, &args).unwrap();
            assert!(
                (tm - thal.tm).abs() < 1e-9,
                "{}: {tm} vs thal {} ({args:?})",
                String::from_utf8_lossy(seq),
                thal.tm
            );
        }
    }
}

#[test]
fn bound_percent_from_thal_dimer_equals_duplex_annealing() {
    // thal's ΔH and (salt-corrected) ΔS for a primer with its complement give
    // the same binding as the duplex model, at any temperature
    let mut seqs = random_seqs(200, 13);
    seqs.extend(palindromes());
    for args in conditions() {
        for seq in &seqs {
            let thal = calc_heterodimer(seq, rc(seq), &args, false).unwrap();
            let symmetric = nn::is_self_complementary(seq);
            for t in [45.0, 60.0, 72.0] {
                let bound = nn::bound_percent(thal.dh, thal.ds, &args, symmetric, t);
                let duplex = nn::duplex_annealing(seq, &args, t).unwrap();
                assert!(
                    (bound - duplex).abs() < 1e-9,
                    "{}: {bound} vs {duplex} at {t} °C",
                    String::from_utf8_lossy(seq)
                );
            }
        }
    }
}

#[test]
fn duplex_built_base_by_base_matches_whole_sequence() {
    let args = conditions()[5];
    for seq in random_seqs(100, 3) {
        let mut d = Duplex::new();
        for &b in &seq {
            d.push(b).unwrap();
        }
        let sym = nn::is_self_complementary(&seq);
        assert_eq!(d.len(), seq.len());
        assert_eq!(d.tm(&args, sym), nn::duplex_tm(&seq, &args).unwrap());
        assert_eq!(
            d.annealing(&args, sym, 60.0),
            nn::duplex_annealing(&seq, &args, 60.0).unwrap()
        );
    }
}

#[test]
fn duplex_accepts_lowercase() {
    let args = ThalArgs::default();
    assert_eq!(
        nn::duplex_tm("acgtacgatcgatgcatgca", &args),
        nn::duplex_tm("ACGTACGATCGATGCATGCA", &args)
    );
}

#[test]
fn duplex_rejects_bad_input() {
    let args = ThalArgs::default();
    assert_eq!(nn::duplex_tm("", &args), Err(ThalError::EmptySequence));
    assert_eq!(
        nn::duplex_tm("ACGTNACGT", &args),
        Err(ThalError::InvalidBase('N'))
    );
    assert_eq!(Duplex::new().push(b'-'), Err(ThalError::InvalidBase('-')));
    // The error keeps the input's case
    assert_eq!(
        nn::duplex_tm("acgn", &args),
        Err(ThalError::InvalidBase('n'))
    );
}

#[test]
fn annealing_is_half_bound_at_tm() {
    // Without DMSO/formamide (which shift Tm but not binding), half the
    // oligo is bound at its own Tm.
    let mut seqs = random_seqs(50, 5);
    seqs.extend(palindromes());
    for args in conditions()
        .into_iter()
        .filter(|a| a.formamide_conc == 0.0 && a.dmso_conc == 0.0)
    {
        for seq in &seqs {
            let tm = nn::duplex_tm(seq, &args).unwrap();
            let bound = nn::duplex_annealing(seq, &args, tm).unwrap();
            assert!((bound - 50.0).abs() < 1e-6, "{bound}");
            // and binding falls as temperature rises
            assert!(nn::duplex_annealing(seq, &args, tm + 5.0).unwrap() < bound);
        }
    }
}

#[test]
fn tm_correction_is_primer3s_formula() {
    let args = ThalArgs {
        dmso_conc: 5.0,
        dmso_fact: 0.6,
        formamide_conc: 0.8,
        ..ThalArgs::default()
    };
    let expected = -5.0 * 0.6 + (0.453 * 0.25 - 2.88) * 0.8;
    assert!((nn::tm_correction(0.25, &args) - expected).abs() < 1e-12);
    assert_eq!(nn::tm_correction(0.25, &ThalArgs::default()), 0.0);
}

/// GC fraction of all the given sequences together.
fn gc(seqs: &[&[u8]]) -> f64 {
    let total: usize = seqs.iter().map(|s| s.len()).sum();
    let gc = seqs
        .iter()
        .flat_map(|s| s.iter())
        .filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c'))
        .count();
    gc as f64 / total as f64
}

#[test]
fn dmso_and_formamide_shift_only_tm() {
    let plain = ThalArgs {
        mv: 100.0,
        dv: 2.0,
        dntp: 0.8,
        dna_conc: 15.0,
        ..ThalArgs::default()
    };
    let corrected = ThalArgs {
        dmso_conc: 2.0,
        formamide_conc: 0.8,
        ..plain
    };
    let hairpin = b"CCCCCATCCGATCAGGGGG";
    let (a, b) = (
        b"GTAAAACGACGGCCAGTGAAT".as_slice(),
        b"ACACCTGTGCCTGTTAAACCAT".as_slice(),
    );

    let cases = [
        (
            calc_hairpin(hairpin, &plain, true).unwrap(),
            calc_hairpin(hairpin, &corrected, true).unwrap(),
            gc(&[hairpin]),
        ),
        (
            calc_homodimer(a, &plain, true).unwrap(),
            calc_homodimer(a, &corrected, true).unwrap(),
            gc(&[a]),
        ),
        (
            calc_heterodimer(a, b, &plain, true).unwrap(),
            calc_heterodimer(a, b, &corrected, true).unwrap(),
            gc(&[a, b]),
        ),
    ];
    for (p, c, gc) in cases {
        assert!(p.structure_found);
        let shift = nn::tm_correction(gc, &corrected);
        assert!(shift < -1.0);
        assert!(
            (c.tm - (p.tm + shift)).abs() < 1e-9,
            "{} vs {}",
            c.tm,
            p.tm + shift
        );
        assert_eq!((c.dg, c.dh, c.ds), (p.dg, p.dh, p.ds));
        assert_eq!(c.ascii_structure, p.ascii_structure);
    }

    // No structure: Tm stays 0
    let none = calc_hairpin("AAAAAAAAAA", &corrected, false).unwrap();
    assert!(!none.structure_found);
    assert_eq!(none.tm, 0.0);
}

/// Every ACGT sequence of 1 to 4 nt, then random ones of 5 to 14 nt.
fn short_seqs() -> Vec<Vec<u8>> {
    let mut seqs: Vec<Vec<u8>> = vec![vec![]];
    let mut all = Vec::new();
    for _ in 0..4 {
        seqs = seqs
            .iter()
            .flat_map(|s| {
                b"ACGT".iter().map(move |&b| {
                    let mut s = s.clone();
                    s.push(b);
                    s
                })
            })
            .collect();
        all.extend(seqs.iter().cloned());
    }
    let mut state = 23u64;
    for _ in 0..200 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let len = 5 + (state >> 33) as usize % 10;
        all.extend(random_seqs(1, state).into_iter().map(|s| s[..len].to_vec()));
    }
    all
}

/// Equal, or both NaN.
fn same(a: f64, b: f64, tol: f64) -> bool {
    (a.is_nan() && b.is_nan()) || (a - b).abs() < tol || a == b
}

#[test]
fn short_duplexes_match_thal() {
    // Below primer length the duplex model still gives thal's numbers,
    // however unphysical (a single base pair has a Tm far below 0 °C)
    let seqs = short_seqs();
    assert_eq!(seqs.len(), 4 + 16 + 64 + 256 + 200);
    for args in conditions() {
        for seq in &seqs {
            let thal = calc_heterodimer(seq, rc(seq), &args, false).unwrap();
            let label = String::from_utf8_lossy(seq);
            assert!(thal.structure_found, "{label}");
            let tm = nn::duplex_tm(seq, &args).unwrap();
            assert!(
                same(tm, thal.tm, 1e-9),
                "{label}: {tm} vs thal {} ({args:?})",
                thal.tm
            );
            let symmetric = nn::is_self_complementary(seq);
            let bound = nn::bound_percent(thal.dh, thal.ds, &args, symmetric, 60.0);
            let duplex = nn::duplex_annealing(seq, &args, 60.0).unwrap();
            assert!(same(bound, duplex, 1e-9), "{label}: {bound} vs {duplex}");
        }
    }
}

#[test]
fn single_bases_match_primer3() {
    // primer3-py 2.3.0: calc_heterodimer(seq, reverse_complement(seq)).tm at
    // its defaults. A single base takes the terminal A-T penalty at both ends
    let args = ThalArgs::default();
    for (seq, primer3) in [
        ("A", -437.0719),
        ("C", -277.9276),
        ("AT", -217.0828),
        ("ACGT", -43.5526),
    ] {
        let tm = nn::duplex_tm(seq, &args).unwrap();
        assert!(
            (tm - primer3).abs() < 1e-3,
            "{seq}: {tm} vs primer3 {primer3}"
        );
    }
}

#[test]
fn zero_salt_matches_thal() {
    // No cations: ln(0) in the salt correction. thal (and primer3) give NaN
    // for one base pair, where 0 · -inf appears, and -273.15 °C otherwise
    let args = ThalArgs {
        mv: 0.0,
        dv: 0.0,
        ..ThalArgs::default()
    };
    for seq in ["A", "AC", "ACGTACGTAC"] {
        let s = seq.as_bytes();
        let thal = calc_heterodimer(s, rc(s), &args, false).unwrap().tm;
        let tm = nn::duplex_tm(seq, &args).unwrap();
        assert!(same(tm, thal, 1e-9), "{seq}: {tm} vs thal {thal}");
    }
    assert!(nn::duplex_tm("A", &args).unwrap().is_nan());
    assert_eq!(nn::duplex_tm("AC", &args).unwrap(), -nn::KELVIN_OFFSET);
}

#[test]
fn bound_percent_saturates() {
    let args = conditions()[1];
    for seq in random_seqs(50, 29) {
        let d = Duplex::from_seq(&seq).unwrap();
        let mut last = 100.0;
        for t in 0..=120 {
            let bound = d.annealing(&args, false, f64::from(t));
            assert!((0.0..=100.0).contains(&bound), "{bound} at {t} °C");
            assert!(bound <= last, "not decreasing at {t} °C");
            last = bound;
        }
        assert!(d.annealing(&args, false, -50.0) > 99.999);
        assert!(d.annealing(&args, false, 200.0) < 1e-3);
    }
    // ΔH/ΔS far beyond any oligo still give a finite answer, not NaN
    for (dh, ds) in [(-1e6, -2.5e3), (1e6, 2.5e3), (-1e6, 0.0), (1e6, 0.0)] {
        for t in [-50.0, 37.0, 200.0] {
            let bound = nn::bound_percent(dh, ds, &args, false, t);
            assert!((0.0..=100.0).contains(&bound), "{dh} {ds} at {t}: {bound}");
        }
    }
}
