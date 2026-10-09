use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use primer3_py_rs::{calc_heterodimer, ThalArgs};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct TortureArgs {
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
}

#[derive(Debug, Deserialize)]
struct TortureResult {
    structure_found: bool,
    tm: f64,
    ds: f64,
    dh: f64,
    dg: f64,
}

#[derive(Debug, Deserialize)]
struct HeterodimerRecord {
    func: String,
    seq1: String,
    seq2: String,
    args: TortureArgs,
    result: TortureResult,
}

#[test]
fn heterodimer_torture_data_matches_primer3() {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("torture_heterodimer.jsonl");
    let file = File::open(&path)
        .unwrap_or_else(|e| panic!("failed to open {}: {e}", path.display()));
    let reader = BufReader::new(file);

    let tm_tol = 0.05;
    let ds_tol = 0.2;
    let dh_tol = 1.0;
    let dg_tol = 0.5;

    let mut count = 0usize;
    for line in reader.lines() {
        let line = line.expect("read line");
        if line.trim().is_empty() {
            continue;
        }
        let rec: HeterodimerRecord = serde_json::from_str(&line).expect("parse heterodimer record");
        assert_eq!(rec.func, "heterodimer", "unexpected func in heterodimer data");

        let args = ThalArgs {
            max_loop: rec.args.max_loop,
            mv: rec.args.mv_conc,
            dv: rec.args.dv_conc,
            dntp: rec.args.dntp_conc,
            dna_conc: rec.args.dna_conc,
            temp: rec.args.temp_c + 273.15,
        };

        let res = calc_heterodimer(&rec.seq1, &rec.seq2, &args, rec.args.output_structure).unwrap();

        assert_eq!(
            res.structure_found, rec.result.structure_found,
            "structure_found mismatch"
        );
        assert!(
            (res.tm - rec.result.tm).abs() <= tm_tol,
            "tm mismatch: expected {}, got {}",
            rec.result.tm,
            res.tm
        );
        assert!(
            (res.ds - rec.result.ds).abs() <= ds_tol,
            "ds mismatch: expected {}, got {}",
            rec.result.ds,
            res.ds
        );
        assert!(
            (res.dh - rec.result.dh).abs() <= dh_tol,
            "dh mismatch: expected {}, got {}",
            rec.result.dh,
            res.dh
        );
        assert!(
            (res.dg - rec.result.dg).abs() <= dg_tol,
            "dg mismatch: expected {}, got {}",
            rec.result.dg,
            res.dg
        );

        count += 1;
    }

    assert!(count > 0, "torture_heterodimer.jsonl had no records");
}
