use std::fs;
use std::path::PathBuf;

use primer3_py_rs::{ThalArgs, ThalError, ThermoResult};
use serde::{Deserialize, Deserializer};

/// Absolute tolerance for tm, dg, dh, ds against primer3-py.
pub const TOL: f64 = 1e-6;

#[derive(Debug, Deserialize)]
pub struct Args {
    pub mv_conc: f64,
    pub dv_conc: f64,
    pub dntp_conc: f64,
    pub dna_conc: f64,
    pub temp_c: f64,
    pub max_loop: i32,
    pub output_structure: bool,
}

impl Args {
    pub fn thal(&self) -> ThalArgs {
        ThalArgs {
            max_loop: self.max_loop,
            mv: self.mv_conc,
            dv: self.dv_conc,
            dntp: self.dntp_conc,
            dna_conc: self.dna_conc,
            temp: 0.0,
        }
        .with_temp_c(self.temp_c)
    }
}

/// Fixtures store non-finite values as the strings "inf", "-inf", "nan".
fn num<'de, D: Deserializer<'de>>(d: D) -> Result<f64, D::Error> {
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum Num {
        F(f64),
        S(String),
    }
    match Num::deserialize(d)? {
        Num::F(v) => Ok(v),
        Num::S(s) => s.parse().map_err(serde::de::Error::custom),
    }
}

#[derive(Debug, Deserialize)]
pub struct Expected {
    pub structure_found: bool,
    #[serde(deserialize_with = "num")]
    pub tm: f64,
    #[serde(deserialize_with = "num")]
    pub dg: f64,
    #[serde(deserialize_with = "num")]
    pub dh: f64,
    #[serde(deserialize_with = "num")]
    pub ds: f64,
    pub ascii_structure: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct Case {
    pub seq: String,
    pub seq2: Option<String>,
    pub args: Args,
    pub result: Option<Expected>,
    pub error: Option<String>,
}

pub fn load(name: &str) -> Vec<Case> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(format!("{name}.jsonl"));
    let text = fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
    // First line is a header with the primer3-py version.
    let cases: Vec<Case> = text
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).expect("parse fixture"))
        .collect();
    assert!(!cases.is_empty(), "{name}.jsonl had no cases");
    cases
}

/// Compare every case and panic once with a summary, so a regression shows its
/// full extent instead of the first mismatch.
pub fn check_all(name: &str, cases: &[Case], f: impl Fn(&Case) -> Result<ThermoResult, ThalError>) {
    let mut failures = Vec::new();
    let mut worst = [0.0f64; 4];
    for (n, case) in cases.iter().enumerate() {
        let got = f(case);
        let label = match &case.seq2 {
            Some(s2) => format!("#{n} {:?}/{:?} {:?}", case.seq, s2, case.args),
            None => format!("#{n} {:?} {:?}", case.seq, case.args),
        };
        match (&case.result, &case.error, got) {
            (Some(exp), _, Ok(got)) => {
                if exp.structure_found != got.structure_found {
                    failures.push(format!(
                        "{label}: structure_found expected {} got {}",
                        exp.structure_found, got.structure_found
                    ));
                    continue;
                }
                let pairs = [
                    ("tm", exp.tm, got.tm),
                    ("dg", exp.dg, got.dg),
                    ("dh", exp.dh, got.dh),
                    ("ds", exp.ds, got.ds),
                ];
                for (k, (field, e, g)) in pairs.into_iter().enumerate() {
                    let d = (e - g).abs();
                    let same = d <= TOL || e == g || (e.is_nan() && g.is_nan());
                    if !same {
                        failures.push(format!("{label}: {field} expected {e} got {g}"));
                    } else if d.is_finite() {
                        worst[k] = worst[k].max(d);
                    }
                }
                // primer3-py reports "" where no structure was found; Rust uses None.
                if let Some(exp_s) = &exp.ascii_structure {
                    if got.ascii_structure.as_deref().unwrap_or("") != exp_s {
                        failures.push(format!(
                            "{label}: structure expected {exp_s:?} got {:?}",
                            got.ascii_structure
                        ));
                    }
                }
            }
            (None, Some(err), Err(_)) => {
                let _ = err;
            }
            (Some(_), _, Err(e)) => failures.push(format!("{label}: unexpected error {e}")),
            (None, Some(err), Ok(got)) => {
                failures.push(format!("{label}: expected {err}, got {got:?}"))
            }
            (None, None, _) => panic!("{label}: fixture has neither result nor error"),
        }
    }
    eprintln!(
        "{name}: {} cases, worst |diff| tm={:.1e} dg={:.1e} dh={:.1e} ds={:.1e}",
        cases.len(),
        worst[0],
        worst[1],
        worst[2],
        worst[3]
    );
    if !failures.is_empty() {
        let shown: Vec<_> = failures.iter().take(25).collect();
        panic!(
            "{name}: {} of {} cases failed:\n{}",
            failures.len(),
            cases.len(),
            shown.iter().map(|s| s.as_str()).collect::<Vec<_>>().join("\n")
        );
    }
}
