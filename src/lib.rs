//! Rust port of primer3's `thal` thermodynamic alignment: hairpin,
//! homodimer and heterodimer ΔG/Tm, matching primer3-py's results.

mod params;
#[cfg(feature = "python")]
mod python;
mod thal;
mod thermo;

pub use thal::{calc_hairpin, calc_heterodimer, calc_homodimer};
#[cfg(feature = "parallel")]
pub use thal::{calc_hairpin_batch, calc_heterodimer_batch, calc_homodimer_batch};

/// Smallest hairpin loop thal considers.
const MIN_HRPN_LOOP: usize = 3;
/// Largest bulge/internal loop; the loop tables have this many entries.
const MAX_LOOP: i32 = 30;

const R: f64 = 1.9872; // cal/K/mol
/// Internal-loop asymmetry penalty per unpaired base.
const ILAS: f64 = -300.0 / 310.15;
const ILAH: f64 = 0.0;
/// Terminal AT penalty.
const AT_H: f64 = 2200.0;
const AT_S: f64 = 6.9;
const MIN_ENTROPY_CUTOFF: f64 = -2500.0;
const MIN_ENTROPY: f64 = -3224.0;
const KELVIN_OFFSET: f64 = 273.15;
/// thal ranks candidate structures by ΔG at a fixed 37 °C.
const TEMP_KELVIN: f64 = 310.15;
const SMALL_NON_ZERO: f64 = 0.000001;
/// thal's stand-in for "zero" entropy, avoiding 0/0 in Tm ratios.
const TINY_ENTROPY: f64 = 0.00000000001;
/// Tolerance when matching DP values during traceback.
const TRACE_EPS: f64 = 1e-5;
/// Duplex initiation for dimers.
const DPLX_INIT_H: f64 = 200.0;
const DPLX_INIT_S: f64 = -5.7;
/// Strand-concentration divisors (nM) for self-complementary and other duplexes.
const SYMMETRIC_DNA_DIVISOR: f64 = 1_000_000_000.0;
const DNA_DIVISOR: f64 = 4_000_000_000.0;

/// At least one sequence of a pair must be no longer than this.
const THAL_MAX_ALIGN: usize = 60;
const THAL_MAX_SEQ: usize = 10000;

/// Conditions for a thermodynamic alignment. Concentrations are in mM
/// (`mv`, `dv`, `dntp`) and nM (`dna_conc`); `temp` is in Kelvin.
#[derive(Clone, Copy, Debug)]
pub struct ThalArgs {
    pub max_loop: i32,
    pub mv: f64,
    pub dv: f64,
    pub dntp: f64,
    pub dna_conc: f64,
    pub temp: f64, // Kelvin
}

/// primer3-py defaults (`primer3.argdefaults.Primer3PyArguments`).
impl Default for ThalArgs {
    fn default() -> Self {
        Self {
            max_loop: MAX_LOOP,
            mv: 50.0,
            dv: 1.5,
            dntp: 0.6,
            dna_conc: 50.0,
            temp: TEMP_KELVIN,
        }
    }
}

impl ThalArgs {
    pub fn with_temp_c(mut self, temp_c: f64) -> Self {
        self.temp = temp_c + KELVIN_OFFSET;
        self
    }

    fn validate(&self) -> Result<(), ThalError> {
        if !(0..=MAX_LOOP).contains(&self.max_loop) {
            return Err(ThalError::InvalidMaxLoop(self.max_loop));
        }
        Ok(())
    }
}

/// Inputs primer3 rejects. Mirrors the errors raised by primer3-py.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ThalError {
    EmptySequence,
    TooLong { len: usize, max: usize },
    InvalidMaxLoop(i32),
}

impl std::fmt::Display for ThalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ThalError::EmptySequence => write!(f, "empty sequence"),
            ThalError::TooLong { len, max } => {
                write!(f, "sequence length {len} exceeds maximum of {max}bp")
            }
            ThalError::InvalidMaxLoop(v) => {
                write!(f, "max_loop must be between 0 and {MAX_LOOP}, received {v}")
            }
        }
    }
}

impl std::error::Error for ThalError {}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ThermoResult {
    pub structure_found: bool,
    pub tm: f64,
    pub dg: f64,
    pub dh: f64,
    pub ds: f64,
    pub ascii_structure: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hairpin_matches_primer3_py_defaults() {
        let res = calc_hairpin("CCCCCATCCGATCAGGGGG", &ThalArgs::default(), false).unwrap();
        assert!(res.structure_found);
        assert!((res.tm - 66.74727264598937).abs() < 1e-9);
        assert!((res.dg - (-3176.918686764744)).abs() < 1e-6);
        assert_eq!(res.dh, -36300.0);
    }
}
