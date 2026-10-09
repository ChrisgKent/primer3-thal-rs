#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(dead_code)]

mod params;
#[cfg(feature = "python")]
mod python;
mod thal;

pub use thal::{calc_hairpin, calc_heterodimer};
#[cfg(feature = "parallel")]
pub use thal::{calc_hairpin_batch, calc_heterodimer_batch};

const MIN_HRPN_LOOP: i32 = 3;
const MAX_LOOP: i32 = 30;

const R: f64 = 1.9872; // cal/Kmol
const ILAS: f64 = -300.0 / 310.15;
const ILAH: f64 = 0.0;
const AT_H: f64 = 2200.0;
const AT_S: f64 = 6.9;
const MIN_ENTROPY_CUTOFF: f64 = -2500.0;
const MIN_ENTROPY: f64 = -3224.0;
const G2: f64 = 0.0;
const ABSOLUTE_ZERO: f64 = 273.15;
const TEMP_KELVIN: f64 = 310.15;
const SMALL_NON_ZERO: f64 = 0.000001;

const THAL_MAX_ALIGN: i32 = 60;
const THAL_MAX_SEQ: i32 = 10000;

const BPI: [[i32; 5]; 5] = [
    [0, 0, 0, 1, 0],
    [0, 0, 1, 0, 0],
    [0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
];

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
        self.temp = temp_c + ABSOLUTE_ZERO;
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

#[derive(Clone, Copy)]
struct Triloop {
    loop_vals: [u8; 5],
    value: f64,
}

impl Default for Triloop {
    fn default() -> Self {
        Self {
            loop_vals: [0; 5],
            value: 0.0,
        }
    }
}

#[derive(Clone, Copy)]
struct Tetraloop {
    loop_vals: [u8; 6],
    value: f64,
}

impl Default for Tetraloop {
    fn default() -> Self {
        Self {
            loop_vals: [0; 6],
            value: 0.0,
        }
    }
}

struct ThermoParams {
    num_triloops: usize,
    num_tetraloops: usize,
    atp_s: [[f64; 5]; 5],
    atp_h: [[f64; 5]; 5],
    dangle_entropies_3: [[[f64; 5]; 5]; 5],
    dangle_enthalpies_3: [[[f64; 5]; 5]; 5],
    dangle_entropies_5: [[[f64; 5]; 5]; 5],
    dangle_enthalpies_5: [[[f64; 5]; 5]; 5],
    stack_entropies: [[[[f64; 5]; 5]; 5]; 5],
    stack_enthalpies: [[[[f64; 5]; 5]; 5]; 5],
    stackint2_entropies: [[[[f64; 5]; 5]; 5]; 5],
    stackint2_enthalpies: [[[[f64; 5]; 5]; 5]; 5],
    interior_loop_entropies: [f64; 30],
    interior_loop_enthalpies: [f64; 30],
    bulge_loop_entropies: [f64; 30],
    bulge_loop_enthalpies: [f64; 30],
    hairpin_loop_entropies: [f64; 30],
    hairpin_loop_enthalpies: [f64; 30],
    tstack_entropies: [[[[f64; 5]; 5]; 5]; 5],
    tstack_enthalpies: [[[[f64; 5]; 5]; 5]; 5],
    tstack2_entropies: [[[[f64; 5]; 5]; 5]; 5],
    tstack2_enthalpies: [[[[f64; 5]; 5]; 5]; 5],
    triloop_entropies: [Triloop; 20],
    triloop_enthalpies: [Triloop; 20],
    tetraloop_entropies: [Tetraloop; 146],
    tetraloop_enthalpies: [Tetraloop; 146],
}

impl ThermoParams {
    fn new() -> Self {
        Self {
            num_triloops: 0,
            num_tetraloops: 0,
            atp_s: [[0.0; 5]; 5],
            atp_h: [[0.0; 5]; 5],
            dangle_entropies_3: [[[0.0; 5]; 5]; 5],
            dangle_enthalpies_3: [[[0.0; 5]; 5]; 5],
            dangle_entropies_5: [[[0.0; 5]; 5]; 5],
            dangle_enthalpies_5: [[[0.0; 5]; 5]; 5],
            stack_entropies: [[[[0.0; 5]; 5]; 5]; 5],
            stack_enthalpies: [[[[0.0; 5]; 5]; 5]; 5],
            stackint2_entropies: [[[[0.0; 5]; 5]; 5]; 5],
            stackint2_enthalpies: [[[[0.0; 5]; 5]; 5]; 5],
            interior_loop_entropies: [0.0; 30],
            interior_loop_enthalpies: [0.0; 30],
            bulge_loop_entropies: [0.0; 30],
            bulge_loop_enthalpies: [0.0; 30],
            hairpin_loop_entropies: [0.0; 30],
            hairpin_loop_enthalpies: [0.0; 30],
            tstack_entropies: [[[[0.0; 5]; 5]; 5]; 5],
            tstack_enthalpies: [[[[0.0; 5]; 5]; 5]; 5],
            tstack2_entropies: [[[[0.0; 5]; 5]; 5]; 5],
            tstack2_enthalpies: [[[[0.0; 5]; 5]; 5]; 5],
            triloop_entropies: [Triloop::default(); 20],
            triloop_enthalpies: [Triloop::default(); 20],
            tetraloop_entropies: [Tetraloop::default(); 146],
            tetraloop_enthalpies: [Tetraloop::default(); 146],
        }
    }
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
