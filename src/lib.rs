//! A pure-Rust port of primer3's `thal` thermodynamic alignment: hairpin,
//! homodimer and heterodimer ΔG, ΔH, ΔS and Tm for DNA oligos.
//!
//! **This crate is an LLM transpilation of
//! [primer3-py v2.3.0](https://github.com/libnano/primer3-py/tree/v2.3.0)**
//! (its bundled `thal.c`, parameter tables and `thermoanalysis.pyx`
//! defaults), restructured and optimised in Rust. It is not an official
//! primer3 or primer3-py project, and it reproduces that version only.
//!
//! Results match [primer3-py](https://github.com/libnano/primer3-py)'s
//! `calc_hairpin`, `calc_homodimer` and `calc_heterodimer` (primer3 2.x
//! `thal.c`): the test suite checks about 18,000 generated cases against
//! primer3-py 2.3.0, and all agree to within 1e-10 (test tolerance 1e-6).
//! The nearest-neighbour parameter tables are compiled in, so there is no
//! runtime configuration or file access.
//!
//! ```
//! use primer3_thal_rs::{ThalArgs, calc_hairpin, calc_heterodimer};
//!
//! // primer3-py's default conditions; temperatures are set in °C.
//! let args = ThalArgs::default().with_temp_c(37.0);
//!
//! let hairpin = calc_hairpin("CCCCCATCCGATCAGGGGG", &args, false)?;
//! assert!(hairpin.structure_found);
//!
//! // Sequences can be &str, String, &[u8] or Vec<u8>, given 5'->3'.
//! let dimer = calc_heterodimer(b"GTAAAACGACGGCCAGT", "ACTGGCCGTCGTTTTAC", &args, true)?;
//! println!("Tm {:.1} °C, ΔG {:.0} cal/mol", dimer.tm, dimer.dg);
//! println!("{}", dimer.ascii_structure.unwrap());
//! # Ok::<(), primer3_thal_rs::ThalError>(())
//! ```
//!
//! The [`nn`] module exposes thal's nearest-neighbour model for a perfectly
//! matched duplex (Tm and annealing on exactly thal's assumptions, without an
//! alignment). `ThalArgs` can also carry DMSO and formamide, which correct Tm
//! as primer3's oligotm does. Both are additions beyond primer3, and the
//! defaults leave results unchanged.
//!
//! # Features
//!
//! - `parallel` (default): `calc_*_batch` functions that run many
//!   calculations in parallel with rayon.
//! - `python`: PyO3 bindings, built with maturin as the `primer3_thal_rs`
//!   Python module. Not needed when using the crate from Rust.
//!
//! # Licence
//!
//! GPL-2.0-or-later, as primer3's `thal.c` and parameter tables, from which
//! this is derived.

#![warn(missing_docs)]

pub mod nn;
mod params;
#[cfg(feature = "python")]
mod python;
mod thal;
mod thermo;

use nn::{AT_H, AT_S, DPLX_INIT_H, DPLX_INIT_S, KELVIN_OFFSET};
pub use thal::{calc_hairpin, calc_heterodimer, calc_homodimer};
#[cfg(feature = "parallel")]
pub use thal::{calc_hairpin_batch, calc_heterodimer_batch, calc_homodimer_batch};

/// Smallest hairpin loop thal considers.
const MIN_HRPN_LOOP: usize = 3;
/// Largest bulge/internal loop; the loop tables have this many entries.
const MAX_LOOP: i32 = 30;

/// Internal-loop asymmetry penalty per unpaired base.
const ILAS: f64 = -300.0 / 310.15;
const ILAH: f64 = 0.0;
const MIN_ENTROPY_CUTOFF: f64 = -2500.0;
const MIN_ENTROPY: f64 = -3224.0;
/// thal ranks candidate structures by ΔG at a fixed 37 °C.
const TEMP_KELVIN: f64 = 310.15;
const SMALL_NON_ZERO: f64 = 0.000001;
/// thal's stand-in for "zero" entropy, avoiding 0/0 in Tm ratios.
const TINY_ENTROPY: f64 = 0.00000000001;
/// Tolerance when matching DP values during traceback.
const TRACE_EPS: f64 = 1e-5;

/// At least one sequence of a pair must be no longer than this.
const THAL_MAX_ALIGN: usize = 60;
const THAL_MAX_SEQ: usize = 10000;

/// Reaction conditions for a calculation.
///
/// Build once and reuse: it is `Copy`, and every `calc_*` function takes it
/// by reference. Start from [`ThalArgs::default`] (primer3-py's defaults)
/// and override fields as needed.
///
/// ```
/// use primer3_thal_rs::ThalArgs;
///
/// let args = ThalArgs { mv: 100.0, dv: 2.0, dntp: 0.8, ..ThalArgs::default() }.with_temp_c(60.0);
/// assert_eq!(args.temp, 333.15);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct ThalArgs {
    /// Largest bulge/internal loop considered, 0 to 30 (default 30).
    pub max_loop: i32,
    /// Monovalent cation concentration, mM (default 50).
    pub mv: f64,
    /// Divalent cation concentration, mM (default 1.5).
    pub dv: f64,
    /// dNTP concentration, mM (default 0.6).
    pub dntp: f64,
    /// DNA (oligo) concentration, nM (default 50).
    pub dna_conc: f64,
    /// Temperature for the reported ΔG, in **Kelvin** (default 310.15,
    /// i.e. 37 °C). Set it in °C with [`ThalArgs::with_temp_c`].
    pub temp: f64,
    /// DMSO concentration, % (default 0). Lowers Tm by
    /// `dmso_conc · dmso_fact`. Not in primer3's thal: an extension using
    /// primer3 oligotm's correction ([`nn::tm_correction`]).
    pub dmso_conc: f64,
    /// Tm reduction per % DMSO, °C (default 0.6, primer3's).
    pub dmso_fact: f64,
    /// Formamide concentration, mol/l (default 0). Shifts Tm by
    /// `(0.453 · GC fraction - 2.88) · formamide_conc`. Not in primer3's
    /// thal: an extension using primer3 oligotm's correction.
    pub formamide_conc: f64,
}

/// primer3-py's defaults (`primer3.argdefaults`): 50 mM monovalent,
/// 1.5 mM divalent, 0.6 mM dNTP, 50 nM DNA, 37 °C, `max_loop` 30, and no
/// DMSO or formamide.
impl Default for ThalArgs {
    fn default() -> Self {
        Self {
            max_loop: MAX_LOOP,
            mv: 50.0,
            dv: 1.5,
            dntp: 0.6,
            dna_conc: 50.0,
            temp: TEMP_KELVIN,
            dmso_conc: 0.0,
            dmso_fact: 0.6,
            formamide_conc: 0.0,
        }
    }
}

impl ThalArgs {
    /// These conditions at `temp_c` degrees Celsius.
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

/// Inputs primer3 rejects; primer3-py raises the equivalent errors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ThalError {
    /// A sequence was empty.
    EmptySequence,
    /// A sequence is longer than thal allows: 60 bases for hairpins and
    /// homodimers, 60 for the shorter sequence of a heterodimer, and 10,000
    /// for any sequence.
    TooLong {
        /// Length of the offending sequence.
        len: usize,
        /// The limit it exceeds.
        max: usize,
    },
    /// [`ThalArgs::max_loop`] was outside 0 to 30.
    InvalidMaxLoop(i32),
    /// A base other than A/C/G/T where only those are allowed (the
    /// [`nn`] duplex model; thal itself treats other bases as unpairable).
    InvalidBase(char),
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
            ThalError::InvalidBase(b) => write!(f, "invalid base {b:?}: only A, C, G and T"),
        }
    }
}

impl std::error::Error for ThalError {}

/// The most stable structure found and its thermodynamics, as primer3-py's
/// `ThermoResult`.
///
/// When no structure forms, `structure_found` is false and the numbers are
/// zero. Under extreme conditions (e.g. zero salt) primer3 can report
/// non-finite values; they are passed through unchanged.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ThermoResult {
    /// Whether any structure was found.
    pub structure_found: bool,
    /// Melting temperature, °C.
    pub tm: f64,
    /// Gibbs free energy at [`ThalArgs::temp`], cal/mol.
    pub dg: f64,
    /// Enthalpy, cal/mol.
    pub dh: f64,
    /// Entropy (salt-corrected), cal/K/mol.
    pub ds: f64,
    /// primer3's ASCII drawing of the structure, if requested and found.
    pub ascii_structure: Option<String>,
}

/// Compiles and runs the README's example as a doctest.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

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
