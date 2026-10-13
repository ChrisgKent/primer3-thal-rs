//! thal's thermodynamic model for a perfectly matched duplex, exposed so
//! other code can compute Tm and annealing on exactly the same assumptions.
//!
//! thal uses these constants and corrections itself, and the Tm of a
//! [`Duplex`] equals thal's Tm for
//! [`calc_heterodimer`](crate::calc_heterodimer)`(seq, reverse_complement(seq))`
//! (tested to 1e-9 °C). Unlike thal it needs no alignment, so it is cheap
//! enough to call while extending a primer one base at a time.
//!
//! ```
//! use primer3_thal_rs::{ThalArgs, nn};
//!
//! let args = ThalArgs { mv: 100.0, dv: 2.0, dntp: 0.8, dna_conc: 15.0, ..ThalArgs::default() };
//! let tm = nn::duplex_tm("ACCAACGATGGTGTGTCCAT", &args)?;
//! // Half of the oligo is bound at its Tm.
//! let bound = nn::duplex_annealing("ACCAACGATGGTGTGTCCAT", &args, tm)?;
//! assert!((bound - 50.0).abs() < 1e-6);
//! # Ok::<(), primer3_thal_rs::ThalError>(())
//! ```

use std::sync::OnceLock;

use crate::params::{Base, params};
use crate::{ThalArgs, ThalError};

/// Gas constant, cal/K/mol.
pub const R: f64 = 1.9872;
/// 0 °C in Kelvin.
pub const KELVIN_OFFSET: f64 = 273.15;
/// Duplex initiation enthalpy, cal/mol.
pub const DPLX_INIT_H: f64 = 200.0;
/// Duplex initiation entropy, cal/K/mol.
pub const DPLX_INIT_S: f64 = -5.7;
/// Terminal A-T pair penalty enthalpy, cal/mol.
pub const AT_H: f64 = 2200.0;
/// Terminal A-T pair penalty entropy, cal/K/mol.
pub const AT_S: f64 = 6.9;
/// Strand-concentration divisor (nM) for a self-complementary duplex.
pub const SYMMETRIC_DNA_DIVISOR: f64 = 1_000_000_000.0;
/// Strand-concentration divisor (nM) for any other duplex.
pub const DNA_DIVISOR: f64 = 4_000_000_000.0;

/// thal's entropy salt correction per base pair beyond the first,
/// cal/K/mol (SantaLucia), with divalent cations converted using the
/// free-Mg2+ term `dv - dntp`.
pub fn salt_correction(args: &ThalArgs) -> f64 {
    let dntp = if args.dv <= 0.0 { args.dv } else { args.dntp };
    0.368 * ((args.mv + 120.0 * f64::max(0.0, args.dv - dntp).sqrt()) / 1000.0).ln()
}

/// Effective strand concentration of a duplex: `dna_conc` divided by
/// [`SYMMETRIC_DNA_DIVISOR`] or [`DNA_DIVISOR`].
pub fn strand_conc(dna_conc: f64, symmetric: bool) -> f64 {
    let divisor = if symmetric {
        SYMMETRIC_DNA_DIVISOR
    } else {
        DNA_DIVISOR
    };
    dna_conc / divisor
}

/// The concentration term in a duplex Tm, `R·ln(strand_conc)`.
pub fn dna_conc_term(dna_conc: f64, symmetric: bool) -> f64 {
    R * strand_conc(dna_conc, symmetric).ln()
}

/// Percentage of oligo bound at `temp_c` °C, from the two-state equilibrium
/// of a duplex with enthalpy `dh` (cal/mol) and salt-corrected entropy `ds`
/// (cal/K/mol), as thal reports them in [`ThermoResult`](crate::ThermoResult).
/// `symmetric` selects the strand concentration as in [`strand_conc`]. 50 at
/// the duplex's Tm before any DMSO/formamide correction, which does not
/// enter binding.
pub fn bound_percent(dh: f64, ds: f64, args: &ThalArgs, symmetric: bool, temp_c: f64) -> f64 {
    let t = temp_c + KELVIN_OFFSET;
    let ka = (-(dh - t * ds) / (R * t)).exp();
    100.0 / (1.0 + (1.0 / (strand_conc(args.dna_conc, symmetric) * ka)).sqrt())
}

/// DMSO and formamide Tm correction, °C, as primer3's oligotm:
/// `-dmso_conc·dmso_fact + (0.453·gc_fraction - 2.88)·formamide_conc`.
/// Zero with the default [`ThalArgs`].
pub fn tm_correction(gc_fraction: f64, args: &ThalArgs) -> f64 {
    -(args.dmso_conc * args.dmso_fact) + (0.453 * gc_fraction - 2.88) * args.formamide_conc
}

/// thal.c `symmetry_thermo`: an even-length sequence equal to its own
/// reverse complement, case-insensitively. Non-ACGT characters only match
/// each other.
pub fn is_self_complementary(seq: impl AsRef<[u8]>) -> bool {
    fn complement(b: u8) -> Option<u8> {
        match b {
            b'A' => Some(b'T'),
            b'T' => Some(b'A'),
            b'C' => Some(b'G'),
            b'G' => Some(b'C'),
            _ => None,
        }
    }
    let seq = seq.as_ref();
    seq.len().is_multiple_of(2)
        && seq
            .iter()
            .zip(seq.iter().rev())
            .take(seq.len() / 2)
            .all(|(&s, &e)| {
                let (s, e) = (s.to_ascii_uppercase(), e.to_ascii_uppercase());
                match (complement(s), complement(e)) {
                    (Some(c), _) => c == e,
                    (None, other) => other.is_none(),
                }
            })
}

/// A base of a duplex, or `InvalidBase` for anything but A/C/G/T.
fn base(b: u8) -> Result<Base, ThalError> {
    match b.to_ascii_uppercase() {
        b'A' => Ok(Base::A),
        b'C' => Ok(Base::C),
        b'G' => Ok(Base::G),
        b'T' => Ok(Base::T),
        _ => Err(ThalError::InvalidBase(b as char)),
    }
}

/// Watson-Crick nearest-neighbour stacks `[a][b]` = (ΔH, ΔS) of 5'-ab-3'
/// on its complement, read from thal's own stack table.
fn wc_stacks() -> &'static [[(f64, f64); 4]; 4] {
    static STACKS: OnceLock<[[(f64, f64); 4]; 4]> = OnceLock::new();
    STACKS.get_or_init(|| {
        let p = params();
        let bases = [Base::A, Base::C, Base::G, Base::T];
        let mut out = [[(0.0, 0.0); 4]; 4];
        for (i, &a) in bases.iter().enumerate() {
            for (j, &b) in bases.iter().enumerate() {
                let t = p.stack(a, b, a.complement(), b.complement());
                out[i][j] = (t.h, t.s);
            }
        }
        out
    })
}

/// A perfectly matched duplex of an oligo with its complement, built one
/// base at a time (5'->3'): the running nearest-neighbour sums, so adding a
/// base is O(1).
#[derive(Clone, Copy, Debug, Default)]
pub struct Duplex {
    dh: f64,
    ds: f64,
    len: usize,
    gc: usize,
    first: Option<Base>,
    last: Option<Base>,
}

impl Duplex {
    /// An empty duplex.
    pub fn new() -> Self {
        Self::default()
    }

    /// The duplex of a whole sequence.
    pub fn from_seq(seq: impl AsRef<[u8]>) -> Result<Self, ThalError> {
        let mut d = Self::new();
        for &b in seq.as_ref() {
            d.push(b)?;
        }
        Ok(d)
    }

    /// Extend the oligo by one base at its 3' end. Only A/C/G/T are
    /// accepted (case-insensitively); the duplex is unchanged on error.
    pub fn push(&mut self, b: u8) -> Result<(), ThalError> {
        let b = base(b)?;
        if let Some(prev) = self.last {
            let (h, s) = wc_stacks()[prev as usize][b as usize];
            self.dh += h;
            self.ds += s;
        }
        if matches!(b, Base::C | Base::G) {
            self.gc += 1;
        }
        self.first.get_or_insert(b);
        self.last = Some(b);
        self.len += 1;
        Ok(())
    }

    /// Number of bases.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether no bases have been added.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Total ΔH (cal/mol) and salt-corrected ΔS (cal/K/mol): the stacks,
    /// initiation and terminal A-T penalties, as thal.
    fn thermo(&self, args: &ThalArgs) -> (f64, f64) {
        let mut dh = self.dh + DPLX_INIT_H;
        let mut ds = self.ds + DPLX_INIT_S;
        for end in [self.first, self.last].into_iter().flatten() {
            if matches!(end, Base::A | Base::T) {
                dh += AT_H;
                ds += AT_S;
            }
        }
        let pairs = self.len.saturating_sub(1) as f64;
        (dh, ds + pairs * salt_correction(args))
    }

    /// Melting temperature, °C, including the DMSO/formamide correction.
    /// `symmetric` is whether the oligo is self-complementary
    /// ([`is_self_complementary`]), which halves the strand concentration
    /// divisor as in thal. NaN if empty.
    pub fn tm(&self, args: &ThalArgs, symmetric: bool) -> f64 {
        if self.is_empty() {
            return f64::NAN;
        }
        let (dh, ds) = self.thermo(args);
        let gc_fraction = self.gc as f64 / self.len as f64;
        dh / (ds + dna_conc_term(args.dna_conc, symmetric)) - KELVIN_OFFSET
            + tm_correction(gc_fraction, args)
    }

    /// Percentage of the oligo bound to its complement at `temp_c` °C,
    /// from the two-state equilibrium at the same ΔH/ΔS and strand
    /// concentration as [`Duplex::tm`] (50 at the uncorrected Tm). DMSO and
    /// formamide do not enter it. NaN if empty.
    pub fn annealing(&self, args: &ThalArgs, symmetric: bool, temp_c: f64) -> f64 {
        if self.is_empty() {
            return f64::NAN;
        }
        let (dh, ds) = self.thermo(args);
        bound_percent(dh, ds, args, symmetric, temp_c)
    }
}

/// Melting temperature (°C) of `seq` paired with its complement: the same
/// value as thal's `calc_heterodimer(seq, reverse_complement(seq)).tm`,
/// including the DMSO/formamide correction.
pub fn duplex_tm(seq: impl AsRef<[u8]>, args: &ThalArgs) -> Result<f64, ThalError> {
    let seq = seq.as_ref();
    let d = non_empty(seq)?;
    Ok(d.tm(args, is_self_complementary(seq)))
}

/// Percentage of `seq` bound to its complement at `temp_c` °C (see
/// [`Duplex::annealing`]).
pub fn duplex_annealing(
    seq: impl AsRef<[u8]>,
    args: &ThalArgs,
    temp_c: f64,
) -> Result<f64, ThalError> {
    let seq = seq.as_ref();
    let d = non_empty(seq)?;
    Ok(d.annealing(args, is_self_complementary(seq), temp_c))
}

fn non_empty(seq: &[u8]) -> Result<Duplex, ThalError> {
    if seq.is_empty() {
        return Err(ThalError::EmptySequence);
    }
    Duplex::from_seq(seq)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wc_stacks_are_thals_stacks() {
        // AA/TT and CG/GC from SantaLucia 1998, as in thal's stack table
        assert_eq!(wc_stacks()[0][0], (-7900.0, -22.2));
        assert_eq!(wc_stacks()[1][2], (-10600.0, -27.2));
        let p = params();
        let bases = [Base::A, Base::C, Base::G, Base::T];
        for (i, &a) in bases.iter().enumerate() {
            for (j, &b) in bases.iter().enumerate() {
                let t = p.stack(a, b, a.complement(), b.complement());
                assert_eq!(wc_stacks()[i][j], (t.h, t.s));
                assert!(t.h.is_finite() && t.s.is_finite());
            }
        }
    }
}
