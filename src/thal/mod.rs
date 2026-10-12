//! Dynamic programming from primer3's thal.c.
//!
//! Indices follow thal.c: sequences are 1-based, padded with `Base::N` at
//! position 0 and `len + 1`, and the DP matrices are row-major with a
//! mode-dependent row stride (see [`Context::new`]). Every comparison and
//! every floating-point expression keeps thal.c's order so results match
//! primer3 bit-for-bit (to within last-ulp FMA effects in the C build).

mod dimer;
mod draw;
mod hairpin;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::nn::{salt_correction, tm_correction};
use crate::params::{Base, ThermoParams, params};
use crate::thermo::{Init, Thermo};
use crate::{ThalArgs, ThalError, ThermoResult};

/// How a loop-energy routine is being used.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    /// Filling the DP table: keep a candidate only if it beats the cell.
    Fill,
    /// Re-deriving a filled cell during traceback.
    Trace,
}

/// Encode `seq` for the DP: 1-based, with an `N` sentinel at each end,
/// optionally reversed (the second strand of a dimer is aligned 3'->5').
fn encode(seq: &[u8], reverse: bool) -> Vec<Base> {
    let mut bases = Vec::with_capacity(seq.len() + 2);
    bases.push(Base::N);
    if reverse {
        bases.extend(seq.iter().rev().map(|&b| Base::from_ascii(b)));
    } else {
        bases.extend(seq.iter().map(|&b| Base::from_ascii(b)));
    }
    bases.push(Base::N);
    bases
}

/// Uppercased (and optionally reversed) bytes, as thal.c prints them in
/// structure drawings. Only built when a structure is requested.
fn display_bytes(seq: &[u8], reverse: bool) -> Vec<u8> {
    let mut out: Vec<u8> = seq.iter().map(u8::to_ascii_uppercase).collect();
    if reverse {
        out.reverse();
    }
    out
}

/// GC fraction of all the given sequences together, for the DMSO/formamide
/// Tm correction.
fn gc_fraction(seqs: &[&[u8]]) -> f64 {
    let len: usize = seqs.iter().map(|s| s.len()).sum();
    let gc = seqs
        .iter()
        .flat_map(|s| s.iter())
        .filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / len as f64
}

/// The two sequences and the entropy/enthalpy DP tables thal fills.
struct Context {
    p: &'static ThermoParams,
    s1: Vec<Base>,
    s2: Vec<Base>,
    len1: usize,
    len2: usize,
    /// Row stride of the DP tables: `len2` for dimers, `len2 - 1` for
    /// hairpins. The hairpin stride makes the last cell of row `i` alias the
    /// first of row `i + 1`; thal.c does the same and results depend on it.
    stride: usize,
    entropy: Vec<f64>,
    enthalpy: Vec<f64>,
    init: Init,
    /// Dimers only: `right_end(i, j)` for every cell, row stride `len2`.
    right_ends: Vec<Thermo>,
    /// Dimers only: for each anti-diagonal `ii + jj`, the ascending `ii` of
    /// cells whose bases can pair, flattened; `pairable_off[sum]` indexes
    /// the start of diagonal `sum`.
    pairable_ii: Vec<u32>,
    pairable_off: Vec<u32>,
}

impl Context {
    fn new(s1: Vec<Base>, s2: Vec<Base>, stride: usize, init: Init) -> Self {
        let len1 = s1.len() - 2;
        let len2 = s2.len() - 2;
        // Zero-initialised: hairpin cells below the diagonal are never
        // written and must read as 0, as in thal.c.
        Self {
            p: params(),
            s1,
            s2,
            len1,
            len2,
            stride,
            entropy: vec![0.0; len1 * len2],
            enthalpy: vec![0.0; len1 * len2],
            init,
            right_ends: Vec::new(),
            pairable_ii: Vec::new(),
            pairable_off: Vec::new(),
        }
    }

    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        (i - 1) * self.stride + j - 1
    }

    #[inline]
    fn dp(&self, i: usize, j: usize) -> Thermo {
        let k = self.idx(i, j);
        Thermo::new(self.entropy[k], self.enthalpy[k])
    }

    #[inline]
    fn set_dp(&mut self, i: usize, j: usize, v: Thermo) {
        let k = self.idx(i, j);
        self.entropy[k] = v.s;
        self.enthalpy[k] = v.h;
    }

    /// thal.c `RSH`: the helix end at pair (i, j) facing (i + 1, j + 1).
    fn right_end(&self, i: usize, j: usize) -> Thermo {
        self.terminal(self.s1[i], self.s1[i + 1], self.s2[j], self.s2[j + 1])
    }

    /// thal.c `LSH`: the helix end at pair (i, j) facing (i - 1, j - 1).
    fn left_end(&self, i: usize, j: usize) -> Thermo {
        self.terminal(self.s2[j], self.s2[j - 1], self.s1[i], self.s1[i - 1])
    }

    /// Best of terminal mismatch, dangling end(s) or a bare AT penalty for
    /// the closing pair `a`-`b` with neighbours `a1`/`b1` (thal.c RSH/LSH).
    fn terminal(&self, a: Base, a1: Base, b: Base, b1: Base) -> Thermo {
        if !a.pairs_with(b) {
            return Thermo::INVALID;
        }
        let p = self.p;
        let init = &self.init;
        let atp = p.atp(a, b);

        let mut best = atp + p.tstack2(a, a1, b, b1);
        let mut g1 = best.dg37();
        if !best.is_finite() || g1 > 0.0 {
            best = Thermo::INVALID;
            g1 = 1.0;
        }
        let mut t1 = f64::NEG_INFINITY;

        if !a1.pairs_with(b1) {
            let d3 = p.dangle3(a, a1, b);
            let d5 = p.dangle5(a, b, b1);
            let dangle = match (d3.is_finite(), d5.is_finite()) {
                (true, true) => Some(atp + d3 + d5),
                (true, false) => Some(atp + d3),
                (false, true) => Some(atp + d5),
                (false, false) => None,
            };
            if let Some(mut cand) = dangle {
                let mut g2 = cand.dg37();
                if !cand.is_finite() || g2 > 0.0 {
                    cand = Thermo::INVALID;
                    g2 = 1.0;
                }
                let t2 = cand.tm(init);
                if best.is_finite() && g1 < 0.0 {
                    t1 = best.tm(init);
                    if t1 < t2 && g2 < 0.0 {
                        best = cand;
                        t1 = t2;
                    }
                } else if g2 < 0.0 {
                    best = cand;
                    t1 = t2;
                }
            }
        }

        // `t1` is still -inf if no dangle was considered, so a finite
        // terminal mismatch loses to the bare AT penalty unless the latter's
        // Tm is NaN (as for GC pairs in hairpins). That is thal.c's logic.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if best.is_finite() && !(t1 < atp.tm(init)) {
            best
        } else {
            atp
        }
    }
}

/// Result when no structure forms: `structure_found` false, all zero.
fn no_structure() -> ThermoResult {
    ThermoResult::default()
}

/// Hairpin (monomer secondary structure) thermodynamics of `seq`, as
/// primer3-py's `calc_hairpin` (thal.c `thal_hairpin` mode).
///
/// `seq` is any byte sequence, e.g. `&str`, `String`, `&[u8]` or `Vec<u8>`;
/// it is read 5'->3', case-insensitively, and any byte other than A/C/G/T
/// (or the digits 0-3, as in primer3) is treated as a base that pairs with
/// nothing. With `output_structure`, [`ThermoResult::ascii_structure`] holds
/// primer3's two-line drawing.
///
/// # Errors
///
/// [`ThalError::EmptySequence`], [`ThalError::TooLong`] if `seq` is longer
/// than 60 bases, or [`ThalError::InvalidMaxLoop`].
///
/// # Example
///
/// ```
/// use primer3_thal_rs::{ThalArgs, calc_hairpin};
///
/// let res = calc_hairpin("CCCCCATCCGATCAGGGGG", &ThalArgs::default(), false)?;
/// assert!(res.structure_found);
/// assert!((res.tm - 66.747).abs() < 1e-3);
/// # Ok::<(), primer3_thal_rs::ThalError>(())
/// ```
pub fn calc_hairpin(
    seq: impl AsRef<[u8]>,
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    hairpin_bytes(seq.as_ref(), args, output_structure)
}

/// Non-generic body of [`calc_hairpin`]: input checks, then the DP.
fn hairpin_bytes(
    seq: &[u8],
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    args.validate()?;
    if seq.is_empty() {
        return Err(ThalError::EmptySequence);
    }
    if seq.len() > crate::THAL_MAX_ALIGN {
        return Err(ThalError::TooLong {
            len: seq.len(),
            max: crate::THAL_MAX_ALIGN,
        });
    }
    Ok(hairpin::calc(seq, args, output_structure))
}

/// Heterodimer thermodynamics of `seq1` binding `seq2`, as primer3-py's
/// `calc_heterodimer` (thal.c `thal_any` mode).
///
/// Both sequences are given 5'->3', exactly as the primers are written; the
/// second is reversed internally. Sequence handling is as for
/// [`calc_hairpin`]. With `output_structure`,
/// [`ThermoResult::ascii_structure`] holds primer3's four-line duplex drawing.
///
/// # Errors
///
/// [`ThalError::EmptySequence`]; [`ThalError::TooLong`] if *both* sequences
/// are longer than 60 bases (one may be longer, up to 10,000), as in primer3;
/// or [`ThalError::InvalidMaxLoop`].
///
/// # Example
///
/// ```
/// use primer3_thal_rs::{ThalArgs, calc_heterodimer};
///
/// // Two primers with complementary 3' ends.
/// let args = ThalArgs::default();
/// let res = calc_heterodimer("GTAAAACGACGGCCAGT", b"CAGGAAACAGCTATGACACTGG", &args, false)?;
/// assert!(res.structure_found);
/// assert!(res.dg < 0.0);
/// # Ok::<(), primer3_thal_rs::ThalError>(())
/// ```
pub fn calc_heterodimer(
    seq1: impl AsRef<[u8]>,
    seq2: impl AsRef<[u8]>,
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    heterodimer_bytes(seq1.as_ref(), seq2.as_ref(), args, output_structure)
}

/// Non-generic body of [`calc_heterodimer`]: input checks, then the DP.
fn heterodimer_bytes(
    seq1: &[u8],
    seq2: &[u8],
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    args.validate()?;
    if seq1.is_empty() || seq2.is_empty() {
        return Err(ThalError::EmptySequence);
    }
    // As in thal.c: one sequence may be longer than THAL_MAX_ALIGN, not both.
    if seq1.len() > crate::THAL_MAX_ALIGN && seq2.len() > crate::THAL_MAX_ALIGN {
        return Err(ThalError::TooLong {
            len: seq1.len().min(seq2.len()),
            max: crate::THAL_MAX_ALIGN,
        });
    }
    if let Some(&len) = [seq1.len(), seq2.len()]
        .iter()
        .find(|&&len| len > crate::THAL_MAX_SEQ)
    {
        return Err(ThalError::TooLong {
            len,
            max: crate::THAL_MAX_SEQ,
        });
    }
    Ok(dimer::calc(seq1, seq2, args, output_structure))
}

/// Homodimer (self-dimer) thermodynamics of `seq`, as primer3-py's
/// `calc_homodimer`, which runs the heterodimer alignment of `seq` with
/// itself.
///
/// # Errors
///
/// As [`calc_heterodimer`]; a homodimer of more than 60 bases is
/// [`ThalError::TooLong`].
pub fn calc_homodimer(
    seq: impl AsRef<[u8]>,
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    let seq = seq.as_ref();
    heterodimer_bytes(seq, seq, args, output_structure)
}

#[cfg(feature = "parallel")]
fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync + Send) -> Vec<R> {
    items.par_iter().map(f).collect()
}

/// [`calc_hairpin`] for many sequences in parallel (rayon), in input order.
#[cfg(feature = "parallel")]
pub fn calc_hairpin_batch<T: AsRef<[u8]> + Sync>(
    seqs: &[T],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    par_map(seqs, |s| calc_hairpin(s.as_ref(), args, output_structure))
}

/// [`calc_homodimer`] for many sequences in parallel (rayon), in input order.
#[cfg(feature = "parallel")]
pub fn calc_homodimer_batch<T: AsRef<[u8]> + Sync>(
    seqs: &[T],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    par_map(seqs, |s| calc_homodimer(s.as_ref(), args, output_structure))
}

/// [`calc_heterodimer`] for many pairs in parallel (rayon), in input order.
#[cfg(feature = "parallel")]
pub fn calc_heterodimer_batch<A: AsRef<[u8]> + Sync, B: AsRef<[u8]> + Sync>(
    pairs: &[(A, B)],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    par_map(pairs, |(a, b)| {
        calc_heterodimer(a, b, args, output_structure)
    })
}
