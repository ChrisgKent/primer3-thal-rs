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
    /// Hairpin traceback: the loop's own energy, without the enclosed cell.
    LoopOnly,
}

/// Uppercased bytes (for drawing) and the padded, encoded sequence.
fn encode(seq: &str, reverse: bool) -> (Vec<u8>, Vec<Base>) {
    let mut oligo: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    if reverse {
        oligo.reverse();
    }
    let mut bases = Vec::with_capacity(oligo.len() + 2);
    bases.push(Base::N);
    bases.extend(oligo.iter().map(|&b| Base::from_ascii(b)));
    bases.push(Base::N);
    (oligo, bases)
}

fn salt_correction(args: &ThalArgs) -> f64 {
    let dntp = if args.dv <= 0.0 { args.dv } else { args.dntp };
    0.368 * ((args.mv + 120.0 * f64::max(0.0, args.dv - dntp).sqrt()) / 1000.0).ln()
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

fn no_structure() -> ThermoResult {
    ThermoResult::default()
}

/// Hairpin (monomer secondary structure), primer3's `thal_hairpin` mode.
pub fn calc_hairpin(
    seq: &str,
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

/// Heterodimer, primer3's `thal_any` mode.
pub fn calc_heterodimer(
    seq1: &str,
    seq2: &str,
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

/// Homodimer: primer3-py runs `thal(seq, seq, thal_any)`.
pub fn calc_homodimer(
    seq: &str,
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    calc_heterodimer(seq, seq, args, output_structure)
}

#[cfg(feature = "parallel")]
fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync + Send) -> Vec<R> {
    items.par_iter().map(f).collect()
}

#[cfg(feature = "parallel")]
pub fn calc_hairpin_batch<T: AsRef<str> + Sync>(
    seqs: &[T],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    par_map(seqs, |s| calc_hairpin(s.as_ref(), args, output_structure))
}

#[cfg(feature = "parallel")]
pub fn calc_homodimer_batch<T: AsRef<str> + Sync>(
    seqs: &[T],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    par_map(seqs, |s| calc_homodimer(s.as_ref(), args, output_structure))
}

#[cfg(feature = "parallel")]
pub fn calc_heterodimer_batch(
    pairs: &[(String, String)],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    par_map(pairs, |(a, b)| calc_heterodimer(a, b, args, output_structure))
}
