//! Hairpin alignment (thal.c `thal_hairpin`): fillMatrix2, calc_terminal_bp,
//! tracebacku, calcHairpin.

use super::{Context, Mode, draw, encode, no_structure, salt_correction};
use crate::thermo::{Init, Thermo};
use crate::{ILAH, ILAS, KELVIN_OFFSET, MIN_HRPN_LOOP, SMALL_NON_ZERO, TINY_ENTROPY, ThalArgs, ThermoResult};

pub(super) fn calc(seq: &str, args: &ThalArgs, output_structure: bool) -> ThermoResult {
    let (oligo, s1) = encode(seq, false);
    let init = Init {
        s: -TINY_ENTROPY,
        h: 0.0,
        rc: 0.0,
    };
    let len = s1.len() - 2;
    let mut ctx = Context::new(s1.clone(), s1, len - 1, init);
    let max_loop = args.max_loop as usize;

    ctx.init_hairpin_matrix();
    ctx.fill_hairpin_matrix(max_loop);
    let (ends, bases) = ctx.terminal_bp(args.temp);

    let best = ends[len];
    if !best.is_finite() {
        return no_structure();
    }
    let bp = ctx.traceback_hairpin(&ends, &bases, max_loop);

    // calcHairpin: thal.c counts paired bases among the first len - 1.
    let mut result = ThermoResult {
        structure_found: true,
        ascii_structure: output_structure.then(|| draw::hairpin(&bp, &oligo)),
        ..ThermoResult::default()
    };
    if best.s.is_finite() {
        let pairs = bp[..len - 1].iter().filter(|&&p| p > 0).count();
        let salt = ((pairs / 2) as f64 - 1.0) * salt_correction(args);
        result.tm = (best.h / (best.s + salt)) - KELVIN_OFFSET;
        result.dg = best.h - (args.temp * (best.s + salt));
        result.dh = best.h;
        result.ds = best.s + salt;
    }
    result
}

/// The four ways a helix can close the 5' end at position `i`
/// (thal.c `END5_1` .. `END5_4`).
#[derive(Clone, Copy)]
enum End5 {
    /// Pair (k+1, i), no dangles.
    Plain,
    /// Pair (k+2, i) with a 5' dangle.
    Dangle5,
    /// Pair (k+1, i-1) with a 3' dangle.
    Dangle3,
    /// Pair (k+2, i-1) with a terminal mismatch.
    Mismatch,
}

impl End5 {
    const ALL: [End5; 4] = [End5::Plain, End5::Dangle5, End5::Dangle3, End5::Mismatch];

    /// Values of `k` thal.c tries for position `i`.
    fn ks(self, i: usize) -> std::ops::Range<usize> {
        let skip = match self {
            End5::Plain => MIN_HRPN_LOOP + 1,
            End5::Dangle5 | End5::Dangle3 => MIN_HRPN_LOOP + 2,
            End5::Mismatch => MIN_HRPN_LOOP + 3,
        };
        0..i.saturating_sub(skip)
    }

    /// The closing pair for position `i` and split point `k`.
    fn pair(self, i: usize, k: usize) -> (usize, usize) {
        match self {
            End5::Plain => (k + 1, i),
            End5::Dangle5 => (k + 2, i),
            End5::Dangle3 => (k + 1, i - 1),
            End5::Mismatch => (k + 2, i - 1),
        }
    }
}

/// Inner pairs (ii, jj) that close a bulge/internal loop with (i, j), in
/// thal.c `CBI` order: widest inner pair first.
fn inner_pairs(i: usize, j: usize, max_loop: usize) -> impl Iterator<Item = (usize, usize)> {
    let span = j - i;
    let min_d = (MIN_HRPN_LOOP + 1).max(span.saturating_sub(2).saturating_sub(max_loop));
    (min_d..=span.saturating_sub(3))
        .rev()
        .flat_map(move |d| (i + 1..j - d).map(move |ii| (ii, ii + d)))
}

/// Pair-only table terms of a bulge/internal loop's closing pair.
struct Closing {
    atp: Thermo,
    stack_mm: Thermo,
    tstack: Thermo,
}

/// thal.c's test for keeping a bulge/internal candidate with Tm `t1` over a
/// cell with Tm `t2`. 1x1 loops within SMALL_NON_ZERO are skipped when filling.
fn accept(mode: Mode, is_1x1: bool, t1: f64, t2: f64) -> bool {
    if is_1x1 && mode == Mode::Fill && (t1 - t2) < SMALL_NON_ZERO {
        return false;
    }
    match mode {
        Mode::Fill => t1 > t2,
        Mode::Trace => t1 >= t2,
    }
}

/// Index of the largest of five Tm values, earliest on ties (thal.c `max5`).
fn argmax5(t: [f64; 5]) -> usize {
    let [a, b, c, d, e] = t;
    if a > b && a > c && a > d && a > e {
        0
    } else if b > c && b > d && b > e {
        1
    } else if c > d && c > e {
        2
    } else if d > e {
        3
    } else {
        4
    }
}

impl Context {
    fn init_hairpin_matrix(&mut self) {
        for i in 1..=self.len1 {
            for j in i..=self.len2 {
                let v = if j - i < MIN_HRPN_LOOP + 1 || !self.s1[i].pairs_with(self.s1[j]) {
                    Thermo::INVALID
                } else {
                    Thermo::PAIRABLE
                };
                self.set_dp(i, j, v);
            }
        }
    }

    fn fill_hairpin_matrix(&mut self, max_loop: usize) {
        for j in 2..=self.len2 {
            for i in (1..=j.saturating_sub(MIN_HRPN_LOOP + 1)).rev() {
                if !self.dp(i, j).is_finite() {
                    continue;
                }
                self.max_tm2(i, j);
                self.cbi_fill(i, j, max_loop);
                let hp = self.hairpin_loop(i, j, Mode::Fill);
                if hp.is_finite() {
                    self.set_dp(i, j, hp.clamp_min_entropy());
                }
            }
        }
    }

    /// Stack of pair (i, j) on (i + 1, j - 1) (thal.c `Ss`/`Hs` with k = 2).
    fn inner_stack(&self, i: usize, j: usize) -> Thermo {
        if i >= j || i == self.len1 || j == self.len2 + 1 {
            return Thermo::INVALID;
        }
        let s = &self.s1;
        self.p.stack(s[i], s[i + 1], s[j], s[j - 1])
    }

    /// thal.c `maxTM2`.
    fn max_tm2(&mut self, i: usize, j: usize) {
        let mut c0 = self.dp(i, j);
        let t0 = c0.tm(&self.init);
        let mut c1 = if c0.is_finite() {
            self.dp(i + 1, j - 1) + self.inner_stack(i, j)
        } else {
            Thermo::INVALID
        };
        let t1 = c1.tm(&self.init);
        c1 = c1.clamp_min_entropy();
        c0 = c0.clamp_min_entropy();
        self.set_dp(i, j, if t1 > t0 { c1 } else { c0 });
    }

    /// Fill (i, j) with its best bulge/internal loop (thal.c `CBI`, traceback 0).
    fn cbi_fill(&mut self, i: usize, j: usize, max_loop: usize) {
        let closing = self.closing(i, j);
        // The cell and its Tm only change when the cell is written.
        let mut cell = self.dp(i, j);
        let mut t_cell = cell.tm(&self.init);
        for (ii, jj) in inner_pairs(i, j, max_loop) {
            if ii > self.len1 || !(self.dp(ii, jj).is_finite() && cell.is_finite()) {
                continue;
            }
            let (c, is_1x1) = self.loop_energy2(i, j, ii, jj, &closing, true);
            if accept(Mode::Fill, is_1x1, c.tm(&self.init), t_cell) && c.is_finite() {
                cell = c.clamp_min_entropy();
                self.set_dp(i, j, cell);
                t_cell = cell.tm(&self.init);
            }
        }
    }

    /// The loop energy traceback compares against (thal.c `CBI`, traceback 2):
    /// the last candidate that ties or beats the cell.
    fn cbi_trace(&self, i: usize, j: usize, max_loop: usize) -> Thermo {
        let closing = self.closing(i, j);
        let cell = self.dp(i, j);
        let t_cell = cell.tm(&self.init);
        let mut acc = Thermo::INVALID;
        for (ii, jj) in inner_pairs(i, j, max_loop) {
            if ii > self.len1 || !(self.dp(ii, jj).is_finite() && cell.is_finite()) {
                continue;
            }
            let (c, is_1x1) = self.loop_energy2(i, j, ii, jj, &closing, true);
            if accept(Mode::Trace, is_1x1, c.tm(&self.init), t_cell) {
                acc = c;
            }
            if acc.is_finite() {
                acc = acc.clamp_min_entropy();
            }
        }
        acc
    }

    /// Hairpin loop closed by (i, j) (thal.c `calc_hairpin`). In `Fill` mode
    /// the current cell is kept instead if it has lower ΔG.
    fn hairpin_loop(&self, i: usize, j: usize, mode: Mode) -> Thermo {
        let loop_size = j - i - 1;
        if loop_size < MIN_HRPN_LOOP {
            return Thermo::INVALID;
        }
        let (p, s) = (self.p, &self.s1);
        let mut c = p.hairpin_loop[loop_size.min(30) - 1];
        if loop_size > 3 {
            c = c + p.tstack2(s[i], s[i + 1], s[j], s[j - 1]);
        } else {
            c = c + p.atp(s[i], s[j]);
        }
        let bonus = match loop_size {
            3 => p.triloop(&s[i..i + 5]),
            4 => p.tetraloop(&s[i..i + 6]),
            _ => None,
        };
        if let Some(b) = bonus {
            c = c + b;
        }
        if !c.is_finite() {
            c = Thermo::INVALID;
        }
        let cur = self.dp(i, j);
        // Negated comparisons: a NaN cell counts as "not positive", as in thal.c.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if c.h > 0.0 && c.s > 0.0 && (!(cur.h > 0.0) || !(cur.s > 0.0)) {
            c = Thermo::INVALID;
        }
        if mode == Mode::Fill {
            let end = self.right_end(i, j);
            if (cur + end).dg37() < (c + end).dg37() {
                return cur;
            }
        }
        c
    }

    /// Table terms of a loop's closing pair (i, j) that don't depend on
    /// the inner pair.
    fn closing(&self, i: usize, j: usize) -> Closing {
        let (p, s) = (self.p, &self.s1);
        Closing {
            atp: p.atp(s[i], s[j]),
            stack_mm: p.stack_mm(s[i], s[i + 1], s[j], s[j - 1]),
            tstack: p.tstack(s[i], s[i + 1], s[j], s[j - 1]),
        }
    }

    /// Bulge or internal loop between closing pair (i, j) and inner pair
    /// (ii, jj), plus the inner cell if `with_inner` (thal.c
    /// `calc_bulge_internal2` before its Tm test). Also reports a 1x1 loop.
    fn loop_energy2(
        &self,
        i: usize,
        j: usize,
        ii: usize,
        jj: usize,
        closing: &Closing,
        with_inner: bool,
    ) -> (Thermo, bool) {
        let (l1, l2) = (ii - i - 1, j - jj - 1);
        let loop_size = l1 + l2 - 1;
        let (p, s) = (self.p, &self.s1);
        let with_inner = |c: Thermo| if with_inner { c + self.dp(ii, jj) } else { c };

        let mut is_1x1 = false;
        let c = if (l1 == 0 && l2 > 0) || (l2 == 0 && l1 > 0) {
            if l1 == 1 || l2 == 1 {
                with_inner(p.bulge_loop[loop_size] + p.stack(s[i], s[ii], s[j], s[jj]))
            } else {
                with_inner(p.bulge_loop[loop_size] + closing.atp + p.atp(s[ii], s[jj]))
            }
        } else if l1 == 1 && l2 == 1 {
            is_1x1 = true;
            with_inner(closing.stack_mm + p.stack_mm(s[jj], s[jj + 1], s[ii], s[ii - 1]))
        } else {
            let asym = l1.abs_diff(l2) as f64;
            with_inner(
                p.interior_loop[loop_size]
                    + closing.tstack
                    + p.tstack(s[jj], s[jj + 1], s[ii], s[ii - 1])
                    + Thermo::new(ILAS * asym, ILAH * asym),
            )
        };
        let c = if c.is_finite() { c } else { Thermo::INVALID };
        (c, is_1x1)
    }

    /// Best 5' end energy at each position (thal.c `calc_terminal_bp`).
    ///
    /// Also returns, per position `k`, the value thal.c builds on when
    /// extending past `k`: `ends[k]`, or zero if that has lower Tm.
    fn terminal_bp(&self, temp: f64) -> (Vec<Thermo>, Vec<Thermo>) {
        let init = self.init;
        let zero_tm = Thermo::ZERO.tm(&init);
        let base = |e: Thermo| if e.tm(&init) >= zero_tm { e } else { Thermo::ZERO };
        let mut ends = vec![Thermo::PAIRABLE; self.len1 + 1];
        ends[0] = Thermo::INVALID;
        ends[1] = Thermo::INVALID;
        let mut bases = Vec::with_capacity(self.len1 + 1);
        bases.push(base(ends[0]));
        for i in 2..=self.len1 {
            // end5(i, ..) reads bases[..i - 4], all final by now.
            bases.push(base(ends[i - 1]));
            let cands = End5::ALL.map(|kind| self.end5(i, kind, &bases));
            let tms = [
                ends[i - 1].tm(&init),
                cands[0].tm(&init),
                cands[1].tm(&init),
                cands[2].tm(&init),
                cands[3].tm(&init),
            ];
            ends[i] = match argmax5(tms) {
                0 => ends[i - 1],
                k => {
                    let c = cands[k - 1];
                    if c.h - (temp * c.s) < 0.0 { c } else { ends[i - 1] }
                }
            };
        }
        bases.push(base(ends[self.len1]));
        (ends, bases)
    }

    /// Candidate energy for closing at `i` with split `k`, on top of `base`.
    fn end5_candidate(&self, i: usize, kind: End5, k: usize, base: Thermo) -> Thermo {
        let s = &self.s1;
        let (lo, hi) = kind.pair(i, k);
        let mut c = base + self.p.atp(s[lo], s[hi]);
        match kind {
            End5::Plain => {}
            End5::Dangle5 => c = c + self.p.dangle5(s[i], s[lo], s[lo - 1]),
            End5::Dangle3 => c = c + self.p.dangle3(s[i - 1], s[i], s[lo]),
            End5::Mismatch => c = c + self.p.tstack2(s[i - 1], s[i], s[lo], s[lo - 1]),
        }
        c + self.dp(lo, hi)
    }

    /// Highest-Tm way of closing at `i` with `kind` (thal.c `END5_*`).
    fn end5(&self, i: usize, kind: End5, bases: &[Thermo]) -> Thermo {
        let mut best = Thermo::INVALID;
        let mut best_tm = f64::NEG_INFINITY;
        for k in kind.ks(i) {
            let mut c = self.end5_candidate(i, kind, k, bases[k]);
            if !c.is_finite() || c.h > 0.0 || c.s > 0.0 {
                c = Thermo::INVALID;
            }
            let tm = c.tm(&self.init);
            if best_tm < tm && c.s > crate::MIN_ENTROPY_CUTOFF {
                best = c;
                best_tm = tm;
            }
        }
        best
    }

    /// Paired partner (1-based, 0 = unpaired) of each base (thal.c `tracebacku`).
    fn traceback_hairpin(&self, ends: &[Thermo], bases: &[Thermo], max_loop: usize) -> Vec<usize> {
        enum Frame {
            End5(usize),
            Pair(usize, usize),
        }
        let mut bp = vec![0; self.len1];
        let mut stack = vec![Frame::End5(self.len1)];
        while let Some(frame) = stack.pop() {
            match frame {
                Frame::End5(mut i) => {
                    while i > 0 && ends[i].approx_eq(ends[i - 1]) {
                        i -= 1;
                    }
                    if i == 0 {
                        continue;
                    }
                    let Some(kind) = End5::ALL
                        .into_iter()
                        .find(|&kind| ends[i].approx_eq(self.end5(i, kind, bases)))
                    else {
                        continue;
                    };
                    for k in kind.ks(i) {
                        let (lo, hi) = kind.pair(i, k);
                        if ends[i].approx_eq(self.end5_candidate(i, kind, k, Thermo::ZERO)) {
                            stack.push(Frame::Pair(lo, hi));
                            break;
                        }
                        if ends[i].approx_eq(self.end5_candidate(i, kind, k, ends[k])) {
                            stack.push(Frame::Pair(lo, hi));
                            stack.push(Frame::End5(k));
                            break;
                        }
                    }
                }
                Frame::Pair(i, j) => {
                    bp[i - 1] = j;
                    bp[j - 1] = i;
                    let cur = self.dp(i, j);
                    if cur.approx_eq(self.inner_stack(i, j) + self.dp(i + 1, j - 1)) {
                        stack.push(Frame::Pair(i + 1, j - 1));
                    } else if cur.approx_eq(self.hairpin_loop(i, j, Mode::Trace)) {
                        // Hairpin loop: nothing enclosed.
                    } else if cur.approx_eq(self.cbi_trace(i, j, max_loop)) {
                        let closing = self.closing(i, j);
                        let inner = inner_pairs(i, j, max_loop).find(|&(ii, jj)| {
                            let (c, _) = self.loop_energy2(i, j, ii, jj, &closing, false);
                            cur.approx_eq(c + self.dp(ii, jj))
                        });
                        if let Some((ii, jj)) = inner {
                            stack.push(Frame::Pair(ii, jj));
                        }
                    }
                }
            }
        }
        bp
    }
}
