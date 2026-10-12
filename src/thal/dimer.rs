//! Dimer alignment (thal.c `thal_any`): fillMatrix, traceback, calcDimer.

use super::{
    Context, display_bytes, draw, encode, gc_fraction, no_structure, salt_correction, tm_correction,
};
use crate::nn::{dna_conc_term, is_self_complementary};
use crate::thermo::{Init, Thermo};
use crate::{
    DPLX_INIT_H, DPLX_INIT_S, ILAH, ILAS, KELVIN_OFFSET, SMALL_NON_ZERO, TEMP_KELVIN, ThalArgs,
    ThermoResult,
};

pub(super) fn calc(
    seq1: &[u8],
    seq2: &[u8],
    args: &ThalArgs,
    output_structure: bool,
) -> ThermoResult {
    let s1 = encode(seq1, false);
    // seq2 is aligned 3'->5' against seq1.
    let s2 = encode(seq2, true);
    let symmetric = is_self_complementary(seq1) && is_self_complementary(seq2);
    let rc = dna_conc_term(args.dna_conc, symmetric);
    let init = Init {
        s: DPLX_INIT_S,
        h: DPLX_INIT_H,
        rc,
    };
    let stride = s2.len() - 2;
    let mut ctx = Context::new(s1, s2, stride, init);
    let max_loop = args.max_loop as usize;

    ctx.cache_right_ends();
    ctx.cache_pairable();
    ctx.init_dimer_matrix();
    ctx.fill_dimer_matrix(max_loop);

    let (bi, bj) = ctx.best_any().unwrap_or((1, 1));
    let best = ctx.dp(bi, bj);
    if !best.is_finite() {
        return no_structure();
    }
    let end = ctx.end_at(bi, bj);
    let dh = best.h + end.h + init.h;
    let ds = best.s + end.s + init.s;

    let (ps1, ps2) = ctx.traceback_dimer(bi, bj, max_loop);
    let ascii_structure = output_structure.then(|| {
        draw::dimer(
            &ps1,
            &ps2,
            &display_bytes(seq1, false),
            &display_bytes(seq2, true),
        )
    });

    if !args.temp.is_finite() {
        return ThermoResult {
            structure_found: true,
            ascii_structure,
            ..ThermoResult::default()
        };
    }
    // calcDimer: salt correction counts one per base pair beyond the first.
    let pairs = ps1.iter().chain(&ps2).filter(|&&p| p > 0).count() as i32;
    let n = (pairs / 2 - 1) as f64;
    let salt = salt_correction(args);
    ThermoResult {
        structure_found: true,
        tm: (dh / (ds + (n * salt) + rc)) - KELVIN_OFFSET
            + tm_correction(gc_fraction(&[seq1, seq2]), args),
        dg: dh - (args.temp * (ds + (n * salt))),
        dh,
        ds: ds + (n * salt),
        ascii_structure,
    }
}

/// Cells (ii, jj) on thal.c's anti-diagonal `d` below (i, j), nearest first.
fn diagonal(i: usize, j: usize, d: usize) -> impl Iterator<Item = (usize, usize)> {
    let (i, j, d) = (i as isize, j as isize, d as isize);
    let (mut ii, mut jj) = (i - 1, j + 1 - d);
    if jj < 1 {
        ii -= 1 - jj;
        jj = 1;
    }
    std::iter::from_fn(move || {
        (ii > 0 && jj < j).then(|| {
            let cell = (ii as usize, jj as usize);
            ii -= 1;
            jj += 1;
            cell
        })
    })
}

impl Context {
    /// `right_end` depends only on the sequences, but the fill needs it for
    /// every bulge/internal candidate; compute it once per cell.
    fn cache_right_ends(&mut self) {
        let mut ends = Vec::with_capacity(self.len1 * self.len2);
        for i in 1..=self.len1 {
            for j in 1..=self.len2 {
                ends.push(self.right_end(i, j));
            }
        }
        self.right_ends = ends;
    }

    /// Index the cells that can pair, by anti-diagonal. In the dimer DP a
    /// cell is finite exactly when its bases pair: `init_dimer_matrix` makes
    /// the rest `INVALID` and the fill never writes them. So the fill can
    /// skip them without changing any result.
    fn cache_pairable(&mut self) {
        let mut ii_out = Vec::new();
        let mut off = Vec::with_capacity(self.len1 + self.len2 + 2);
        for sum in 0..=self.len1 + self.len2 {
            off.push(ii_out.len() as u32);
            let lo = 1.max(sum.saturating_sub(self.len2));
            let hi = self.len1.min(sum.saturating_sub(1));
            for ii in lo..=hi {
                if self.s1[ii].pairs_with(self.s2[sum - ii]) {
                    ii_out.push(ii as u32);
                }
            }
        }
        off.push(ii_out.len() as u32);
        self.pairable_ii = ii_out;
        self.pairable_off = off;
    }

    /// Pairable `ii` on anti-diagonal `sum`, ascending.
    #[inline]
    fn pairable_on(&self, sum: usize) -> &[u32] {
        &self.pairable_ii[self.pairable_off[sum] as usize..self.pairable_off[sum + 1] as usize]
    }

    #[inline]
    fn end_at(&self, i: usize, j: usize) -> Thermo {
        self.right_ends[(i - 1) * self.len2 + j - 1]
    }

    fn init_dimer_matrix(&mut self) {
        for i in 1..=self.len1 {
            for j in 1..=self.len2 {
                let v = if self.s1[i].pairs_with(self.s2[j]) {
                    Thermo::PAIRABLE
                } else {
                    Thermo::INVALID
                };
                self.set_dp(i, j, v);
            }
        }
    }

    fn fill_dimer_matrix(&mut self, max_loop: usize) {
        for i in 1..=self.len1 {
            for j in 1..=self.len2 {
                if !self.dp(i, j).is_finite() {
                    continue;
                }
                let left = self.left_end(i, j);
                if left.is_finite() {
                    self.set_dp(i, j, left);
                }
                if i > 1 && j > 1 {
                    self.max_tm(i, j);
                    self.fill_bulge_internal(i, j, max_loop);
                }
            }
        }
    }

    /// Keep whichever of the current cell or a stack on (i-1, j-1) melts
    /// higher (thal.c `maxTM`).
    fn max_tm(&mut self, i: usize, j: usize) {
        let init = self.init;
        let mut c0 = self.dp(i, j);
        let end = self.end_at(i, j);
        let t0 = (c0.h + init.h + end.h) / (c0.s + init.s + end.s + init.rc);

        let prev = self.dp(i - 1, j - 1);
        let stack = self
            .p
            .stack(self.s1[i - 1], self.s1[i], self.s2[j - 1], self.s2[j]);
        let (mut c1, t1) = if prev.is_finite() && stack.is_finite() {
            let c1 = prev + stack;
            (
                c1,
                (c1.h + init.h + end.h) / (c1.s + init.s + end.s + init.rc),
            )
        } else {
            (Thermo::INVALID, Thermo::INVALID.tm(&init))
        };

        c1 = c1.clamp_min_entropy();
        c0 = c0.clamp_min_entropy();
        if t1 > t0 {
            self.set_dp(i, j, c1);
        } else if t0 >= t1 {
            self.set_dp(i, j, c0);
        }
    }

    /// Improve (i, j) with every bulge/internal loop that lowers its ΔG
    /// (the candidate loop of thal.c `fillMatrix`).
    fn fill_bulge_internal(&mut self, i: usize, j: usize, max_loop: usize) {
        let outer = self.outer(i, j);
        let end = self.end_at(i, j);
        // ΔG of the cell as thal.c's `calc_bulge_internal` compares it;
        // only changes when the cell is written.
        let mut g_cell = (self.dp(i, j) + end).dg37();
        for d in 3..=max_loop + 2 {
            // The cells of diagonal(i, j, d) that can pair, in the same
            // order (ii descending).
            let Some(sum) = (i + j).checked_sub(d).filter(|&s| s >= 2) else {
                break;
            };
            let lo = 1.max((sum + 1).saturating_sub(j));
            let hi = (i - 1).min(sum - 1);
            if lo > hi {
                continue;
            }
            let diag = self.pairable_on(sum);
            let start = diag.partition_point(|&ii| (ii as usize) < lo);
            let stop = diag.partition_point(|&ii| ii as usize <= hi);
            let base = self.pairable_off[sum] as usize;
            for k in (start..stop).rev() {
                let ii = self.pairable_ii[base + k] as usize;
                let jj = sum - ii;
                if !self.dp(ii, jj).is_finite() {
                    continue;
                }
                let c = self.loop_energy(ii, jj, &outer);
                if (c + end).dg37() < g_cell {
                    let c = c.clamp_min_entropy();
                    if c.is_finite() {
                        self.set_dp(i, j, c);
                        g_cell = (c + end).dg37();
                    }
                }
            }
        }
    }

    /// Table terms of the outer pair (ii, jj) shared by all its loops.
    fn outer(&self, ii: usize, jj: usize) -> Outer {
        let (p, s1, s2) = (self.p, &self.s1, &self.s2);
        Outer {
            ii,
            jj,
            atp: p.atp(s1[ii], s2[jj]),
            stack_mm: p.stack_mm(s2[jj], s2[jj - 1], s1[ii], s1[ii - 1]),
            tstack: p.tstack(s2[jj], s2[jj - 1], s1[ii], s1[ii - 1]),
        }
    }

    /// Bulge or internal loop closed by inner pair (i, j) and `outer`, plus
    /// the inner cell (thal.c `calc_bulge_internal`, before its ΔG test).
    fn loop_energy(&self, i: usize, j: usize, outer: &Outer) -> Thermo {
        let Outer { ii, jj, .. } = *outer;
        let (l1, l2) = (ii - i - 1, jj - j - 1);
        let loop_size = l1 + l2 - 1;
        let (p, s1, s2) = (self.p, &self.s1, &self.s2);
        let inner = self.dp(i, j);

        let mut c;
        if (l1 == 0 && l2 > 0) || (l2 == 0 && l1 > 0) {
            if l1 == 1 || l2 == 1 {
                c = p.bulge_loop[loop_size] + p.stack(s1[i], s1[ii], s2[j], s2[jj]);
                if c.h > 0.0 || c.s > 0.0 {
                    c = Thermo::INVALID;
                }
                c = c + inner;
                if !c.is_finite() {
                    c = Thermo::INVALID;
                }
            } else {
                c = p.bulge_loop[loop_size] + p.atp(s1[i], s2[j]) + outer.atp + inner;
                c = forbid_unfavourable(c);
            }
        } else if l1 == 1 && l2 == 1 {
            c = p.stack_mm(s1[i], s1[i + 1], s2[j], s2[j + 1]) + outer.stack_mm + inner;
            c = forbid_unfavourable(c);
        } else {
            let asym = l1.abs_diff(l2) as f64;
            c = p.interior_loop[loop_size]
                + p.tstack(s1[i], s1[i + 1], s2[j], s2[j + 1])
                + outer.tstack
                + Thermo::new(ILAS * asym, ILAH * asym)
                + inner;
            c = forbid_unfavourable(c);
        }
        c
    }

    /// Cell with the lowest ΔG at 37 °C including its right end.
    fn best_any(&self) -> Option<(usize, usize)> {
        let init = self.init;
        let mut best = None;
        let mut best_g = f64::INFINITY;
        for i in 1..=self.len1 {
            for j in 1..=self.len2 {
                let c = self.dp(i, j);
                let end = self.end_at(i, j);
                let (end_s, end_h) = (end.s + SMALL_NON_ZERO, end.h + SMALL_NON_ZERO);
                let g = (c.h + end_h + init.h) - TEMP_KELVIN * (c.s + end_s + init.s);
                if g < best_g {
                    best_g = g;
                    best = Some((i, j));
                }
            }
        }
        best.filter(|_| best_g.is_finite())
    }

    /// Paired positions (1-based partner, 0 = unpaired) for each strand.
    fn traceback_dimer(
        &self,
        mut i: usize,
        mut j: usize,
        max_loop: usize,
    ) -> (Vec<usize>, Vec<usize>) {
        let mut ps1 = vec![0; self.len1];
        let mut ps2 = vec![0; self.len2];
        ps1[i - 1] = j;
        ps2[j - 1] = i;
        loop {
            let cur = self.dp(i, j);
            if cur.approx_eq(self.left_end(i, j)) {
                break;
            }
            let stacked = i > 1
                && j > 1
                && cur.approx_eq(
                    self.p
                        .stack(self.s1[i - 1], self.s1[i], self.s2[j - 1], self.s2[j])
                        + self.dp(i - 1, j - 1),
                );
            let next = if stacked {
                Some((i - 1, j - 1))
            } else {
                let outer = self.outer(i, j);
                (3..=max_loop + 2)
                    .flat_map(|d| diagonal(i, j, d))
                    .find(|&(ii, jj)| cur.approx_eq(self.loop_energy(ii, jj, &outer)))
            };
            // thal.c would spin forever here; a consistent table never gets here.
            let Some((ni, nj)) = next else { break };
            (i, j) = (ni, nj);
            ps1[i - 1] = j;
            ps2[j - 1] = i;
        }
        (ps1, ps2)
    }
}

/// Outer pair of a bulge/internal loop and its pair-only table terms.
struct Outer {
    ii: usize,
    jj: usize,
    atp: Thermo,
    stack_mm: Thermo,
    tstack: Thermo,
}

/// thal.c rejects non-finite loops, and loops with both ΔH and ΔS positive.
fn forbid_unfavourable(c: Thermo) -> Thermo {
    if !c.is_finite() || (c.h > 0.0 && c.s > 0.0) {
        Thermo::INVALID
    } else {
        c
    }
}
