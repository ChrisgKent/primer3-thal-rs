#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::params::*;
use crate::*;

fn is_positive(x: f64) -> bool {
    x > 0.0
}

fn dbl_eq(x: f64, y: f64) -> i32 {
    if (x - y) < SMALL_NON_ZERO { 1 } else { 2 }
}

fn equal(a: f64, b: f64) -> bool {
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    (a - b).abs() < 1e-5
}

fn bp_indx(a: u8, b: u8) -> i32 {
    BPI[a as usize][b as usize]
}

fn at_penalty_s(a: u8, b: u8, ctx: &Context<'_>) -> f64 {
    ctx.params.atp_s[a as usize][b as usize]
}

fn at_penalty_h(a: u8, b: u8, ctx: &Context<'_>) -> f64 {
    ctx.params.atp_h[a as usize][b as usize]
}

fn salt_correct_s(mv: f64, dv: f64, dntp: f64) -> f64 {
    let mut dntp_val = dntp;
    if dv <= 0.0 {
        dntp_val = dv;
    }
    0.368 * ((mv + 120.0 * (f64::max(0.0, dv - dntp_val)).sqrt()) / 1000.0).ln()
}

fn idx(i: i32, j: i32, len3: i32) -> usize {
    (j + ((i - 1) * len3) - 1) as usize
}

struct Context<'a> {
    params: &'static ThermoParams,
    num_seq1: &'a [u8],
    num_seq2: &'a [u8],
    len1: i32,
    len2: i32,
    len3: i32,
    entropy_dpt: &'a mut [f64],
    enthalpy_dpt: &'a mut [f64],
    dplx_init_s: f64,
    dplx_init_h: f64,
    rc: f64,
}

fn init_matrix(ctx: &mut Context<'_>) {
    let len3 = ctx.len3;
    for i in 1..=ctx.len1 {
        let base = ((i - 1) * len3) as usize;
        for j in 1..=ctx.len2 {
            let id = base + (j - 1) as usize;
            if bp_indx(ctx.num_seq1[i as usize], ctx.num_seq2[j as usize]) == 0 {
                ctx.enthalpy_dpt[id] = f64::INFINITY;
                ctx.entropy_dpt[id] = -1.0;
            } else {
                ctx.enthalpy_dpt[id] = 0.0;
                ctx.entropy_dpt[id] = MIN_ENTROPY;
            }
        }
    }
}

fn symmetry_thermo(seq: &[u8]) -> bool {
    let len = seq.len();
    if len % 2 == 1 {
        return false;
    }
    let mut i = 0usize;
    let mut j = len.saturating_sub(1);
    let half = len / 2;
    while i < half {
        let s = (seq[i] as char).to_ascii_uppercase() as u8;
        let e = (seq[j] as char).to_ascii_uppercase() as u8;
        if (s == b'A' && e != b'T')
            || (s == b'T' && e != b'A')
            || (e == b'A' && s != b'T')
            || (e == b'T' && s != b'A')
        {
            return false;
        }
        if (s == b'C' && e != b'G')
            || (s == b'G' && e != b'C')
            || (e == b'C' && s != b'G')
            || (e == b'G' && s != b'C')
        {
            return false;
        }
        i += 1;
        if j == 0 {
            break;
        }
        j -= 1;
    }
    true
}

impl<'a> Context<'a> {
    fn entropy(&self, i: i32, j: i32) -> f64 {
        self.entropy_dpt[idx(i, j, self.len3)]
    }
    fn enthalpy(&self, i: i32, j: i32) -> f64 {
        self.enthalpy_dpt[idx(i, j, self.len3)]
    }
    fn set_entropy(&mut self, i: i32, j: i32, v: f64) {
        let id = idx(i, j, self.len3);
        self.entropy_dpt[id] = v;
    }
    fn set_enthalpy(&mut self, i: i32, j: i32, v: f64) {
        let id = idx(i, j, self.len3);
        self.enthalpy_dpt[id] = v;
    }
}

/// Hairpin (monomer secondary structure), equivalent to primer3's `thal_hairpin` mode.
pub fn calc_hairpin(
    seq: &str,
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    args.validate()?;
    if seq.is_empty() {
        return Err(ThalError::EmptySequence);
    }
    if seq.len() > THAL_MAX_ALIGN as usize {
        return Err(ThalError::TooLong {
            len: seq.len(),
            max: THAL_MAX_ALIGN as usize,
        });
    }

    let mut result = ThermoResult::default();
    let len1 = seq.len() as i32;

    let params = params();

    let oligo1: Vec<u8> = seq
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();

    let len2 = oligo1.len() as i32;
    let len3 = len2 - 1;

    let mut num_seq1 = vec![4u8; (len1 as usize) + 2];
    let mut num_seq2 = vec![4u8; (len2 as usize) + 2];

    for i in 1..=len1 {
        num_seq1[i as usize] = str2int(oligo1[(i - 1) as usize]);
    }
    for i in 1..=len2 {
        num_seq2[i as usize] = str2int(oligo1[(i - 1) as usize]);
    }

    let salt_correction = salt_correct_s(args.mv, args.dv, args.dntp);

    let mut entropy_dpt = vec![0.0f64; (len1 * len2) as usize];
    let mut enthalpy_dpt = vec![0.0f64; (len1 * len2) as usize];

    let dplx_init_h = 0.0;
    let dplx_init_s = -0.00000000001;
    let rc = 0.0;

    let mut ctx = Context {
        params,
        num_seq1: &num_seq1,
        num_seq2: &num_seq2,
        len1,
        len2,
        len3,
        entropy_dpt: &mut entropy_dpt,
        enthalpy_dpt: &mut enthalpy_dpt,
        dplx_init_s,
        dplx_init_h,
        rc,
    };

    init_matrix2(&mut ctx);
    fill_matrix2(args.max_loop, &mut ctx);

    let mut send5 = vec![0.0f64; (len1 as usize) + 1];
    let mut hend5 = vec![0.0f64; (len1 as usize) + 1];
    calc_terminal_bp(args.temp, &mut send5, &mut hend5, &mut ctx);

    let mh = hend5[len1 as usize];
    let ms = send5[len1 as usize];

    if !mh.is_finite() {
        return Ok(result);
    }

    let mut bp = vec![0i32; len1 as usize];
    tracebacku(&mut bp, args.max_loop, &mut ctx, &send5, &hend5);

    if output_structure {
        let ascii = draw_hairpin(&bp, mh, ms, args.temp, &oligo1, salt_correction);
        result.ascii_structure = Some(ascii);
    }

    let (tm, ds, dh, dg) = calc_hairpin_values(&bp, mh, ms, args.temp, salt_correction);

    result.structure_found = true;
    result.tm = tm;
    result.ds = ds;
    result.dh = dh;
    result.dg = dg;
    Ok(result)
}

#[cfg(feature = "parallel")]
pub fn calc_hairpin_batch<T: AsRef<str> + Sync>(
    seqs: &[T],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    seqs.par_iter()
        .map(|s| calc_hairpin(s.as_ref(), args, output_structure))
        .collect()
}

struct PreparedDimer {
    oligo1: Vec<u8>,
    oligo2: Vec<u8>,
    num_seq1: Vec<u8>,
    num_seq2: Vec<u8>,
    len1: i32,
    len2: i32,
    len3: i32,
}

fn prepare_dimer_sequences(seq1: &str, seq2: &str) -> Result<PreparedDimer, ThalError> {
    if seq1.is_empty() || seq2.is_empty() {
        return Err(ThalError::EmptySequence);
    }
    // As in thal.c: one sequence may be longer than THAL_MAX_ALIGN, not both.
    if seq1.len() > THAL_MAX_ALIGN as usize && seq2.len() > THAL_MAX_ALIGN as usize {
        return Err(ThalError::TooLong {
            len: seq1.len().min(seq2.len()),
            max: THAL_MAX_ALIGN as usize,
        });
    }
    for len in [seq1.len(), seq2.len()] {
        if len > THAL_MAX_SEQ as usize {
            return Err(ThalError::TooLong {
                len,
                max: THAL_MAX_SEQ as usize,
            });
        }
    }
    let len1 = seq1.len() as i32;
    let len2 = seq2.len() as i32;

    let oligo1: Vec<u8> = seq1
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let mut oligo2: Vec<u8> = seq2
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    // Match thalflex: reverse oligo2 for 3'->5' direction in dimer DP.
    oligo2.reverse();

    let mut num_seq1 = vec![4u8; (len1 as usize) + 2];
    let mut num_seq2 = vec![4u8; (len2 as usize) + 2];
    for (i, base) in oligo1.iter().enumerate() {
        num_seq1[i + 1] = str2int(*base);
    }
    for (i, base) in oligo2.iter().enumerate() {
        num_seq2[i + 1] = str2int(*base);
    }

    Ok(PreparedDimer {
        oligo1,
        oligo2,
        num_seq1,
        num_seq2,
        len1,
        len2,
        len3: len2,
    })
}

fn find_best_any(ctx: &Context<'_>) -> Option<(i32, i32)> {
    let mut best_i = 0;
    let mut best_j = 0;
    let mut best_g = f64::INFINITY;
    let mut sh = [0.0f64; 2];
    let temp = TEMP_KELVIN;
    let dplx_h = ctx.dplx_init_h;
    let dplx_s = ctx.dplx_init_s;

    for i in 1..=ctx.len1 {
        for j in 1..=ctx.len2 {
            let enthalpy = ctx.enthalpy(i, j);
            let entropy = ctx.entropy(i, j);
            rsh(i, j, &mut sh, ctx);
            sh[0] += SMALL_NON_ZERO;
            sh[1] += SMALL_NON_ZERO;
            let g1 = (enthalpy + sh[1] + dplx_h) - temp * (entropy + sh[0] + dplx_s);
            if g1 < best_g {
                best_g = g1;
                best_i = i;
                best_j = j;
            }
        }
    }

    if best_g.is_finite() {
        Some((best_i, best_j))
    } else {
        None
    }
}

/// Heterodimer, equivalent to primer3's `thal_any` mode.
pub fn calc_heterodimer(
    seq1: &str,
    seq2: &str,
    args: &ThalArgs,
    output_structure: bool,
) -> Result<ThermoResult, ThalError> {
    args.validate()?;
    let prepared = prepare_dimer_sequences(seq1, seq2)?;
    let mut result = ThermoResult::default();

    let params = params();
    let salt_correction = salt_correct_s(args.mv, args.dv, args.dntp);
    let dplx_init_h = 200.0;
    let dplx_init_s = -5.7;
    let rc = if symmetry_thermo(seq1.as_bytes()) && symmetry_thermo(seq2.as_bytes()) {
        R * (args.dna_conc / 1_000_000_000.0).ln()
    } else {
        R * (args.dna_conc / 4_000_000_000.0).ln()
    };

    let mut entropy_dpt = vec![0.0f64; (prepared.len1 * prepared.len2) as usize];
    let mut enthalpy_dpt = vec![0.0f64; (prepared.len1 * prepared.len2) as usize];

    let mut ctx = Context {
        params,
        num_seq1: &prepared.num_seq1,
        num_seq2: &prepared.num_seq2,
        len1: prepared.len1,
        len2: prepared.len2,
        len3: prepared.len3,
        entropy_dpt: &mut entropy_dpt,
        enthalpy_dpt: &mut enthalpy_dpt,
        dplx_init_s,
        dplx_init_h,
        rc,
    };

    init_matrix(&mut ctx);
    fill_matrix(args.max_loop, &mut ctx);

    let mut sh = [0.0f64; 2];
    let (best_i, best_j) = find_best_any(&ctx).unwrap_or((1, 1));

    rsh(best_i, best_j, &mut sh, &ctx);
    let dh = ctx.enthalpy(best_i, best_j) + sh[1] + ctx.dplx_init_h;
    let ds = ctx.entropy(best_i, best_j) + sh[0] + ctx.dplx_init_s;

    if !ctx.enthalpy(best_i, best_j).is_finite() {
        return Ok(result);
    }

    let mut ps1 = vec![0i32; prepared.len1 as usize];
    let mut ps2 = vec![0i32; prepared.len2 as usize];
    traceback_dimer(best_i, best_j, &mut ps1, &mut ps2, args.max_loop, &mut ctx);

    let (tm, ds_adj, dh_adj, dg) =
        calc_dimer_values(&ps1, &ps2, dh, ds, args.temp, salt_correction, rc);

    result.structure_found = true;
    result.tm = tm;
    result.ds = ds_adj;
    result.dh = dh_adj;
    result.dg = dg;
    if output_structure {
        let ascii = draw_dimer(&ps1, &ps2, &prepared.oligo1, &prepared.oligo2);
        result.ascii_structure = Some(ascii);
    }
    Ok(result)
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
pub fn calc_homodimer_batch<T: AsRef<str> + Sync>(
    seqs: &[T],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    seqs.par_iter()
        .map(|s| calc_homodimer(s.as_ref(), args, output_structure))
        .collect()
}

#[cfg(feature = "parallel")]
pub fn calc_heterodimer_batch(
    pairs: &[(String, String)],
    args: &ThalArgs,
    output_structure: bool,
) -> Vec<Result<ThermoResult, ThalError>> {
    pairs
        .par_iter()
        .map(|(a, b)| calc_heterodimer(a, b, args, output_structure))
        .collect()
}

fn init_matrix2(ctx: &mut Context<'_>) {
    let len3 = ctx.len3;
    for i in 1..=ctx.len1 {
        let base = ((i - 1) * len3) as usize;
        for j in i..=ctx.len2 {
            let id = base + (j - 1) as usize;
            if (j - i) < MIN_HRPN_LOOP + 1
                || bp_indx(ctx.num_seq1[i as usize], ctx.num_seq1[j as usize]) == 0
            {
                ctx.enthalpy_dpt[id] = f64::INFINITY;
                ctx.entropy_dpt[id] = -1.0;
            } else {
                ctx.enthalpy_dpt[id] = 0.0;
                ctx.entropy_dpt[id] = MIN_ENTROPY;
            }
        }
    }
}

fn fill_matrix2(max_loop: i32, ctx: &mut Context<'_>) {
    let mut sh = [0.0f64; 2];
    let len3 = ctx.len3;
    for j in 2..=ctx.len2 {
        let mut i = j - MIN_HRPN_LOOP - 1;
        while i >= 1 {
            let id = ((i - 1) * len3 + (j - 1)) as usize;
            if ctx.enthalpy_dpt[id].is_finite() {
                sh[0] = -1.0;
                sh[1] = f64::INFINITY;
                max_tm2(i, j, ctx);
                cbi(i, j, &mut sh, 0, max_loop, ctx);
                sh[0] = -1.0;
                sh[1] = f64::INFINITY;
                calc_hairpin_inner(i, j, &mut sh, 0, ctx);
                if sh[1].is_finite() {
                    if sh[0] < MIN_ENTROPY_CUTOFF {
                        sh[0] = MIN_ENTROPY;
                        sh[1] = 0.0;
                    }
                    ctx.entropy_dpt[id] = sh[0];
                    ctx.enthalpy_dpt[id] = sh[1];
                }
            }
            i -= 1;
        }
    }
}

fn max_tm(i: i32, j: i32, ctx: &mut Context<'_>) {
    let mut sh = [0.0f64; 2];

    let mut s0 = ctx.entropy(i, j);
    let mut h0 = ctx.enthalpy(i, j);
    rsh(i, j, &mut sh, ctx);
    let t0 = (h0 + ctx.dplx_init_h + sh[1]) / (s0 + ctx.dplx_init_s + sh[0] + ctx.rc);

    let (mut s1, mut h1, t1) =
        if ctx.enthalpy(i - 1, j - 1).is_finite() && hs(i - 1, j - 1, 1, ctx).is_finite() {
            let s1 = ctx.entropy(i - 1, j - 1) + ss(i - 1, j - 1, 1, ctx);
            let h1 = ctx.enthalpy(i - 1, j - 1) + hs(i - 1, j - 1, 1, ctx);
            let t1 = (h1 + ctx.dplx_init_h + sh[1]) / (s1 + ctx.dplx_init_s + sh[0] + ctx.rc);
            (s1, h1, t1)
        } else {
            let s1 = -1.0;
            let h1 = f64::INFINITY;
            let t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);
            (s1, h1, t1)
        };

    if s1 < MIN_ENTROPY_CUTOFF {
        s1 = MIN_ENTROPY;
        h1 = 0.0;
    }
    if s0 < MIN_ENTROPY_CUTOFF {
        s0 = MIN_ENTROPY;
        h0 = 0.0;
    }

    if t1 > t0 {
        ctx.set_entropy(i, j, s1);
        ctx.set_enthalpy(i, j, h1);
    } else if t0 >= t1 {
        ctx.set_entropy(i, j, s0);
        ctx.set_enthalpy(i, j, h0);
    }
}

fn calc_bulge_internal(
    i: i32,
    j: i32,
    ii: i32,
    jj: i32,
    entropy_enthalpy: &mut [f64; 2],
    traceback: i32,
    max_loop: i32,
    ctx: &mut Context<'_>,
) {
    let loop_size1 = ii - i - 1;
    let loop_size2 = jj - j - 1;
    if loop_size1 + loop_size2 > max_loop {
        entropy_enthalpy[0] = -1.0;
        entropy_enthalpy[1] = f64::INFINITY;
        return;
    }
    let loop_size = loop_size1 + loop_size2 - 1;
    let params = ctx.params;

    let mut s;
    let mut h;

    if (loop_size1 == 0 && loop_size2 > 0) || (loop_size2 == 0 && loop_size1 > 0) {
        if loop_size2 == 1 || loop_size1 == 1 {
            h = params.bulge_loop_enthalpies[loop_size as usize]
                + params.stack_enthalpies[ctx.num_seq1[i as usize] as usize]
                    [ctx.num_seq1[ii as usize] as usize][ctx.num_seq2[j as usize] as usize]
                    [ctx.num_seq2[jj as usize] as usize];
            s = params.bulge_loop_entropies[loop_size as usize]
                + params.stack_entropies[ctx.num_seq1[i as usize] as usize]
                    [ctx.num_seq1[ii as usize] as usize][ctx.num_seq2[j as usize] as usize]
                    [ctx.num_seq2[jj as usize] as usize];
            if is_positive(h) || is_positive(s) {
                h = f64::INFINITY;
                s = -1.0;
            }
            h += ctx.enthalpy(i, j);
            s += ctx.entropy(i, j);
            if !h.is_finite() {
                h = f64::INFINITY;
                s = -1.0;
            }
            let mut sh = [0.0f64; 2];
            rsh(ii, jj, &mut sh, ctx);
            let g1 = h + sh[1] - TEMP_KELVIN * (s + sh[0]);
            let g2 = ctx.enthalpy(ii, jj) + sh[1] - TEMP_KELVIN * (ctx.entropy(ii, jj) + sh[0]);
            if (g1 < g2) || (traceback == 1) {
                entropy_enthalpy[0] = s;
                entropy_enthalpy[1] = h;
            }
        } else {
            h = params.bulge_loop_enthalpies[loop_size as usize]
                + at_penalty_h(ctx.num_seq1[i as usize], ctx.num_seq2[j as usize], ctx)
                + at_penalty_h(ctx.num_seq1[ii as usize], ctx.num_seq2[jj as usize], ctx);
            h += ctx.enthalpy(i, j);
            s = params.bulge_loop_entropies[loop_size as usize]
                + at_penalty_s(ctx.num_seq1[i as usize], ctx.num_seq2[j as usize], ctx)
                + at_penalty_s(ctx.num_seq1[ii as usize], ctx.num_seq2[jj as usize], ctx);
            s += ctx.entropy(i, j);
            if !h.is_finite() {
                h = f64::INFINITY;
                s = -1.0;
            }
            if is_positive(h) && is_positive(s) {
                h = f64::INFINITY;
                s = -1.0;
            }
            let mut sh = [0.0f64; 2];
            rsh(ii, jj, &mut sh, ctx);
            let g1 = h + sh[1] - TEMP_KELVIN * (s + sh[0]);
            let g2 = ctx.enthalpy(ii, jj) + sh[1] - TEMP_KELVIN * (ctx.entropy(ii, jj) + sh[0]);
            if (g1 < g2) || (traceback == 1) {
                entropy_enthalpy[0] = s;
                entropy_enthalpy[1] = h;
            }
        }
    } else if loop_size1 == 1 && loop_size2 == 1 {
        s = params.stackint2_entropies[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j + 1) as usize] as usize]
            + params.stackint2_entropies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj - 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize];
        s += ctx.entropy(i, j);
        h = params.stackint2_enthalpies[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j + 1) as usize] as usize]
            + params.stackint2_enthalpies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj - 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize];
        h += ctx.enthalpy(i, j);
        if !h.is_finite() {
            h = f64::INFINITY;
            s = -1.0;
        }
        if is_positive(h) && is_positive(s) {
            h = f64::INFINITY;
            s = -1.0;
        }
        let mut sh = [0.0f64; 2];
        rsh(ii, jj, &mut sh, ctx);
        let g1 = h + sh[1] - TEMP_KELVIN * (s + sh[0]);
        let g2 = ctx.enthalpy(ii, jj) + sh[1] - TEMP_KELVIN * (ctx.entropy(ii, jj) + sh[0]);
        if (g1 < g2) || (traceback == 1) {
            entropy_enthalpy[0] = s;
            entropy_enthalpy[1] = h;
        }
        return;
    } else {
        h = params.interior_loop_enthalpies[loop_size as usize]
            + params.tstack_enthalpies[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j + 1) as usize] as usize]
            + params.tstack_enthalpies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj - 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize]
            + (ILAH * (loop_size1 - loop_size2).abs() as f64);
        h += ctx.enthalpy(i, j);

        s = params.interior_loop_entropies[loop_size as usize]
            + params.tstack_entropies[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j + 1) as usize] as usize]
            + params.tstack_entropies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj - 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize]
            + (ILAS * (loop_size1 - loop_size2).abs() as f64);
        s += ctx.entropy(i, j);
        if !h.is_finite() {
            h = f64::INFINITY;
            s = -1.0;
        }
        if is_positive(h) && is_positive(s) {
            h = f64::INFINITY;
            s = -1.0;
        }
        let mut sh = [0.0f64; 2];
        rsh(ii, jj, &mut sh, ctx);
        let g1 = h + sh[1] - TEMP_KELVIN * (s + sh[0]);
        let g2 = ctx.enthalpy(ii, jj) + sh[1] - TEMP_KELVIN * (ctx.entropy(ii, jj) + sh[0]);
        if (g1 < g2) || (traceback == 1) {
            entropy_enthalpy[0] = s;
            entropy_enthalpy[1] = h;
        }
    }
}

fn fill_matrix(max_loop: i32, ctx: &mut Context<'_>) {
    let mut sh = [0.0f64; 2];
    for i in 1..=ctx.len1 {
        for j in 1..=ctx.len2 {
            if ctx.enthalpy(i, j).is_finite() {
                sh[0] = -1.0;
                sh[1] = f64::INFINITY;
                lsh(i, j, &mut sh, ctx);
                if sh[1].is_finite() {
                    ctx.set_entropy(i, j, sh[0]);
                    ctx.set_enthalpy(i, j, sh[1]);
                }
                if i > 1 && j > 1 {
                    max_tm(i, j, ctx);
                    for d in 3..=max_loop + 2 {
                        let mut ii = i - 1;
                        let mut jj = -ii - d + (j + i);
                        if jj < 1 {
                            ii -= (jj - 1).abs();
                            jj = 1;
                        }
                        while ii > 0 && jj < j {
                            if ctx.enthalpy(ii, jj).is_finite() {
                                sh[0] = -1.0;
                                sh[1] = f64::INFINITY;
                                calc_bulge_internal(ii, jj, i, j, &mut sh, 0, max_loop, ctx);
                                if sh[0] < MIN_ENTROPY_CUTOFF {
                                    sh[0] = MIN_ENTROPY;
                                    sh[1] = 0.0;
                                }
                                if sh[1].is_finite() {
                                    ctx.set_enthalpy(i, j, sh[1]);
                                    ctx.set_entropy(i, j, sh[0]);
                                }
                            }
                            ii -= 1;
                            jj += 1;
                        }
                    }
                }
            }
        }
    }
}

fn traceback_dimer(
    i0: i32,
    j0: i32,
    ps1: &mut [i32],
    ps2: &mut [i32],
    max_loop: i32,
    ctx: &mut Context<'_>,
) {
    let mut i = i0;
    let mut j = j0;
    let mut sh = [0.0f64; 2];
    ps1[(i - 1) as usize] = j;
    ps2[(j - 1) as usize] = i;
    loop {
        sh[0] = -1.0;
        sh[1] = f64::INFINITY;
        lsh(i, j, &mut sh, ctx);
        if equal(ctx.entropy(i, j), sh[0]) && equal(ctx.enthalpy(i, j), sh[1]) {
            break;
        }
        let mut done = false;
        if i > 1
            && j > 1
            && equal(
                ctx.entropy(i, j),
                ss(i - 1, j - 1, 1, ctx) + ctx.entropy(i - 1, j - 1),
            )
            && equal(
                ctx.enthalpy(i, j),
                hs(i - 1, j - 1, 1, ctx) + ctx.enthalpy(i - 1, j - 1),
            )
        {
            i -= 1;
            j -= 1;
            ps1[(i - 1) as usize] = j;
            ps2[(j - 1) as usize] = i;
            done = true;
        }
        for d in 3..=max_loop + 2 {
            if done {
                break;
            }
            let mut ii = i - 1;
            let mut jj = -ii - d + (j + i);
            if jj < 1 {
                ii -= (jj - 1).abs();
                jj = 1;
            }
            while !done && ii > 0 && jj < j {
                sh[0] = -1.0;
                sh[1] = f64::INFINITY;
                calc_bulge_internal(ii, jj, i, j, &mut sh, 1, max_loop, ctx);
                if equal(ctx.entropy(i, j), sh[0]) && equal(ctx.enthalpy(i, j), sh[1]) {
                    i = ii;
                    j = jj;
                    ps1[(i - 1) as usize] = j;
                    ps2[(j - 1) as usize] = i;
                    done = true;
                    break;
                }
                ii -= 1;
                jj += 1;
            }
        }
    }
}
fn max_tm2(i: i32, j: i32, ctx: &mut Context<'_>) {
    let mut s0 = ctx.entropy(i, j);
    let mut h0 = ctx.enthalpy(i, j);
    let t0 = (h0 + ctx.dplx_init_h) / (s0 + ctx.dplx_init_s + ctx.rc);

    let (mut s1, mut h1) = if ctx.enthalpy(i, j).is_finite() {
        (
            ctx.entropy(i + 1, j - 1) + ss(i, j, 2, ctx),
            ctx.enthalpy(i + 1, j - 1) + hs(i, j, 2, ctx),
        )
    } else {
        (-1.0, f64::INFINITY)
    };
    let t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);

    if s1 < MIN_ENTROPY_CUTOFF {
        s1 = MIN_ENTROPY;
        h1 = 0.0;
    }
    if s0 < MIN_ENTROPY_CUTOFF {
        s0 = MIN_ENTROPY;
        h0 = 0.0;
    }
    if t1 > t0 {
        ctx.set_entropy(i, j, s1);
        ctx.set_enthalpy(i, j, h1);
    } else {
        ctx.set_entropy(i, j, s0);
        ctx.set_enthalpy(i, j, h0);
    }
}

fn ss(i: i32, j: i32, k: i32, ctx: &Context<'_>) -> f64 {
    let params = ctx.params;
    if k == 2 {
        if i >= j {
            return -1.0;
        }
        if i == ctx.len1 || j == ctx.len2 + 1 {
            return -1.0;
        }
        let mut ii = i;
        let mut jj = j;
        if ii > ctx.len1 {
            ii -= ctx.len1;
        }
        if jj > ctx.len2 {
            jj -= ctx.len2;
        }
        return params.stack_entropies[ctx.num_seq1[ii as usize] as usize]
            [ctx.num_seq1[(ii + 1) as usize] as usize][ctx.num_seq2[jj as usize] as usize]
            [ctx.num_seq2[(jj - 1) as usize] as usize];
    }
    params.stack_entropies[ctx.num_seq1[i as usize] as usize]
        [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
        [ctx.num_seq2[(j + 1) as usize] as usize]
}

fn hs(i: i32, j: i32, k: i32, ctx: &Context<'_>) -> f64 {
    let params = ctx.params;
    if k == 2 {
        if i >= j {
            return f64::INFINITY;
        }
        if i == ctx.len1 || j == ctx.len2 + 1 {
            return f64::INFINITY;
        }
        let mut ii = i;
        let mut jj = j;
        if ii > ctx.len1 {
            ii -= ctx.len1;
        }
        if jj > ctx.len2 {
            jj -= ctx.len2;
        }
        let check = params.stack_enthalpies[ctx.num_seq1[ii as usize] as usize]
            [ctx.num_seq1[(ii + 1) as usize] as usize][ctx.num_seq2[jj as usize] as usize]
            [ctx.num_seq2[(jj - 1) as usize] as usize];
        if check.is_finite() {
            check
        } else {
            f64::INFINITY
        }
    } else {
        params.stack_enthalpies[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j + 1) as usize] as usize]
    }
}

fn cbi(
    i: i32,
    j: i32,
    entropy_enthalpy: &mut [f64; 2],
    traceback: i32,
    max_loop: i32,
    ctx: &mut Context<'_>,
) {
    let mut d = j - i - 3;
    while d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - max_loop {
        let mut ii = i + 1;
        while ii < j - d && ii <= ctx.len1 {
            let jj = d + ii;
            if traceback == 0 {
                entropy_enthalpy[0] = -1.0;
                entropy_enthalpy[1] = f64::INFINITY;
            }
            if ctx.enthalpy(ii, jj).is_finite() && ctx.enthalpy(i, j).is_finite() {
                calc_bulge_internal2(i, j, ii, jj, entropy_enthalpy, traceback, max_loop, ctx);
                if entropy_enthalpy[1].is_finite() {
                    if entropy_enthalpy[0] < MIN_ENTROPY_CUTOFF {
                        entropy_enthalpy[0] = MIN_ENTROPY;
                        entropy_enthalpy[1] = 0.0;
                    }
                    if traceback == 0 {
                        ctx.set_enthalpy(i, j, entropy_enthalpy[1]);
                        ctx.set_entropy(i, j, entropy_enthalpy[0]);
                    }
                }
            }
            ii += 1;
        }
        d -= 1;
    }
}

fn calc_hairpin_inner(
    i: i32,
    j: i32,
    entropy_enthalpy: &mut [f64; 2],
    traceback: i32,
    ctx: &mut Context<'_>,
) {
    let loop_size = j - i - 1;
    let mut sh = [0.0f64; 2];
    sh[0] = -1.0;
    sh[1] = f64::INFINITY;

    if loop_size < MIN_HRPN_LOOP {
        entropy_enthalpy[0] = -1.0;
        entropy_enthalpy[1] = f64::INFINITY;
        return;
    }
    if i <= ctx.len1 && ctx.len2 < j {
        entropy_enthalpy[0] = -1.0;
        entropy_enthalpy[1] = f64::INFINITY;
        return;
    }

    let params = ctx.params;
    if loop_size <= 30 {
        entropy_enthalpy[1] = params.hairpin_loop_enthalpies[(loop_size - 1) as usize];
        entropy_enthalpy[0] = params.hairpin_loop_entropies[(loop_size - 1) as usize];
    } else {
        entropy_enthalpy[1] = params.hairpin_loop_enthalpies[29];
        entropy_enthalpy[0] = params.hairpin_loop_entropies[29];
    }

    if loop_size > 3 {
        let ni = ctx.num_seq1[i as usize] as usize;
        let ni1 = ctx.num_seq1[(i + 1) as usize] as usize;
        let nj = ctx.num_seq1[j as usize] as usize;
        let nj1 = ctx.num_seq1[(j - 1) as usize] as usize;
        entropy_enthalpy[1] += params.tstack2_enthalpies[ni][ni1][nj][nj1];
        entropy_enthalpy[0] += params.tstack2_entropies[ni][ni1][nj][nj1];
    } else if loop_size == 3 {
        entropy_enthalpy[1] +=
            params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq1[j as usize] as usize];
        entropy_enthalpy[0] +=
            params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq1[j as usize] as usize];
    }

    if loop_size == 3 {
        if params.num_triloops > 0 {
            if let Some(v) = find_triloop(
                &params.triloop_enthalpies,
                params.num_triloops,
                &ctx.num_seq1[i as usize..(i + 5) as usize],
            ) {
                entropy_enthalpy[1] += v;
            }
            if let Some(v) = find_triloop(
                &params.triloop_entropies,
                params.num_triloops,
                &ctx.num_seq1[i as usize..(i + 5) as usize],
            ) {
                entropy_enthalpy[0] += v;
            }
        }
    } else if loop_size == 4 {
        if params.num_tetraloops > 0 {
            if let Some(v) = find_tetraloop(
                &params.tetraloop_enthalpies,
                params.num_tetraloops,
                &ctx.num_seq1[i as usize..(i + 6) as usize],
            ) {
                entropy_enthalpy[1] += v;
            }
            if let Some(v) = find_tetraloop(
                &params.tetraloop_entropies,
                params.num_tetraloops,
                &ctx.num_seq1[i as usize..(i + 6) as usize],
            ) {
                entropy_enthalpy[0] += v;
            }
        }
    }

    if !entropy_enthalpy[1].is_finite() {
        entropy_enthalpy[1] = f64::INFINITY;
        entropy_enthalpy[0] = -1.0;
    }
    if is_positive(entropy_enthalpy[1])
        && is_positive(entropy_enthalpy[0])
        && (!is_positive(ctx.enthalpy(i, j)) || !is_positive(ctx.entropy(i, j)))
    {
        entropy_enthalpy[1] = f64::INFINITY;
        entropy_enthalpy[0] = -1.0;
    }

    rsh(i, j, &mut sh, ctx);
    let g1 = entropy_enthalpy[1] + sh[1] - TEMP_KELVIN * (entropy_enthalpy[0] + sh[0]);
    let g2 = ctx.enthalpy(i, j) + sh[1] - TEMP_KELVIN * (ctx.entropy(i, j) + sh[0]);
    if g2 < g1 && traceback == 0 {
        entropy_enthalpy[0] = ctx.entropy(i, j);
        entropy_enthalpy[1] = ctx.enthalpy(i, j);
    }
}

fn calc_bulge_internal2(
    i: i32,
    j: i32,
    ii: i32,
    jj: i32,
    entropy_enthalpy: &mut [f64; 2],
    traceback: i32,
    max_loop: i32,
    ctx: &mut Context<'_>,
) {
    let loop_size1 = ii - i - 1;
    let loop_size2 = j - jj - 1;
    if loop_size1 + loop_size2 > max_loop {
        entropy_enthalpy[0] = -1.0;
        entropy_enthalpy[1] = f64::INFINITY;
        return;
    }
    let loop_size = loop_size1 + loop_size2 - 1;
    let params = ctx.params;

    if (loop_size1 == 0 && loop_size2 > 0) || (loop_size2 == 0 && loop_size1 > 0) {
        if loop_size2 == 1 || loop_size1 == 1 {
            let mut h = params.bulge_loop_enthalpies[loop_size as usize]
                + params.stack_enthalpies[ctx.num_seq1[i as usize] as usize]
                    [ctx.num_seq1[ii as usize] as usize][ctx.num_seq2[j as usize] as usize]
                    [ctx.num_seq2[jj as usize] as usize];
            let mut s = params.bulge_loop_entropies[loop_size as usize]
                + params.stack_entropies[ctx.num_seq1[i as usize] as usize]
                    [ctx.num_seq1[ii as usize] as usize][ctx.num_seq2[j as usize] as usize]
                    [ctx.num_seq2[jj as usize] as usize];

            if traceback != 1 {
                h += ctx.enthalpy(ii, jj);
                s += ctx.entropy(ii, jj);
            }
            if !h.is_finite() {
                h = f64::INFINITY;
                s = -1.0;
            }
            let t1 = (h + ctx.dplx_init_h) / ((s + ctx.dplx_init_s) + ctx.rc);
            let t2 = (ctx.enthalpy(i, j) + ctx.dplx_init_h)
                / ((ctx.entropy(i, j)) + ctx.dplx_init_s + ctx.rc);
            if (t1 > t2) || ((traceback != 0 && t1 >= t2) || traceback == 1) {
                entropy_enthalpy[0] = s;
                entropy_enthalpy[1] = h;
            }
        } else {
            let mut h = params.bulge_loop_enthalpies[loop_size as usize]
                + params.atp_h[ctx.num_seq1[i as usize] as usize]
                    [ctx.num_seq2[j as usize] as usize]
                + params.atp_h[ctx.num_seq1[ii as usize] as usize]
                    [ctx.num_seq2[jj as usize] as usize];
            if traceback != 1 {
                h += ctx.enthalpy(ii, jj);
            }
            let mut s = params.bulge_loop_entropies[loop_size as usize]
                + params.atp_s[ctx.num_seq1[i as usize] as usize]
                    [ctx.num_seq2[j as usize] as usize]
                + params.atp_s[ctx.num_seq1[ii as usize] as usize]
                    [ctx.num_seq2[jj as usize] as usize];
            if traceback != 1 {
                s += ctx.entropy(ii, jj);
            }
            if !h.is_finite() {
                h = f64::INFINITY;
                s = -1.0;
            }
            let t1 = (h + ctx.dplx_init_h) / ((s + ctx.dplx_init_s) + ctx.rc);
            let t2 = (ctx.enthalpy(i, j) + ctx.dplx_init_h)
                / (ctx.entropy(i, j) + ctx.dplx_init_s + ctx.rc);
            if (t1 > t2) || ((traceback != 0 && t1 >= t2) || (traceback == 1)) {
                entropy_enthalpy[0] = s;
                entropy_enthalpy[1] = h;
            }
        }
    } else if loop_size1 == 1 && loop_size2 == 1 {
        let mut s = params.stackint2_entropies[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j - 1) as usize] as usize]
            + params.stackint2_entropies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj + 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize];
        if traceback != 1 {
            s += ctx.entropy(ii, jj);
        }
        let mut h = params.stackint2_enthalpies[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j - 1) as usize] as usize]
            + params.stackint2_enthalpies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj + 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize];
        if traceback != 1 {
            h += ctx.enthalpy(ii, jj);
        }
        if !h.is_finite() {
            h = f64::INFINITY;
            s = -1.0;
        }
        let t1 = (h + ctx.dplx_init_h) / ((s + ctx.dplx_init_s) + ctx.rc);
        let t2 =
            (ctx.enthalpy(i, j) + ctx.dplx_init_h) / (ctx.entropy(i, j) + ctx.dplx_init_s + ctx.rc);
        if (dbl_eq(t1, t2) == 2) || traceback != 0 {
            if (t1 > t2) || ((traceback != 0 && t1 >= t2) || traceback == 1) {
                entropy_enthalpy[0] = s;
                entropy_enthalpy[1] = h;
            }
        }
        return;
    } else {
        let mut h = params.interior_loop_enthalpies[loop_size as usize]
            + params.tstack_enthalpies[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j - 1) as usize] as usize]
            + params.tstack_enthalpies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj + 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize]
            + (ILAH * (loop_size1 - loop_size2).abs() as f64);
        if traceback != 1 {
            h += ctx.enthalpy(ii, jj);
        }
        let mut s = params.interior_loop_entropies[loop_size as usize]
            + params.tstack_entropies[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j - 1) as usize] as usize]
            + params.tstack_entropies[ctx.num_seq2[jj as usize] as usize]
                [ctx.num_seq2[(jj + 1) as usize] as usize][ctx.num_seq1[ii as usize] as usize]
                [ctx.num_seq1[(ii - 1) as usize] as usize]
            + (ILAS * (loop_size1 - loop_size2).abs() as f64);
        if traceback != 1 {
            s += ctx.entropy(ii, jj);
        }
        if !h.is_finite() {
            h = f64::INFINITY;
            s = -1.0;
        }
        let t1 = (h + ctx.dplx_init_h) / ((s + ctx.dplx_init_s) + ctx.rc);
        let t2 = (ctx.enthalpy(i, j) + ctx.dplx_init_h)
            / ((ctx.entropy(i, j)) + ctx.dplx_init_s + ctx.rc);
        if (t1 > t2) || ((traceback != 0 && t1 >= t2) || (traceback == 1)) {
            entropy_enthalpy[0] = s;
            entropy_enthalpy[1] = h;
        }
    }
}

fn rsh(i: i32, j: i32, entropy_enthalpy: &mut [f64; 2], ctx: &Context<'_>) {
    let params = ctx.params;
    let mut s1;
    let mut h1;
    let mut s2;
    let mut h2;
    let mut t1 = f64::NEG_INFINITY;
    let mut t2;

    if bp_indx(ctx.num_seq1[i as usize], ctx.num_seq2[j as usize]) == 0 {
        entropy_enthalpy[0] = -1.0;
        entropy_enthalpy[1] = f64::INFINITY;
        return;
    }

    s1 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
        + params.tstack2_entropies[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j + 1) as usize] as usize];
    h1 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
        + params.tstack2_enthalpies[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j + 1) as usize] as usize];
    let mut g1 = h1 - TEMP_KELVIN * s1;
    if !h1.is_finite() || g1 > 0.0 {
        h1 = f64::INFINITY;
        s1 = -1.0;
        g1 = 1.0;
    }

    if bp_indx(
        ctx.num_seq1[(i + 1) as usize],
        ctx.num_seq2[(j + 1) as usize],
    ) == 0
        && params.dangle_enthalpies_3[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            .is_finite()
        && params.dangle_enthalpies_5[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq2[j as usize] as usize][ctx.num_seq2[(j + 1) as usize] as usize]
            .is_finite()
    {
        s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_entropies_3[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_entropies_5[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq2[j as usize] as usize][ctx.num_seq2[(j + 1) as usize] as usize];
        h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_enthalpies_3[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_enthalpies_5[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq2[j as usize] as usize][ctx.num_seq2[(j + 1) as usize] as usize];
        let mut g2 = h2 - TEMP_KELVIN * s2;
        if !h2.is_finite() || g2 > 0.0 {
            h2 = f64::INFINITY;
            s2 = -1.0;
            g2 = 1.0;
        }
        t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
        if h1.is_finite() && g1 < 0.0 {
            t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);
            if (t1 < t2) && (g2 < 0.0) {
                s1 = s2;
                h1 = h2;
                t1 = t2;
            }
        } else if g2 < 0.0 {
            s1 = s2;
            h1 = h2;
            t1 = t2;
        }
    } else if bp_indx(
        ctx.num_seq1[(i + 1) as usize],
        ctx.num_seq2[(j + 1) as usize],
    ) == 0
        && params.dangle_enthalpies_3[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize]
            .is_finite()
    {
        s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_entropies_3[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize];
        h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_enthalpies_3[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq2[j as usize] as usize];
        let mut g2 = h2 - TEMP_KELVIN * s2;
        if !h2.is_finite() || g2 > 0.0 {
            h2 = f64::INFINITY;
            s2 = -1.0;
            g2 = 1.0;
        }
        t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
        if h1.is_finite() && g1 < 0.0 {
            t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);
            if (t1 < t2) && (g2 < 0.0) {
                s1 = s2;
                h1 = h2;
                t1 = t2;
            }
        } else if g2 < 0.0 {
            s1 = s2;
            h1 = h2;
            t1 = t2;
        }
    } else if bp_indx(
        ctx.num_seq1[(i + 1) as usize],
        ctx.num_seq2[(j + 1) as usize],
    ) == 0
        && params.dangle_enthalpies_5[ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq2[j as usize] as usize][ctx.num_seq2[(j + 1) as usize] as usize]
            .is_finite()
    {
        s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_entropies_5[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq2[j as usize] as usize][ctx.num_seq2[(j + 1) as usize] as usize];
        h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_enthalpies_5[ctx.num_seq1[i as usize] as usize]
                [ctx.num_seq2[j as usize] as usize][ctx.num_seq2[(j + 1) as usize] as usize];
        let mut g2 = h2 - TEMP_KELVIN * s2;
        if !h2.is_finite() || g2 > 0.0 {
            h2 = f64::INFINITY;
            s2 = -1.0;
            g2 = 1.0;
        }
        t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
        if h1.is_finite() && g1 < 0.0 {
            t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);
            if (t1 < t2) && (g2 < 0.0) {
                s1 = s2;
                h1 = h2;
                t1 = t2;
            }
        } else if g2 < 0.0 {
            s1 = s2;
            h1 = h2;
            t1 = t2;
        }
    }

    s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize];
    h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize];
    t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
    let _g1n = h1 - TEMP_KELVIN * s1;
    let g2n = h2 - TEMP_KELVIN * s2;
    if h1.is_finite() {
        if t1 < t2 {
            entropy_enthalpy[0] = s2;
            entropy_enthalpy[1] = h2;
        } else {
            entropy_enthalpy[0] = s1;
            entropy_enthalpy[1] = h1;
        }
    } else if g2n.is_finite() {
        entropy_enthalpy[0] = s2;
        entropy_enthalpy[1] = h2;
    }
}

fn lsh(i: i32, j: i32, entropy_enthalpy: &mut [f64; 2], ctx: &Context<'_>) {
    let params = ctx.params;
    let mut s1;
    let mut h1;
    let mut s2;
    let mut h2;
    let mut t1 = f64::NEG_INFINITY;
    let mut t2;

    if bp_indx(ctx.num_seq1[i as usize], ctx.num_seq2[j as usize]) == 0 {
        entropy_enthalpy[0] = -1.0;
        entropy_enthalpy[1] = f64::INFINITY;
        return;
    }

    s1 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
        + params.tstack2_entropies[ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i - 1) as usize] as usize];
    h1 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
        + params.tstack2_enthalpies[ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize]
            [ctx.num_seq1[(i - 1) as usize] as usize];
    let mut g1 = h1 - TEMP_KELVIN * s1;
    if !h1.is_finite() || g1 > 0.0 {
        h1 = f64::INFINITY;
        s1 = -1.0;
        g1 = 1.0;
    }

    if bp_indx(
        ctx.num_seq1[(i - 1) as usize],
        ctx.num_seq2[(j - 1) as usize],
    ) != 1
        && params.dangle_enthalpies_3[ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize]
            .is_finite()
        && params.dangle_enthalpies_5[ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq1[i as usize] as usize][ctx.num_seq1[(i - 1) as usize] as usize]
            .is_finite()
    {
        s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_entropies_3[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize]
            + params.dangle_entropies_5[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq1[i as usize] as usize][ctx.num_seq1[(i - 1) as usize] as usize];
        h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_enthalpies_3[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize]
            + params.dangle_enthalpies_5[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq1[i as usize] as usize][ctx.num_seq1[(i - 1) as usize] as usize];
        let mut g2 = h2 - TEMP_KELVIN * s2;
        if !h2.is_finite() || g2 > 0.0 {
            h2 = f64::INFINITY;
            s2 = -1.0;
            g2 = 1.0;
        }
        t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
        if h1.is_finite() && g1 < 0.0 {
            t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);
            if (t1 < t2) && (g2 < 0.0) {
                s1 = s2;
                h1 = h2;
                t1 = t2;
            }
        } else if g2 < 0.0 {
            s1 = s2;
            h1 = h2;
            t1 = t2;
        }
    } else if bp_indx(
        ctx.num_seq1[(i - 1) as usize],
        ctx.num_seq2[(j - 1) as usize],
    ) != 1
        && params.dangle_enthalpies_3[ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize]
            .is_finite()
    {
        s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_entropies_3[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize];
        h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_enthalpies_3[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq2[(j - 1) as usize] as usize][ctx.num_seq1[i as usize] as usize];
        let mut g2 = h2 - TEMP_KELVIN * s2;
        if !h2.is_finite() || g2 > 0.0 {
            h2 = f64::INFINITY;
            s2 = -1.0;
            g2 = 1.0;
        }
        t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
        if h1.is_finite() && g1 < 0.0 {
            t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);
            if (t1 < t2) && (g2 < 0.0) {
                s1 = s2;
                h1 = h2;
                t1 = t2;
            }
        } else if g2 < 0.0 {
            s1 = s2;
            h1 = h2;
            t1 = t2;
        }
    } else if bp_indx(
        ctx.num_seq1[(i - 1) as usize],
        ctx.num_seq2[(j - 1) as usize],
    ) != 1
        && params.dangle_enthalpies_5[ctx.num_seq2[j as usize] as usize]
            [ctx.num_seq1[i as usize] as usize][ctx.num_seq1[(i - 1) as usize] as usize]
            .is_finite()
    {
        s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_entropies_5[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq1[i as usize] as usize][ctx.num_seq1[(i - 1) as usize] as usize];
        h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize]
            + params.dangle_enthalpies_5[ctx.num_seq2[j as usize] as usize]
                [ctx.num_seq1[i as usize] as usize][ctx.num_seq1[(i - 1) as usize] as usize];
        let mut g2 = h2 - TEMP_KELVIN * s2;
        if !h2.is_finite() || g2 > 0.0 {
            h2 = f64::INFINITY;
            s2 = -1.0;
            g2 = 1.0;
        }
        t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
        if h1.is_finite() && g1 < 0.0 {
            t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);
            if (t1 < t2) && (g2 < 0.0) {
                s1 = s2;
                h1 = h2;
                t1 = t2;
            }
        } else if g2 < 0.0 {
            s1 = s2;
            h1 = h2;
            t1 = t2;
        }
    }

    s2 = params.atp_s[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize];
    h2 = params.atp_h[ctx.num_seq1[i as usize] as usize][ctx.num_seq2[j as usize] as usize];
    t2 = (h2 + ctx.dplx_init_h) / (s2 + ctx.dplx_init_s + ctx.rc);
    let _g1n = h1 - TEMP_KELVIN * s1;
    let g2n = h2 - TEMP_KELVIN * s2;
    if h1.is_finite() {
        if t1 < t2 {
            entropy_enthalpy[0] = s2;
            entropy_enthalpy[1] = h2;
        } else {
            entropy_enthalpy[0] = s1;
            entropy_enthalpy[1] = h1;
        }
    } else if g2n.is_finite() {
        entropy_enthalpy[0] = s2;
        entropy_enthalpy[1] = h2;
    }
}

fn calc_terminal_bp(temp: f64, send5: &mut [f64], hend5: &mut [f64], ctx: &mut Context<'_>) {
    send5[0] = -1.0;
    send5[1] = -1.0;
    hend5[0] = f64::INFINITY;
    hend5[1] = f64::INFINITY;
    for i in 2..=ctx.len1 as usize {
        send5[i] = MIN_ENTROPY;
        hend5[i] = 0.0;
    }

    for i in 2..=ctx.len1 {
        let t1 = (hend5[(i - 1) as usize] + ctx.dplx_init_h)
            / (send5[(i - 1) as usize] + ctx.dplx_init_s + ctx.rc);
        let t2 = (end5_1(i, 1, send5, hend5, ctx) + ctx.dplx_init_h)
            / (end5_1(i, 2, send5, hend5, ctx) + ctx.dplx_init_s + ctx.rc);
        let t3 = (end5_2(i, 1, send5, hend5, ctx) + ctx.dplx_init_h)
            / (end5_2(i, 2, send5, hend5, ctx) + ctx.dplx_init_s + ctx.rc);
        let t4 = (end5_3(i, 1, send5, hend5, ctx) + ctx.dplx_init_h)
            / (end5_3(i, 2, send5, hend5, ctx) + ctx.dplx_init_s + ctx.rc);
        let t5 = (end5_4(i, 1, send5, hend5, ctx) + ctx.dplx_init_h)
            / (end5_4(i, 2, send5, hend5, ctx) + ctx.dplx_init_s + ctx.rc);
        let max = max5(t1, t2, t3, t4, t5);
        match max {
            1 => {
                send5[i as usize] = send5[(i - 1) as usize];
                hend5[i as usize] = hend5[(i - 1) as usize];
            }
            2 => {
                let g = end5_1(i, 1, send5, hend5, ctx) - (temp * end5_1(i, 2, send5, hend5, ctx));
                if g < G2 {
                    send5[i as usize] = end5_1(i, 2, send5, hend5, ctx);
                    hend5[i as usize] = end5_1(i, 1, send5, hend5, ctx);
                } else {
                    send5[i as usize] = send5[(i - 1) as usize];
                    hend5[i as usize] = hend5[(i - 1) as usize];
                }
            }
            3 => {
                let g = end5_2(i, 1, send5, hend5, ctx) - (temp * end5_2(i, 2, send5, hend5, ctx));
                if g < G2 {
                    send5[i as usize] = end5_2(i, 2, send5, hend5, ctx);
                    hend5[i as usize] = end5_2(i, 1, send5, hend5, ctx);
                } else {
                    send5[i as usize] = send5[(i - 1) as usize];
                    hend5[i as usize] = hend5[(i - 1) as usize];
                }
            }
            4 => {
                let g = end5_3(i, 1, send5, hend5, ctx) - (temp * end5_3(i, 2, send5, hend5, ctx));
                if g < G2 {
                    send5[i as usize] = end5_3(i, 2, send5, hend5, ctx);
                    hend5[i as usize] = end5_3(i, 1, send5, hend5, ctx);
                } else {
                    send5[i as usize] = send5[(i - 1) as usize];
                    hend5[i as usize] = hend5[(i - 1) as usize];
                }
            }
            5 => {
                let g = end5_4(i, 1, send5, hend5, ctx) - (temp * end5_4(i, 2, send5, hend5, ctx));
                if g < G2 {
                    send5[i as usize] = end5_4(i, 2, send5, hend5, ctx);
                    hend5[i as usize] = end5_4(i, 1, send5, hend5, ctx);
                } else {
                    send5[i as usize] = send5[(i - 1) as usize];
                    hend5[i as usize] = hend5[(i - 1) as usize];
                }
            }
            _ => {}
        }
    }
}

fn end5_1(i: i32, hs: i32, send5: &[f64], hend5: &[f64], ctx: &Context<'_>) -> f64 {
    let mut h_max = f64::INFINITY;
    let mut s_max = -1.0;
    let mut max_tm = f64::NEG_INFINITY;
    for k in 0..=i - MIN_HRPN_LOOP - 2 {
        let t1 =
            (hend5[k as usize] + ctx.dplx_init_h) / (send5[k as usize] + ctx.dplx_init_s + ctx.rc);
        let t2 = (0.0 + ctx.dplx_init_h) / (0.0 + ctx.dplx_init_s + ctx.rc);
        let (mut h, mut s) = if t1 >= t2 {
            (
                hend5[k as usize]
                    + ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                        [ctx.num_seq1[i as usize] as usize]
                    + ctx.enthalpy(k + 1, i),
                send5[k as usize]
                    + ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                        [ctx.num_seq1[i as usize] as usize]
                    + ctx.entropy(k + 1, i),
            )
        } else {
            (
                0.0 + ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                    [ctx.num_seq1[i as usize] as usize]
                    + ctx.enthalpy(k + 1, i),
                0.0 + ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                    [ctx.num_seq1[i as usize] as usize]
                    + ctx.entropy(k + 1, i),
            )
        };
        if !h.is_finite() || h > 0.0 || s > 0.0 {
            h = f64::INFINITY;
            s = -1.0;
        }
        let t1n = (h + ctx.dplx_init_h) / (s + ctx.dplx_init_s + ctx.rc);
        if max_tm < t1n {
            if s > MIN_ENTROPY_CUTOFF {
                h_max = h;
                s_max = s;
                max_tm = t1n;
            }
        }
    }
    if hs == 1 { h_max } else { s_max }
}

fn end5_2(i: i32, hs: i32, send5: &[f64], hend5: &[f64], ctx: &Context<'_>) -> f64 {
    let mut h_max = f64::INFINITY;
    let mut s_max = -1.0;
    let mut max_tm = f64::NEG_INFINITY;
    for k in 0..=i - MIN_HRPN_LOOP - 3 {
        let t1 =
            (hend5[k as usize] + ctx.dplx_init_h) / (send5[k as usize] + ctx.dplx_init_s + ctx.rc);
        let t2 = (0.0 + ctx.dplx_init_h) / (0.0 + ctx.dplx_init_s + ctx.rc);
        let (mut h, mut s) = if t1 >= t2 {
            (
                hend5[k as usize]
                    + ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                        [ctx.num_seq1[i as usize] as usize]
                    + hd5(i, k + 2, ctx)
                    + ctx.enthalpy(k + 2, i),
                send5[k as usize]
                    + ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                        [ctx.num_seq1[i as usize] as usize]
                    + sd5(i, k + 2, ctx)
                    + ctx.entropy(k + 2, i),
            )
        } else {
            (
                0.0 + ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                    [ctx.num_seq1[i as usize] as usize]
                    + hd5(i, k + 2, ctx)
                    + ctx.enthalpy(k + 2, i),
                0.0 + ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                    [ctx.num_seq1[i as usize] as usize]
                    + sd5(i, k + 2, ctx)
                    + ctx.entropy(k + 2, i),
            )
        };
        if !h.is_finite() || h > 0.0 || s > 0.0 {
            h = f64::INFINITY;
            s = -1.0;
        }
        let t1n = (h + ctx.dplx_init_h) / (s + ctx.dplx_init_s + ctx.rc);
        if max_tm < t1n {
            if s > MIN_ENTROPY_CUTOFF {
                h_max = h;
                s_max = s;
                max_tm = t1n;
            }
        }
    }
    if hs == 1 { h_max } else { s_max }
}

fn end5_3(i: i32, hs: i32, send5: &[f64], hend5: &[f64], ctx: &Context<'_>) -> f64 {
    let mut h_max = f64::INFINITY;
    let mut s_max = -1.0;
    let mut max_tm = f64::NEG_INFINITY;
    for k in 0..=i - MIN_HRPN_LOOP - 3 {
        let t1 =
            (hend5[k as usize] + ctx.dplx_init_h) / (send5[k as usize] + ctx.dplx_init_s + ctx.rc);
        let t2 = (0.0 + ctx.dplx_init_h) / (0.0 + ctx.dplx_init_s + ctx.rc);
        let (mut h, mut s) = if t1 >= t2 {
            (
                hend5[k as usize]
                    + ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                        [ctx.num_seq1[(i - 1) as usize] as usize]
                    + hd3(i - 1, k + 1, ctx)
                    + ctx.enthalpy(k + 1, i - 1),
                send5[k as usize]
                    + ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                        [ctx.num_seq1[(i - 1) as usize] as usize]
                    + sd3(i - 1, k + 1, ctx)
                    + ctx.entropy(k + 1, i - 1),
            )
        } else {
            (
                0.0 + ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                    [ctx.num_seq1[(i - 1) as usize] as usize]
                    + hd3(i - 1, k + 1, ctx)
                    + ctx.enthalpy(k + 1, i - 1),
                0.0 + ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                    [ctx.num_seq1[(i - 1) as usize] as usize]
                    + sd3(i - 1, k + 1, ctx)
                    + ctx.entropy(k + 1, i - 1),
            )
        };
        if !h.is_finite() || h > 0.0 || s > 0.0 {
            h = f64::INFINITY;
            s = -1.0;
        }
        let t1n = (h + ctx.dplx_init_h) / (s + ctx.dplx_init_s + ctx.rc);
        if max_tm < t1n {
            if s > MIN_ENTROPY_CUTOFF {
                h_max = h;
                s_max = s;
                max_tm = t1n;
            }
        }
    }
    if hs == 1 { h_max } else { s_max }
}

fn end5_4(i: i32, hs: i32, send5: &[f64], hend5: &[f64], ctx: &Context<'_>) -> f64 {
    let mut h_max = f64::INFINITY;
    let mut s_max = -1.0;
    let mut max_tm = f64::NEG_INFINITY;
    for k in 0..=i - MIN_HRPN_LOOP - 4 {
        let t1 =
            (hend5[k as usize] + ctx.dplx_init_h) / (send5[k as usize] + ctx.dplx_init_s + ctx.rc);
        let t2 = (0.0 + ctx.dplx_init_h) / (0.0 + ctx.dplx_init_s + ctx.rc);
        let (mut h, mut s) = if t1 >= t2 {
            (
                hend5[k as usize]
                    + ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                        [ctx.num_seq1[(i - 1) as usize] as usize]
                    + htstack(i - 1, k + 2, ctx)
                    + ctx.enthalpy(k + 2, i - 1),
                send5[k as usize]
                    + ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                        [ctx.num_seq1[(i - 1) as usize] as usize]
                    + ststack(i - 1, k + 2, ctx)
                    + ctx.entropy(k + 2, i - 1),
            )
        } else {
            (
                0.0 + ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                    [ctx.num_seq1[(i - 1) as usize] as usize]
                    + htstack(i - 1, k + 2, ctx)
                    + ctx.enthalpy(k + 2, i - 1),
                0.0 + ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                    [ctx.num_seq1[(i - 1) as usize] as usize]
                    + ststack(i - 1, k + 2, ctx)
                    + ctx.entropy(k + 2, i - 1),
            )
        };
        if !h.is_finite() || h > 0.0 || s > 0.0 {
            h = f64::INFINITY;
            s = -1.0;
        }
        let t1n = (h + ctx.dplx_init_h) / (s + ctx.dplx_init_s + ctx.rc);
        if max_tm < t1n {
            if s > MIN_ENTROPY_CUTOFF {
                h_max = h;
                s_max = s;
                max_tm = t1n;
            }
        }
    }
    if hs == 1 { h_max } else { s_max }
}

fn sd5(i: i32, j: i32, ctx: &Context<'_>) -> f64 {
    ctx.params.dangle_entropies_5[ctx.num_seq1[i as usize] as usize]
        [ctx.num_seq1[j as usize] as usize][ctx.num_seq1[(j - 1) as usize] as usize]
}

fn hd5(i: i32, j: i32, ctx: &Context<'_>) -> f64 {
    ctx.params.dangle_enthalpies_5[ctx.num_seq1[i as usize] as usize]
        [ctx.num_seq1[j as usize] as usize][ctx.num_seq1[(j - 1) as usize] as usize]
}

fn sd3(i: i32, j: i32, ctx: &Context<'_>) -> f64 {
    ctx.params.dangle_entropies_3[ctx.num_seq1[i as usize] as usize]
        [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq1[j as usize] as usize]
}

fn hd3(i: i32, j: i32, ctx: &Context<'_>) -> f64 {
    ctx.params.dangle_enthalpies_3[ctx.num_seq1[i as usize] as usize]
        [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq1[j as usize] as usize]
}

fn ststack(i: i32, j: i32, ctx: &Context<'_>) -> f64 {
    ctx.params.tstack2_entropies[ctx.num_seq1[i as usize] as usize]
        [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq1[j as usize] as usize]
        [ctx.num_seq1[(j - 1) as usize] as usize]
}

fn htstack(i: i32, j: i32, ctx: &Context<'_>) -> f64 {
    ctx.params.tstack2_enthalpies[ctx.num_seq1[i as usize] as usize]
        [ctx.num_seq1[(i + 1) as usize] as usize][ctx.num_seq1[j as usize] as usize]
        [ctx.num_seq1[(j - 1) as usize] as usize]
}

fn max5(a: f64, b: f64, c: f64, d: f64, e: f64) -> i32 {
    if a > b && a > c && a > d && a > e {
        1
    } else if b > c && b > d && b > e {
        2
    } else if c > d && c > e {
        3
    } else if d > e {
        4
    } else {
        5
    }
}

fn calc_hairpin_values(
    bp: &[i32],
    mh: f64,
    ms: f64,
    temp: f64,
    salt_correction: f64,
) -> (f64, f64, f64, f64) {
    if !mh.is_finite() || !ms.is_finite() {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut n = 0;
    for i in 1..bp.len() {
        if bp[i - 1] > 0 {
            n += 1;
        }
    }
    // Match C integer division for N/2 (truncates for odd N).
    let n_half = (n / 2) as f64;
    let salt_term = (n_half - 1.0) * salt_correction;
    let t = (mh / (ms + salt_term)) - ABSOLUTE_ZERO;
    let mg = mh - (temp * (ms + salt_term));
    let ms_adj = ms + salt_term;
    (t, ms_adj, mh, mg)
}

fn calc_dimer_values(
    ps1: &[i32],
    ps2: &[i32],
    h: f64,
    s: f64,
    temp: f64,
    salt_correction: f64,
    rc: f64,
) -> (f64, f64, f64, f64) {
    if !temp.is_finite() {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut n = 0i32;
    for v in ps1 {
        if *v > 0 {
            n += 1;
        }
    }
    for v in ps2 {
        if *v > 0 {
            n += 1;
        }
    }
    let n = (n / 2) - 1;
    let n_f = n as f64;
    let t = (h / (s + (n_f * salt_correction) + rc)) - ABSOLUTE_ZERO;
    let g = h - (temp * (s + (n_f * salt_correction)));
    let s_adj = s + (n_f * salt_correction);
    (t, s_adj, h, g)
}

fn trim_trailing_whitespace(buf: &mut Vec<u8>) {
    while matches!(buf.last(), Some(b) if b.is_ascii_whitespace()) {
        buf.pop();
    }
}

fn draw_dimer(ps1: &[i32], ps2: &[i32], oligo1: &[u8], oligo2: &[u8]) -> String {
    let len1 = oligo1.len();
    let len2 = oligo2.len();
    let cap = len1 + len2 + 1;

    let mut duplex_0: Vec<u8> = Vec::with_capacity(cap);
    let mut duplex_1: Vec<u8> = Vec::with_capacity(cap);
    let mut duplex_2: Vec<u8> = Vec::with_capacity(cap);
    let mut duplex_3: Vec<u8> = Vec::with_capacity(cap);

    let num_ss1 = ps1.iter().take_while(|v| **v == 0).count();
    let num_ss2 = ps2.iter().take_while(|v| **v == 0).count();

    if num_ss1 >= num_ss2 {
        for i in 0..num_ss1 {
            duplex_0.push(oligo1[i]);
            duplex_1.push(b' ');
            duplex_2.push(b' ');
        }
        for _ in 0..(num_ss1 - num_ss2) {
            duplex_3.push(b' ');
        }
        for j in 0..num_ss2 {
            duplex_3.push(oligo2[j]);
        }
    } else {
        for j in 0..num_ss2 {
            duplex_3.push(oligo2[j]);
            duplex_1.push(b' ');
            duplex_2.push(b' ');
        }
        for _ in 0..(num_ss2 - num_ss1) {
            duplex_0.push(b' ');
        }
        for i in 0..num_ss1 {
            duplex_0.push(oligo1[i]);
        }
    }

    let mut i = num_ss1 + 1; // 1-based
    let mut j = num_ss2 + 1; // 1-based

    while i <= len1 {
        while i <= len1
            && j <= len2
            && ps1[i - 1] != 0
            && ps2[j - 1] != 0
        {
            duplex_0.push(b' ');
            duplex_1.push(oligo1[i - 1]);
            duplex_2.push(oligo2[j - 1]);
            duplex_3.push(b' ');
            i += 1;
            j += 1;
        }

        let mut gap1 = 0usize;
        while i <= len1 && ps1[i - 1] == 0 {
            duplex_0.push(oligo1[i - 1]);
            duplex_1.push(b' ');
            gap1 += 1;
            i += 1;
        }

        let mut gap2 = 0usize;
        while j <= len2 && ps2[j - 1] == 0 {
            duplex_2.push(b' ');
            duplex_3.push(oligo2[j - 1]);
            gap2 += 1;
            j += 1;
        }

        if gap1 < gap2 {
            for _ in 0..(gap2 - gap1) {
                duplex_0.push(b'-');
                duplex_1.push(b' ');
            }
        } else if gap1 > gap2 {
            for _ in 0..(gap1 - gap2) {
                duplex_2.push(b' ');
                duplex_3.push(b'-');
            }
        }
    }

    trim_trailing_whitespace(&mut duplex_0);
    trim_trailing_whitespace(&mut duplex_1);
    trim_trailing_whitespace(&mut duplex_2);
    trim_trailing_whitespace(&mut duplex_3);

    let s0 = String::from_utf8(duplex_0).unwrap_or_default();
    let s1 = String::from_utf8(duplex_1).unwrap_or_default();
    let s2 = String::from_utf8(duplex_2).unwrap_or_default();
    let s3 = String::from_utf8(duplex_3).unwrap_or_default();

    format!("SEQ\t{s0}\nSEQ\t{s1}\nSTR\t{s2}\nSTR\t{s3}\n")
}

fn draw_hairpin(
    bp: &[i32],
    mh: f64,
    ms: f64,
    temp: f64,
    oligo1: &[u8],
    salt_correction: f64,
) -> String {
    let (t, ds, dh, dg) = calc_hairpin_values(bp, mh, ms, temp, salt_correction);
    let mut ascii_row = vec![b'-'; bp.len()];
    for i in 1..=bp.len() {
        let bpv = bp[i - 1];
        if bpv > 0 {
            let idx = (bpv - 1) as usize;
            if bpv as usize > i - 1 {
                ascii_row[idx] = b'\\';
            } else {
                ascii_row[idx] = b'/';
            }
        }
    }

    let mut out = String::new();
    let _ = (t, ds, dh, dg); // keep for parity with C output side effects
    out.push_str("SEQ\t");
    out.push_str(&String::from_utf8_lossy(&ascii_row));
    out.push_str("\nSTR\t");
    out.push_str(&String::from_utf8_lossy(oligo1));
    out.push('\n');
    out
}

#[derive(Clone, Copy)]
struct Tracer {
    i: i32,
    j: i32,
    mtrx: i32,
}

fn tracebacku(bp: &mut [i32], max_loop: i32, ctx: &mut Context<'_>, send5: &[f64], hend5: &[f64]) {
    let mut stack: Vec<Tracer> = Vec::new();
    stack.push(Tracer {
        i: ctx.len1,
        j: 0,
        mtrx: 1,
    });
    let mut sh1 = [0.0f64; 2];
    let mut sh2 = [0.0f64; 2];
    let mut entropy_enthalpy = [0.0f64; 2];

    while let Some(top) = stack.pop() {
        let mut i = top.i;
        let j = top.j;
        if top.mtrx == 1 {
            while i > 0
                && equal(send5[i as usize], send5[(i - 1) as usize])
                && equal(hend5[i as usize], hend5[(i - 1) as usize])
            {
                i -= 1;
            }
            if i == 0 {
                continue;
            }
            if equal(send5[i as usize], end5_1(i, 2, send5, hend5, ctx))
                && equal(hend5[i as usize], end5_1(i, 1, send5, hend5, ctx))
            {
                for k in 0..=i - MIN_HRPN_LOOP - 2 {
                    if equal(
                        send5[i as usize],
                        ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                            [ctx.num_seq1[i as usize] as usize]
                            + ctx.entropy(k + 1, i),
                    ) && equal(
                        hend5[i as usize],
                        ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                            [ctx.num_seq1[i as usize] as usize]
                            + ctx.enthalpy(k + 1, i),
                    ) {
                        stack.push(Tracer {
                            i: k + 1,
                            j: i,
                            mtrx: 0,
                        });
                        break;
                    } else if equal(
                        send5[i as usize],
                        send5[k as usize]
                            + ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                                [ctx.num_seq1[i as usize] as usize]
                            + ctx.entropy(k + 1, i),
                    ) && equal(
                        hend5[i as usize],
                        hend5[k as usize]
                            + ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                                [ctx.num_seq1[i as usize] as usize]
                            + ctx.enthalpy(k + 1, i),
                    ) {
                        stack.push(Tracer {
                            i: k + 1,
                            j: i,
                            mtrx: 0,
                        });
                        stack.push(Tracer {
                            i: k,
                            j: 0,
                            mtrx: 1,
                        });
                        break;
                    }
                }
            } else if equal(send5[i as usize], end5_2(i, 2, send5, hend5, ctx))
                && equal(hend5[i as usize], end5_2(i, 1, send5, hend5, ctx))
            {
                for k in 0..=i - MIN_HRPN_LOOP - 3 {
                    if equal(
                        send5[i as usize],
                        ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                            [ctx.num_seq1[i as usize] as usize]
                            + sd5(i, k + 2, ctx)
                            + ctx.entropy(k + 2, i),
                    ) && equal(
                        hend5[i as usize],
                        ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                            [ctx.num_seq1[i as usize] as usize]
                            + hd5(i, k + 2, ctx)
                            + ctx.enthalpy(k + 2, i),
                    ) {
                        stack.push(Tracer {
                            i: k + 2,
                            j: i,
                            mtrx: 0,
                        });
                        break;
                    } else if equal(
                        send5[i as usize],
                        send5[k as usize]
                            + ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                                [ctx.num_seq1[i as usize] as usize]
                            + sd5(i, k + 2, ctx)
                            + ctx.entropy(k + 2, i),
                    ) && equal(
                        hend5[i as usize],
                        hend5[k as usize]
                            + ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                                [ctx.num_seq1[i as usize] as usize]
                            + hd5(i, k + 2, ctx)
                            + ctx.enthalpy(k + 2, i),
                    ) {
                        stack.push(Tracer {
                            i: k + 2,
                            j: i,
                            mtrx: 0,
                        });
                        stack.push(Tracer {
                            i: k,
                            j: 0,
                            mtrx: 1,
                        });
                        break;
                    }
                }
            } else if equal(send5[i as usize], end5_3(i, 2, send5, hend5, ctx))
                && equal(hend5[i as usize], end5_3(i, 1, send5, hend5, ctx))
            {
                for k in 0..=i - MIN_HRPN_LOOP - 3 {
                    if equal(
                        send5[i as usize],
                        ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                            [ctx.num_seq1[(i - 1) as usize] as usize]
                            + sd3(i - 1, k + 1, ctx)
                            + ctx.entropy(k + 1, i - 1),
                    ) && equal(
                        hend5[i as usize],
                        ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                            [ctx.num_seq1[(i - 1) as usize] as usize]
                            + hd3(i - 1, k + 1, ctx)
                            + ctx.enthalpy(k + 1, i - 1),
                    ) {
                        stack.push(Tracer {
                            i: k + 1,
                            j: i - 1,
                            mtrx: 0,
                        });
                        break;
                    } else if equal(
                        send5[i as usize],
                        send5[k as usize]
                            + ctx.params.atp_s[ctx.num_seq1[(k + 1) as usize] as usize]
                                [ctx.num_seq1[(i - 1) as usize] as usize]
                            + sd3(i - 1, k + 1, ctx)
                            + ctx.entropy(k + 1, i - 1),
                    ) && equal(
                        hend5[i as usize],
                        hend5[k as usize]
                            + ctx.params.atp_h[ctx.num_seq1[(k + 1) as usize] as usize]
                                [ctx.num_seq1[(i - 1) as usize] as usize]
                            + hd3(i - 1, k + 1, ctx)
                            + ctx.enthalpy(k + 1, i - 1),
                    ) {
                        stack.push(Tracer {
                            i: k + 1,
                            j: i - 1,
                            mtrx: 0,
                        });
                        stack.push(Tracer {
                            i: k,
                            j: 0,
                            mtrx: 1,
                        });
                        break;
                    }
                }
            } else if equal(send5[i as usize], end5_4(i, 2, send5, hend5, ctx))
                && equal(hend5[i as usize], end5_4(i, 1, send5, hend5, ctx))
            {
                for k in 0..=i - MIN_HRPN_LOOP - 4 {
                    if equal(
                        send5[i as usize],
                        ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                            [ctx.num_seq1[(i - 1) as usize] as usize]
                            + ststack(i - 1, k + 2, ctx)
                            + ctx.entropy(k + 2, i - 1),
                    ) && equal(
                        hend5[i as usize],
                        ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                            [ctx.num_seq1[(i - 1) as usize] as usize]
                            + htstack(i - 1, k + 2, ctx)
                            + ctx.enthalpy(k + 2, i - 1),
                    ) {
                        stack.push(Tracer {
                            i: k + 2,
                            j: i - 1,
                            mtrx: 0,
                        });
                        break;
                    } else if equal(
                        send5[i as usize],
                        send5[k as usize]
                            + ctx.params.atp_s[ctx.num_seq1[(k + 2) as usize] as usize]
                                [ctx.num_seq1[(i - 1) as usize] as usize]
                            + ststack(i - 1, k + 2, ctx)
                            + ctx.entropy(k + 2, i - 1),
                    ) && equal(
                        hend5[i as usize],
                        hend5[k as usize]
                            + ctx.params.atp_h[ctx.num_seq1[(k + 2) as usize] as usize]
                                [ctx.num_seq1[(i - 1) as usize] as usize]
                            + htstack(i - 1, k + 2, ctx)
                            + ctx.enthalpy(k + 2, i - 1),
                    ) {
                        stack.push(Tracer {
                            i: k + 2,
                            j: i - 1,
                            mtrx: 0,
                        });
                        stack.push(Tracer {
                            i: k,
                            j: 0,
                            mtrx: 1,
                        });
                        break;
                    }
                }
            }
        } else if top.mtrx == 0 {
            bp[(i - 1) as usize] = j;
            bp[(j - 1) as usize] = i;
            sh1[0] = -1.0;
            sh1[1] = f64::INFINITY;
            calc_hairpin_inner(i, j, &mut sh1, 1, ctx);
            sh2[0] = -1.0;
            sh2[1] = f64::INFINITY;
            cbi(i, j, &mut sh2, 2, max_loop, ctx);
            if equal(
                ctx.entropy(i, j),
                ss(i, j, 2, ctx) + ctx.entropy(i + 1, j - 1),
            ) && equal(
                ctx.enthalpy(i, j),
                hs(i, j, 2, ctx) + ctx.enthalpy(i + 1, j - 1),
            ) {
                stack.push(Tracer {
                    i: i + 1,
                    j: j - 1,
                    mtrx: 0,
                });
            } else if equal(ctx.entropy(i, j), sh1[0]) && equal(ctx.enthalpy(i, j), sh1[1]) {
                // hairpin, no further traceback
            } else if equal(ctx.entropy(i, j), sh2[0]) && equal(ctx.enthalpy(i, j), sh2[1]) {
                let mut done = false;
                let mut d = j - i - 3;
                while d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - max_loop && !done {
                    let mut ii = i + 1;
                    while ii < j - d {
                        let jj = d + ii;
                        entropy_enthalpy[0] = -1.0;
                        entropy_enthalpy[1] = f64::INFINITY;
                        calc_bulge_internal2(i, j, ii, jj, &mut entropy_enthalpy, 1, max_loop, ctx);
                        if equal(ctx.entropy(i, j), entropy_enthalpy[0] + ctx.entropy(ii, jj))
                            && equal(
                                ctx.enthalpy(i, j),
                                entropy_enthalpy[1] + ctx.enthalpy(ii, jj),
                            )
                        {
                            stack.push(Tracer {
                                i: ii,
                                j: jj,
                                mtrx: 0,
                            });
                            done = true;
                            break;
                        }
                        ii += 1;
                    }
                    d -= 1;
                }
            }
        }
    }
}
