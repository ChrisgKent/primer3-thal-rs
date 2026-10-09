use std::cmp::Ordering;
use std::sync::OnceLock;

use crate::*;

pub(crate) static PARAMS: OnceLock<ThermoParams> = OnceLock::new();

pub(crate) fn params() -> &'static ThermoParams {
    PARAMS.get_or_init(|| {
        let mut p = ThermoParams::new();
        load_params(&mut p);
        p
    })
}

const DANGLE_DH: &str =
    include_str!("../primer3_config/dangle.dh");
const DANGLE_DS: &str =
    include_str!("../primer3_config/dangle.ds");
const LOOPS_DH: &str = include_str!("../primer3_config/loops.dh");
const LOOPS_DS: &str = include_str!("../primer3_config/loops.ds");
const STACK_DH: &str = include_str!("../primer3_config/stack.dh");
const STACK_DS: &str = include_str!("../primer3_config/stack.ds");
const STACKMM_DH: &str =
    include_str!("../primer3_config/stackmm.dh");
const STACKMM_DS: &str =
    include_str!("../primer3_config/stackmm.ds");
const TETRALOOP_DH: &str =
    include_str!("../primer3_config/tetraloop.dh");
const TETRALOOP_DS: &str =
    include_str!("../primer3_config/tetraloop.ds");
const TRILOOP_DH: &str =
    include_str!("../primer3_config/triloop.dh");
const TRILOOP_DS: &str =
    include_str!("../primer3_config/triloop.ds");
const TSTACK_TM_INF_DS: &str =
    include_str!("../primer3_config/tstack_tm_inf.ds");
const TSTACK_DH: &str =
    include_str!("../primer3_config/tstack.dh");
const TSTACK2_DH: &str =
    include_str!("../primer3_config/tstack2.dh");
const TSTACK2_DS: &str =
    include_str!("../primer3_config/tstack2.ds");

pub(crate) fn parse_double(token: &str) -> f64 {
    let t = token.trim();
    if t.eq_ignore_ascii_case("inf") {
        f64::INFINITY
    } else {
        t.parse::<f64>().unwrap_or(f64::INFINITY)
    }
}

pub(crate) fn next_line<'a, I: Iterator<Item = &'a str>>(iter: &mut I) -> Option<&'a str> {
    while let Some(line) = iter.next() {
        if !line.trim().is_empty() {
            return Some(line);
        }
    }
    None
}

pub(crate) fn load_params(p: &mut ThermoParams) {
    get_stack(&mut p.stack_entropies, &mut p.stack_enthalpies);
    get_stackint2(&mut p.stackint2_entropies, &mut p.stackint2_enthalpies);
    get_dangle(
        &mut p.dangle_entropies_3,
        &mut p.dangle_enthalpies_3,
        &mut p.dangle_entropies_5,
        &mut p.dangle_enthalpies_5,
    );
    get_loop(
        &mut p.hairpin_loop_entropies,
        &mut p.interior_loop_entropies,
        &mut p.bulge_loop_entropies,
        &mut p.hairpin_loop_enthalpies,
        &mut p.interior_loop_enthalpies,
        &mut p.bulge_loop_enthalpies,
    );
    get_tstack(&mut p.tstack_entropies, &mut p.tstack_enthalpies);
    get_tstack2(&mut p.tstack2_entropies, &mut p.tstack2_enthalpies);
    get_triloop(
        &mut p.triloop_entropies,
        &mut p.triloop_enthalpies,
        &mut p.num_triloops,
    );
    get_tetraloop(
        &mut p.tetraloop_entropies,
        &mut p.tetraloop_enthalpies,
        &mut p.num_tetraloops,
    );
    table_start_at_s(&mut p.atp_s, AT_S);
    table_start_at_h(&mut p.atp_h, AT_H);
}

pub(crate) fn get_stack(stack_s: &mut [[[[f64; 5]; 5]; 5]; 5], stack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
    let mut ds = STACK_DS.lines();
    let mut dh = STACK_DH.lines();
    for i in 0..5 {
        for ii in 0..5 {
            for j in 0..5 {
                for jj in 0..5 {
                    if i == 4 || j == 4 || ii == 4 || jj == 4 {
                        stack_s[i][ii][j][jj] = -1.0;
                        stack_h[i][ii][j][jj] = f64::INFINITY;
                    } else {
                        let s = parse_double(next_line(&mut ds).unwrap_or("inf"));
                        let h = parse_double(next_line(&mut dh).unwrap_or("inf"));
                        if !s.is_finite() || !h.is_finite() {
                            stack_s[i][ii][j][jj] = -1.0;
                            stack_h[i][ii][j][jj] = f64::INFINITY;
                        } else {
                            stack_s[i][ii][j][jj] = s;
                            stack_h[i][ii][j][jj] = h;
                        }
                    }
                }
            }
        }
    }
}

pub(crate) fn get_stackint2(stack_s: &mut [[[[f64; 5]; 5]; 5]; 5], stack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
    let mut ds = STACKMM_DS.lines();
    let mut dh = STACKMM_DH.lines();
    for i in 0..5 {
        for ii in 0..5 {
            for j in 0..5 {
                for jj in 0..5 {
                    if i == 4 || j == 4 || ii == 4 || jj == 4 {
                        stack_s[i][ii][j][jj] = -1.0;
                        stack_h[i][ii][j][jj] = f64::INFINITY;
                    } else {
                        let s = parse_double(next_line(&mut ds).unwrap_or("inf"));
                        let h = parse_double(next_line(&mut dh).unwrap_or("inf"));
                        if !s.is_finite() || !h.is_finite() {
                            stack_s[i][ii][j][jj] = -1.0;
                            stack_h[i][ii][j][jj] = f64::INFINITY;
                        } else {
                            stack_s[i][ii][j][jj] = s;
                            stack_h[i][ii][j][jj] = h;
                        }
                    }
                }
            }
        }
    }
}

pub(crate) fn get_dangle(
    s3: &mut [[[f64; 5]; 5]; 5],
    h3: &mut [[[f64; 5]; 5]; 5],
    s5: &mut [[[f64; 5]; 5]; 5],
    h5: &mut [[[f64; 5]; 5]; 5],
) {
    let mut ds = DANGLE_DS.lines();
    let mut dh = DANGLE_DH.lines();
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                if i == 4 || j == 4 || k == 4 {
                    s3[i][k][j] = -1.0;
                    h3[i][k][j] = f64::INFINITY;
                } else {
                    let s = parse_double(next_line(&mut ds).unwrap_or("inf"));
                    let h = parse_double(next_line(&mut dh).unwrap_or("inf"));
                    if !s.is_finite() || !h.is_finite() {
                        s3[i][k][j] = -1.0;
                        h3[i][k][j] = f64::INFINITY;
                    } else {
                        s3[i][k][j] = s;
                        h3[i][k][j] = h;
                    }
                }
            }
        }
    }
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                if i == 4 || j == 4 || k == 4 {
                    s5[i][j][k] = -1.0;
                    h5[i][j][k] = f64::INFINITY;
                } else {
                    let s = parse_double(next_line(&mut ds).unwrap_or("inf"));
                    let h = parse_double(next_line(&mut dh).unwrap_or("inf"));
                    if !s.is_finite() || !h.is_finite() {
                        s5[i][j][k] = -1.0;
                        h5[i][j][k] = f64::INFINITY;
                    } else {
                        s5[i][j][k] = s;
                        h5[i][j][k] = h;
                    }
                }
            }
        }
    }
}

pub(crate) fn get_loop(
    hairpin_s: &mut [f64; 30],
    interior_s: &mut [f64; 30],
    bulge_s: &mut [f64; 30],
    hairpin_h: &mut [f64; 30],
    interior_h: &mut [f64; 30],
    bulge_h: &mut [f64; 30],
) {
    let mut ds = LOOPS_DS.lines();
    let mut dh = LOOPS_DH.lines();
    for k in 0..30 {
        let line_s = next_line(&mut ds).unwrap_or("0 0 0 0");
        let parts_s: Vec<&str> = line_s.split_whitespace().collect();
        interior_s[k] = parse_double(parts_s.get(1).unwrap_or(&"inf"));
        bulge_s[k] = parse_double(parts_s.get(2).unwrap_or(&"inf"));
        hairpin_s[k] = parse_double(parts_s.get(3).unwrap_or(&"inf"));

        let line_h = next_line(&mut dh).unwrap_or("0 0 0 0");
        let parts_h: Vec<&str> = line_h.split_whitespace().collect();
        interior_h[k] = parse_double(parts_h.get(1).unwrap_or(&"inf"));
        bulge_h[k] = parse_double(parts_h.get(2).unwrap_or(&"inf"));
        hairpin_h[k] = parse_double(parts_h.get(3).unwrap_or(&"inf"));
    }
}

pub(crate) fn get_tstack(tstack_s: &mut [[[[f64; 5]; 5]; 5]; 5], tstack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
    let mut ds = TSTACK_TM_INF_DS.lines();
    let mut dh = TSTACK_DH.lines();
    for i1 in 0..5 {
        for i2 in 0..5 {
            for j1 in 0..5 {
                for j2 in 0..5 {
                    if i1 == 4 || j1 == 4 {
                        tstack_h[i1][i2][j1][j2] = f64::INFINITY;
                        tstack_s[i1][i2][j1][j2] = -1.0;
                    } else if i2 == 4 || j2 == 4 {
                        tstack_s[i1][i2][j1][j2] = 0.00000000001;
                        tstack_h[i1][i2][j1][j2] = 0.0;
                    } else {
                        let s = parse_double(next_line(&mut ds).unwrap_or("inf"));
                        let h = parse_double(next_line(&mut dh).unwrap_or("inf"));
                        if !s.is_finite() || !h.is_finite() {
                            tstack_s[i1][i2][j1][j2] = -1.0;
                            tstack_h[i1][i2][j1][j2] = f64::INFINITY;
                        } else {
                            tstack_s[i1][i2][j1][j2] = s;
                            tstack_h[i1][i2][j1][j2] = h;
                        }
                    }
                }
            }
        }
    }
}

pub(crate) fn get_tstack2(tstack_s: &mut [[[[f64; 5]; 5]; 5]; 5], tstack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
    let mut ds = TSTACK2_DS.lines();
    let mut dh = TSTACK2_DH.lines();
    for i1 in 0..5 {
        for i2 in 0..5 {
            for j1 in 0..5 {
                for j2 in 0..5 {
                    if i1 == 4 || j1 == 4 {
                        tstack_h[i1][i2][j1][j2] = f64::INFINITY;
                        tstack_s[i1][i2][j1][j2] = -1.0;
                    } else if i2 == 4 || j2 == 4 {
                        tstack_s[i1][i2][j1][j2] = 0.00000000001;
                        tstack_h[i1][i2][j1][j2] = 0.0;
                    } else {
                        let s = parse_double(next_line(&mut ds).unwrap_or("inf"));
                        let h = parse_double(next_line(&mut dh).unwrap_or("inf"));
                        if !s.is_finite() || !h.is_finite() {
                            tstack_s[i1][i2][j1][j2] = -1.0;
                            tstack_h[i1][i2][j1][j2] = f64::INFINITY;
                        } else {
                            tstack_s[i1][i2][j1][j2] = s;
                            tstack_h[i1][i2][j1][j2] = h;
                        }
                    }
                }
            }
        }
    }
}

pub(crate) fn get_triloop(entropies: &mut [Triloop; 20], enthalpies: &mut [Triloop; 20], num: &mut usize) {
    let mut ds = TRILOOP_DS.lines();
    let mut count = 0usize;
    while let Some(line) = next_line(&mut ds) {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        if count >= entropies.len() {
            break;
        }
        let loop_vals = to_loop_vals(parts[0], 5);
        entropies[count] = Triloop {
            loop_vals,
            value: parse_double(parts[1]),
        };
        count += 1;
    }
    entropies[..count].sort_by(|a, b| loop_cmp_5(&a.loop_vals, &b.loop_vals));

    let mut dh = TRILOOP_DH.lines();
    let mut count_h = 0usize;
    while let Some(line) = next_line(&mut dh) {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        if count_h >= enthalpies.len() {
            break;
        }
        let loop_vals = to_loop_vals(parts[0], 5);
        enthalpies[count_h] = Triloop {
            loop_vals,
            value: parse_double(parts[1]),
        };
        count_h += 1;
    }
    enthalpies[..count_h].sort_by(|a, b| loop_cmp_5(&a.loop_vals, &b.loop_vals));
    *num = count_h;
}

pub(crate) fn get_tetraloop(
    entropies: &mut [Tetraloop; 146],
    enthalpies: &mut [Tetraloop; 146],
    num: &mut usize,
) {
    let mut ds = TETRALOOP_DS.lines();
    let mut count = 0usize;
    while let Some(line) = next_line(&mut ds) {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        if count >= entropies.len() {
            break;
        }
        let loop_vals = to_loop_vals_6(parts[0]);
        entropies[count] = Tetraloop {
            loop_vals,
            value: parse_double(parts[1]),
        };
        count += 1;
    }
    entropies[..count].sort_by(|a, b| loop_cmp_6(&a.loop_vals, &b.loop_vals));

    let mut dh = TETRALOOP_DH.lines();
    let mut count_h = 0usize;
    while let Some(line) = next_line(&mut dh) {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        if count_h >= enthalpies.len() {
            break;
        }
        let loop_vals = to_loop_vals_6(parts[0]);
        enthalpies[count_h] = Tetraloop {
            loop_vals,
            value: parse_double(parts[1]),
        };
        count_h += 1;
    }
    enthalpies[..count_h].sort_by(|a, b| loop_cmp_6(&a.loop_vals, &b.loop_vals));
    *num = count_h;
}

pub(crate) fn table_start_at_s(atp: &mut [[f64; 5]; 5], value: f64) {
    for i in 0..5 {
        for j in 0..5 {
            atp[i][j] = 0.00000000001;
        }
    }
    atp[0][3] = value;
    atp[3][0] = value;
}

pub(crate) fn table_start_at_h(atp: &mut [[f64; 5]; 5], value: f64) {
    for i in 0..5 {
        for j in 0..5 {
            atp[i][j] = 0.0;
        }
    }
    atp[0][3] = value;
    atp[3][0] = value;
}

pub(crate) fn loop_cmp_5(a: &[u8; 5], b: &[u8; 5]) -> Ordering {
    for i in 0..5 {
        if a[i] < b[i] {
            return Ordering::Less;
        } else if a[i] > b[i] {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

pub(crate) fn loop_cmp_6(a: &[u8; 6], b: &[u8; 6]) -> Ordering {
    for i in 0..6 {
        if a[i] < b[i] {
            return Ordering::Less;
        } else if a[i] > b[i] {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

pub(crate) fn to_loop_vals(s: &str, len: usize) -> [u8; 5] {
    let mut out = [4u8; 5];
    for (i, b) in s.as_bytes().iter().take(len).enumerate() {
        out[i] = str2int(*b);
    }
    out
}

pub(crate) fn to_loop_vals_6(s: &str) -> [u8; 6] {
    let mut out = [4u8; 6];
    for (i, b) in s.as_bytes().iter().take(6).enumerate() {
        out[i] = str2int(*b);
    }
    out
}

pub(crate) const fn build_str2int_lut() -> [u8; 256] {
    let mut lut = [4u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'0' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'c' as usize] = 1;
    lut[b'1' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'g' as usize] = 2;
    lut[b'2' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b't' as usize] = 3;
    lut[b'3' as usize] = 3;
    lut
}

const STR2INT_LUT: [u8; 256] = build_str2int_lut();

pub(crate) fn str2int(c: u8) -> u8 {
    STR2INT_LUT[c as usize]
}

pub(crate) fn find_triloop(entries: &[Triloop], num: usize, key: &[u8]) -> Option<f64> {
    if key.len() < 5 || num == 0 {
        return None;
    }
    let mut k = [0u8; 5];
    k.copy_from_slice(&key[..5]);
    let slice = &entries[..num];
    match slice.binary_search_by(|e| loop_cmp_5(&e.loop_vals, &k)) {
        Ok(idx) => Some(slice[idx].value),
        Err(_) => None,
    }
}

pub(crate) fn find_tetraloop(entries: &[Tetraloop], num: usize, key: &[u8]) -> Option<f64> {
    if key.len() < 6 || num == 0 {
        return None;
    }
    let mut k = [0u8; 6];
    k.copy_from_slice(&key[..6]);
    let slice = &entries[..num];
    match slice.binary_search_by(|e| loop_cmp_6(&e.loop_vals, &k)) {
        Ok(idx) => Some(slice[idx].value),
        Err(_) => None,
    }
}
