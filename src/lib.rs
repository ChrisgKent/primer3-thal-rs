#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(dead_code)]

use std::cmp::Ordering;
use std::sync::OnceLock;

use pyo3::prelude::*;
use pyo3::types::PyDict;

const MIN_HRPN_LOOP: i32 = 3;
const MAX_LOOP: i32 = 30;

const R: f64 = 1.9872; // cal/Kmol
const ILAS: f64 = (-300.0 / 310.15);
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

#[derive(Clone, Copy, Debug)]
pub struct ThalArgs {
    pub max_loop: i32,
    pub mv: f64,
    pub dv: f64,
    pub dntp: f64,
    pub dna_conc: f64,
    pub temp: f64, // Kelvin
}

impl ThalArgs {
    pub fn hairpin_defaults() -> Self {
        Self {
            max_loop: MAX_LOOP,
            mv: 50.0,
            dv: 0.0,
            dntp: 0.8,
            dna_conc: 50.0,
            temp: TEMP_KELVIN,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ThermoResult {
    pub structure_found: bool,
    pub tm: f64,
    pub dg: f64,
    pub dh: f64,
    pub ds: f64,
    pub ascii_structure: Option<String>,
}

impl Default for ThermoResult {
    fn default() -> Self {
        Self {
            structure_found: false,
            tm: 0.0,
            dg: 0.0,
            dh: 0.0,
            ds: 0.0,
            ascii_structure: None,
        }
    }
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

static PARAMS: OnceLock<ThermoParams> = OnceLock::new();

fn params() -> &'static ThermoParams {
    PARAMS.get_or_init(|| {
        let mut p = ThermoParams::new();
        load_params(&mut p);
        p
    })
}

const DANGLE_DH: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/dangle.dh");
const DANGLE_DS: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/dangle.ds");
const LOOPS_DH: &str = include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/loops.dh");
const LOOPS_DS: &str = include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/loops.ds");
const STACK_DH: &str = include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/stack.dh");
const STACK_DS: &str = include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/stack.ds");
const STACKMM_DH: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/stackmm.dh");
const STACKMM_DS: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/stackmm.ds");
const TETRALOOP_DH: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/tetraloop.dh");
const TETRALOOP_DS: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/tetraloop.ds");
const TRILOOP_DH: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/triloop.dh");
const TRILOOP_DS: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/triloop.ds");
const TSTACK_TM_INF_DS: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/tstack_tm_inf.ds");
const TSTACK_DH: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/tstack.dh");
const TSTACK2_DH: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/tstack2.dh");
const TSTACK2_DS: &str =
    include_str!("../primer3-py/primer3/src/libprimer3/primer3_config/tstack2.ds");

fn parse_double(token: &str) -> f64 {
    let t = token.trim();
    if t.eq_ignore_ascii_case("inf") {
        f64::INFINITY
    } else {
        t.parse::<f64>().unwrap_or(f64::INFINITY)
    }
}

fn next_line<'a, I: Iterator<Item = &'a str>>(iter: &mut I) -> Option<&'a str> {
    while let Some(line) = iter.next() {
        if !line.trim().is_empty() {
            return Some(line);
        }
    }
    None
}

fn load_params(p: &mut ThermoParams) {
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

fn get_stack(stack_s: &mut [[[[f64; 5]; 5]; 5]; 5], stack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
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

fn get_stackint2(stack_s: &mut [[[[f64; 5]; 5]; 5]; 5], stack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
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

fn get_dangle(
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

fn get_loop(
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

fn get_tstack(tstack_s: &mut [[[[f64; 5]; 5]; 5]; 5], tstack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
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

fn get_tstack2(tstack_s: &mut [[[[f64; 5]; 5]; 5]; 5], tstack_h: &mut [[[[f64; 5]; 5]; 5]; 5]) {
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

fn get_triloop(entropies: &mut [Triloop; 20], enthalpies: &mut [Triloop; 20], num: &mut usize) {
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

fn get_tetraloop(
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

fn table_start_at_s(atp: &mut [[f64; 5]; 5], value: f64) {
    for i in 0..5 {
        for j in 0..5 {
            atp[i][j] = 0.0;
        }
    }
    atp[0][3] = value;
    atp[3][0] = value;
}

fn table_start_at_h(atp: &mut [[f64; 5]; 5], value: f64) {
    for i in 0..5 {
        for j in 0..5 {
            atp[i][j] = 0.0;
        }
    }
    atp[0][3] = value;
    atp[3][0] = value;
}

fn loop_cmp_5(a: &[u8; 5], b: &[u8; 5]) -> Ordering {
    for i in 0..5 {
        if a[i] < b[i] {
            return Ordering::Less;
        } else if a[i] > b[i] {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

fn loop_cmp_6(a: &[u8; 6], b: &[u8; 6]) -> Ordering {
    for i in 0..6 {
        if a[i] < b[i] {
            return Ordering::Less;
        } else if a[i] > b[i] {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

fn to_loop_vals(s: &str, len: usize) -> [u8; 5] {
    let mut out = [4u8; 5];
    for (i, b) in s.as_bytes().iter().take(len).enumerate() {
        out[i] = str2int(*b);
    }
    out
}

fn to_loop_vals_6(s: &str) -> [u8; 6] {
    let mut out = [4u8; 6];
    for (i, b) in s.as_bytes().iter().take(6).enumerate() {
        out[i] = str2int(*b);
    }
    out
}

fn str2int(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' | b'0' => 0,
        b'C' | b'1' => 1,
        b'G' | b'2' => 2,
        b'T' | b'3' => 3,
        _ => 4,
    }
}

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

pub fn calc_hairpin(seq: &str, args: &ThalArgs, output_structure: bool) -> ThermoResult {
    let mut result = ThermoResult::default();

    if seq.is_empty() {
        return result;
    }

    let len1 = seq.len() as i32;
    if len1 > THAL_MAX_SEQ {
        return result;
    }
    if len1 > THAL_MAX_ALIGN {
        return result;
    }

    let params = params();

    let mut oligo1: Vec<u8> = seq
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let mut oligo2 = oligo1.clone();

    let len2 = oligo2.len() as i32;
    let len3 = len2 - 1;

    let mut num_seq1 = vec![4u8; (len1 as usize) + 2];
    let mut num_seq2 = vec![4u8; (len2 as usize) + 2];

    for i in 1..=len1 {
        num_seq1[i as usize] = str2int(oligo1[(i - 1) as usize]);
    }
    for i in 1..=len2 {
        num_seq2[i as usize] = str2int(oligo2[(i - 1) as usize]);
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
        result.structure_found = false;
        result.tm = 0.0;
        return result;
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
    result
}

fn init_matrix2(ctx: &mut Context<'_>) {
    for i in 1..=ctx.len1 {
        for j in i..=ctx.len2 {
            if (j - i) < MIN_HRPN_LOOP + 1
                || bp_indx(ctx.num_seq1[i as usize], ctx.num_seq1[j as usize]) == 0
            {
                ctx.set_enthalpy(i, j, f64::INFINITY);
                ctx.set_entropy(i, j, -1.0);
            } else {
                ctx.set_enthalpy(i, j, 0.0);
                ctx.set_entropy(i, j, MIN_ENTROPY);
            }
        }
    }
}

fn fill_matrix2(max_loop: i32, ctx: &mut Context<'_>) {
    let mut sh = [0.0f64; 2];
    for j in 2..=ctx.len2 {
        let mut i = j - MIN_HRPN_LOOP - 1;
        while i >= 1 {
            if ctx.enthalpy(i, j).is_finite() {
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
                    ctx.set_entropy(i, j, sh[0]);
                    ctx.set_enthalpy(i, j, sh[1]);
                }
            }
            i -= 1;
        }
    }
}

fn max_tm2(i: i32, j: i32, ctx: &mut Context<'_>) {
    let mut t0 = f64::NEG_INFINITY;
    let mut t1 = f64::NEG_INFINITY;
    let mut s0 = ctx.entropy(i, j);
    let mut h0 = ctx.enthalpy(i, j);
    t0 = (h0 + ctx.dplx_init_h) / (s0 + ctx.dplx_init_s + ctx.rc);

    let (mut s1, mut h1) = if ctx.enthalpy(i, j).is_finite() {
        (
            ctx.entropy(i + 1, j - 1) + ss(i, j, 2, ctx),
            ctx.enthalpy(i + 1, j - 1) + hs(i, j, 2, ctx),
        )
    } else {
        (-1.0, f64::INFINITY)
    };
    t1 = (h1 + ctx.dplx_init_h) / (s1 + ctx.dplx_init_s + ctx.rc);

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
    let mut s1 = -1.0;
    let mut h1 = f64::INFINITY;
    let mut s2 = -1.0;
    let mut h2 = f64::INFINITY;
    let mut t1 = f64::NEG_INFINITY;
    let mut t2 = f64::NEG_INFINITY;

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
    let g1n = h1 - TEMP_KELVIN * s1;
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
    let mut s1 = -1.0;
    let mut h1 = f64::INFINITY;
    let mut s2 = -1.0;
    let mut h2 = f64::INFINITY;
    let mut t1 = f64::NEG_INFINITY;
    let mut t2 = f64::NEG_INFINITY;

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
    let g1n = h1 - TEMP_KELVIN * s1;
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
        let mut t1 = (hend5[(i - 1) as usize] + ctx.dplx_init_h)
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
        t1 = (hend5[(i - 1) as usize] + ctx.dplx_init_h)
            / (send5[(i - 1) as usize] + ctx.dplx_init_s + ctx.rc);
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
    let t = (mh / (ms + (((n as f64 / 2.0) - 1.0) * salt_correction))) - ABSOLUTE_ZERO;
    let mg = mh - (temp * (ms + (((n as f64 / 2.0) - 1.0) * salt_correction)));
    let ms_adj = ms + (((n as f64 / 2.0) - 1.0) * salt_correction);
    (t, ms_adj, mh, mg)
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

fn find_triloop(entries: &[Triloop], num: usize, key: &[u8]) -> Option<f64> {
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

fn find_tetraloop(entries: &[Tetraloop], num: usize, key: &[u8]) -> Option<f64> {
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

#[pyfunction]
fn calc_hairpin_py(
    py: Python<'_>,
    seq: &str,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
) -> PyResult<Py<PyDict>> {
    let args = ThalArgs {
        max_loop,
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        temp: temp_c + 273.15,
    };
    let res = calc_hairpin(seq, &args, output_structure);
    let out = PyDict::new(py);
    out.set_item("structure_found", res.structure_found)?;
    out.set_item("tm", res.tm)?;
    out.set_item("ds", res.ds)?;
    out.set_item("dh", res.dh)?;
    out.set_item("dg", res.dg)?;
    match res.ascii_structure {
        Some(s) => out.set_item("ascii_structure", s)?,
        None => out.set_item("ascii_structure", py.None())?,
    };
    Ok(out.unbind())
}

#[pymodule]
fn primer3_py_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_hairpin_py, m)?)?;
    m.add("__doc__", "Rust hairpin thermodynamics (PyO3 bindings)")?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hairpin_matches_readme_example() {
        let args = ThalArgs::hairpin_defaults();
        let res = calc_hairpin("CCCCCATCCGATCAGGGGG", &args, false);
        assert!(res.structure_found);
        assert!((res.tm - 34.15).abs() < 0.05);
        assert!((res.dg - 337.09).abs() < 0.5);
        assert!((res.dh - (-36300.00)).abs() < 1.0);
        assert!((res.ds - (-118.13)).abs() < 0.1);
    }
}
