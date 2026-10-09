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

#[derive(Clone, Debug)]
pub struct DimerResult {
    pub structure_found: bool,
    pub tm: f64,
    pub dg: f64,
    pub dh: f64,
    pub ds: f64,
    pub ascii_structure: Option<String>,
}

impl Default for DimerResult {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hairpin_matches_readme_example() {
        let args = ThalArgs::hairpin_defaults();
        let res = calc_hairpin("CCCCCATCCGATCAGGGGG", &args, false);
        assert!(res.structure_found);
        assert!((res.tm - 62.63).abs() < 0.05);
        assert!((res.dg - (-2770.52)).abs() < 0.5);
        assert!((res.dh - (-36300.00)).abs() < 1.0);
        assert!((res.ds - (-108.11)).abs() < 0.1);
    }
}
