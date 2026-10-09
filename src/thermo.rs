//! Entropy/enthalpy pairs: every thal energy term and DP cell is one of these.

use std::ops::Add;

use crate::{MIN_ENTROPY, MIN_ENTROPY_CUTOFF, TEMP_KELVIN, TRACE_EPS};

/// Entropy `s` (cal/K/mol) and enthalpy `h` (cal/mol).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Thermo {
    pub s: f64,
    pub h: f64,
}

/// Duplex initiation terms and the strand-concentration correction `rc`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Init {
    pub s: f64,
    pub h: f64,
    pub rc: f64,
}

impl Thermo {
    /// primer3's marker for "no structure" / forbidden.
    pub const INVALID: Thermo = Thermo {
        s: -1.0,
        h: f64::INFINITY,
    };
    pub const ZERO: Thermo = Thermo { s: 0.0, h: 0.0 };
    /// Starting value of a cell that can pair.
    pub const PAIRABLE: Thermo = Thermo {
        s: MIN_ENTROPY,
        h: 0.0,
    };

    pub const fn new(s: f64, h: f64) -> Self {
        Self { s, h }
    }

    pub fn is_finite(self) -> bool {
        self.h.is_finite()
    }

    /// Melting temperature (K) including duplex initiation.
    pub fn tm(self, init: &Init) -> f64 {
        (self.h + init.h) / (self.s + init.s + init.rc)
    }

    /// Free energy at primer3's fixed ranking temperature (37 °C).
    pub fn dg37(self) -> f64 {
        self.h - TEMP_KELVIN * self.s
    }

    /// primer3 replaces implausibly low entropies by the pairable start value.
    pub fn clamp_min_entropy(self) -> Self {
        if self.s < MIN_ENTROPY_CUTOFF {
            Self::PAIRABLE
        } else {
            self
        }
    }

    /// The tolerance-based comparison thal.c uses during traceback.
    pub fn approx_eq(self, other: Self) -> bool {
        approx_eq(self.s, other.s) && approx_eq(self.h, other.h)
    }
}

/// Sums entropy and enthalpy separately, in the same order as thal.c.
impl Add for Thermo {
    type Output = Thermo;
    fn add(self, o: Thermo) -> Thermo {
        Thermo {
            s: self.s + o.s,
            h: self.h + o.h,
        }
    }
}

fn approx_eq(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() < TRACE_EPS
}
