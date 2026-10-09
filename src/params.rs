//! Nearest-neighbour parameter tables, embedded from primer3's config files.

use std::str::Lines;
use std::sync::OnceLock;

use crate::thermo::Thermo;
use crate::{AT_H, AT_S, TINY_ENTROPY};

/// A nucleotide as thal encodes it. Anything that isn't A/C/G/T is `N`,
/// which pairs with nothing.
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl Base {
    /// thal.c's `str2int`: case-insensitive, and also accepts the digits 0-3.
    pub fn from_ascii(b: u8) -> Base {
        match b {
            b'A' | b'a' | b'0' => Base::A,
            b'C' | b'c' | b'1' => Base::C,
            b'G' | b'g' | b'2' => Base::G,
            b'T' | b't' | b'3' => Base::T,
            _ => Base::N,
        }
    }

    /// Watson-Crick pairing only.
    pub fn pairs_with(self, other: Base) -> bool {
        matches!(
            (self, other),
            (Base::A, Base::T) | (Base::T, Base::A) | (Base::C, Base::G) | (Base::G, Base::C)
        )
    }
}

type Table2 = [[Thermo; 5]; 5];
type Table3 = [[[Thermo; 5]; 5]; 5];
type Table4 = [[[[Thermo; 5]; 5]; 5]; 5];

pub(crate) struct ThermoParams {
    stack: Table4,
    stack_mm: Table4,
    tstack: Table4,
    tstack2: Table4,
    dangle3: Table3,
    dangle5: Table3,
    atp: Table2,
    pub interior_loop: [Thermo; 30],
    pub bulge_loop: [Thermo; 30],
    pub hairpin_loop: [Thermo; 30],
    triloops: Vec<([Base; 5], Thermo)>,
    tetraloops: Vec<([Base; 6], Thermo)>,
}

#[rustfmt::skip]
impl ThermoParams {
    /// Watson-Crick stack of pair a-c on pair b-d.
    pub fn stack(&self, a: Base, b: Base, c: Base, d: Base) -> Thermo {
        self.stack[a as usize][b as usize][c as usize][d as usize]
    }
    /// Single internal mismatch (1x1 loop) stack.
    pub fn stack_mm(&self, a: Base, b: Base, c: Base, d: Base) -> Thermo {
        self.stack_mm[a as usize][b as usize][c as usize][d as usize]
    }
    /// Terminal mismatch inside an internal loop.
    pub fn tstack(&self, a: Base, b: Base, c: Base, d: Base) -> Thermo {
        self.tstack[a as usize][b as usize][c as usize][d as usize]
    }
    /// Terminal mismatch at a helix end or hairpin closing pair.
    pub fn tstack2(&self, a: Base, b: Base, c: Base, d: Base) -> Thermo {
        self.tstack2[a as usize][b as usize][c as usize][d as usize]
    }
    pub fn dangle3(&self, a: Base, b: Base, c: Base) -> Thermo {
        self.dangle3[a as usize][b as usize][c as usize]
    }
    pub fn dangle5(&self, a: Base, b: Base, c: Base) -> Thermo {
        self.dangle5[a as usize][b as usize][c as usize]
    }
    /// Terminal AT penalty for the pair a-b.
    pub fn atp(&self, a: Base, b: Base) -> Thermo {
        self.atp[a as usize][b as usize]
    }
    /// Bonus for a hairpin whose closing pair + 3 nt loop is `key`.
    pub fn triloop(&self, key: &[Base]) -> Option<Thermo> {
        find_loop(&self.triloops, key)
    }
    /// Bonus for a hairpin whose closing pair + 4 nt loop is `key`.
    pub fn tetraloop(&self, key: &[Base]) -> Option<Thermo> {
        find_loop(&self.tetraloops, key)
    }
}

fn find_loop<const N: usize>(table: &[([Base; N], Thermo)], key: &[Base]) -> Option<Thermo> {
    table
        .binary_search_by(|(k, _)| k.as_slice().cmp(key))
        .ok()
        .map(|i| table[i].1)
}

pub(crate) fn params() -> &'static ThermoParams {
    static PARAMS: OnceLock<ThermoParams> = OnceLock::new();
    PARAMS.get_or_init(|| ThermoParams::load().expect("embedded primer3 parameter tables"))
}

macro_rules! config {
    ($name:literal) => {
        include_str!(concat!("../primer3_config/", $name))
    };
}

/// One `.ds` file and its `.dh` partner, read value by value in step.
struct Pair<'a> {
    name: &'static str,
    ds: Lines<'a>,
    dh: Lines<'a>,
}

impl<'a> Pair<'a> {
    fn new(name: &'static str, ds: &'a str, dh: &'a str) -> Self {
        Self {
            name,
            ds: ds.lines(),
            dh: dh.lines(),
        }
    }

    fn next_line(lines: &mut Lines<'a>, name: &str) -> Result<&'a str, String> {
        lines
            .find(|l| !l.trim().is_empty())
            .ok_or_else(|| format!("{name}: unexpected end of table"))
    }

    /// The next entry; non-finite values become `Thermo::INVALID`, as in thal.c.
    fn next(&mut self) -> Result<Thermo, String> {
        let s = parse_value(Self::next_line(&mut self.ds, self.name)?, self.name)?;
        let h = parse_value(Self::next_line(&mut self.dh, self.name)?, self.name)?;
        Ok(if s.is_finite() && h.is_finite() {
            Thermo::new(s, h)
        } else {
            Thermo::INVALID
        })
    }

    /// The next whitespace-separated rows of both files.
    fn next_rows(&mut self) -> Result<Option<(&'a str, &'a str)>, String> {
        let ds = self.ds.find(|l| !l.trim().is_empty());
        let dh = self.dh.find(|l| !l.trim().is_empty());
        match (ds, dh) {
            (Some(s), Some(h)) => Ok(Some((s, h))),
            (None, None) => Ok(None),
            _ => Err(format!("{}: .ds and .dh differ in length", self.name)),
        }
    }
}

fn parse_value(token: &str, name: &str) -> Result<f64, String> {
    let t = token.trim();
    if t.eq_ignore_ascii_case("inf") {
        Ok(f64::INFINITY)
    } else {
        t.parse()
            .map_err(|_| format!("{name}: bad value {t:?}"))
    }
}

/// A 5x5x5x5 table where `edge` supplies fixed values for some entries and
/// the rest are read from the files in index order.
#[allow(clippy::needless_range_loop)] // indices also feed `edge`
fn read_table4(
    mut pair: Pair<'_>,
    edge: impl Fn([usize; 4]) -> Option<Thermo>,
) -> Result<Table4, String> {
    let mut t = [[[[Thermo::ZERO; 5]; 5]; 5]; 5];
    for a in 0..5 {
        for b in 0..5 {
            for c in 0..5 {
                for d in 0..5 {
                    t[a][b][c][d] = match edge([a, b, c, d]) {
                        Some(v) => v,
                        None => pair.next()?,
                    };
                }
            }
        }
    }
    Ok(t)
}

/// Stack tables: anything involving N is forbidden.
fn read_stack(pair: Pair<'_>) -> Result<Table4, String> {
    read_table4(pair, |idx| idx.contains(&4).then_some(Thermo::INVALID))
}

/// Terminal-mismatch tables: an N in a pair is forbidden; an N as the
/// mismatching base contributes (nearly) nothing.
fn read_tstack(pair: Pair<'_>) -> Result<Table4, String> {
    read_table4(pair, |[a, b, c, d]| {
        if a == 4 || c == 4 {
            Some(Thermo::INVALID)
        } else if b == 4 || d == 4 {
            Some(Thermo::new(TINY_ENTROPY, 0.0))
        } else {
            None
        }
    })
}

/// dangle.ds/dh holds the 3' table then the 5' table. The 3' file order is
/// (a, c, b) for entry `[a][b][c]`, as in thal.c.
#[allow(clippy::needless_range_loop)] // the 3' table is filled in permuted order
fn read_dangles(mut pair: Pair<'_>) -> Result<(Table3, Table3), String> {
    let mut d3 = [[[Thermo::ZERO; 5]; 5]; 5];
    let mut d5 = [[[Thermo::ZERO; 5]; 5]; 5];
    for a in 0..5 {
        for b in 0..5 {
            for c in 0..5 {
                d3[a][c][b] = if a == 4 || b == 4 || c == 4 {
                    Thermo::INVALID
                } else {
                    pair.next()?
                };
            }
        }
    }
    for a in 0..5 {
        for b in 0..5 {
            for c in 0..5 {
                d5[a][b][c] = if a == 4 || b == 4 || c == 4 {
                    Thermo::INVALID
                } else {
                    pair.next()?
                };
            }
        }
    }
    Ok((d3, d5))
}

/// loops.ds/dh rows are `size interior bulge hairpin`.
fn read_loops(mut pair: Pair<'_>) -> Result<[[Thermo; 30]; 3], String> {
    let mut out = [[Thermo::ZERO; 30]; 3];
    for k in 0..30 {
        let (ds, dh) = pair
            .next_rows()?
            .ok_or_else(|| format!("{}: expected 30 rows", pair.name))?;
        let ds: Vec<&str> = ds.split_whitespace().collect();
        let dh: Vec<&str> = dh.split_whitespace().collect();
        for (col, table) in out.iter_mut().enumerate() {
            let field = |row: &[&str]| -> Result<f64, String> {
                let tok = row
                    .get(col + 1)
                    .ok_or_else(|| format!("{}: short row", pair.name))?;
                parse_value(tok, pair.name)
            };
            table[k] = Thermo::new(field(&ds)?, field(&dh)?);
        }
    }
    Ok(out)
}

/// Tri/tetraloop bonus tables, sorted by sequence for binary search. The
/// .ds and .dh files must list the same loops in the same order.
fn read_loop_bonus<const N: usize>(mut pair: Pair<'_>) -> Result<Vec<([Base; N], Thermo)>, String> {
    let mut out = Vec::new();
    while let Some((ds, dh)) = pair.next_rows()? {
        let (key_s, s) = ds
            .split_once(char::is_whitespace)
            .ok_or_else(|| format!("{}: bad row {ds:?}", pair.name))?;
        let (key_h, h) = dh
            .split_once(char::is_whitespace)
            .ok_or_else(|| format!("{}: bad row {dh:?}", pair.name))?;
        if key_s != key_h || key_s.len() != N {
            return Err(format!("{}: mismatched loop {key_s} / {key_h}", pair.name));
        }
        let mut key = [Base::N; N];
        for (k, b) in key.iter_mut().zip(key_s.bytes()) {
            *k = Base::from_ascii(b);
        }
        out.push((key, Thermo::new(parse_value(s, pair.name)?, parse_value(h, pair.name)?)));
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(out)
}

impl ThermoParams {
    fn load() -> Result<Self, String> {
        let (dangle3, dangle5) = read_dangles(Pair::new(
            "dangle",
            config!("dangle.ds"),
            config!("dangle.dh"),
        ))?;
        let [interior_loop, bulge_loop, hairpin_loop] =
            read_loops(Pair::new("loops", config!("loops.ds"), config!("loops.dh")))?;

        let mut atp = [[Thermo::new(TINY_ENTROPY, 0.0); 5]; 5];
        let at = Thermo::new(AT_S, AT_H);
        atp[Base::A as usize][Base::T as usize] = at;
        atp[Base::T as usize][Base::A as usize] = at;

        Ok(Self {
            stack: read_stack(Pair::new("stack", config!("stack.ds"), config!("stack.dh")))?,
            stack_mm: read_stack(Pair::new(
                "stackmm",
                config!("stackmm.ds"),
                config!("stackmm.dh"),
            ))?,
            tstack: read_tstack(Pair::new(
                "tstack",
                config!("tstack_tm_inf.ds"),
                config!("tstack.dh"),
            ))?,
            tstack2: read_tstack(Pair::new(
                "tstack2",
                config!("tstack2.ds"),
                config!("tstack2.dh"),
            ))?,
            dangle3,
            dangle5,
            atp,
            interior_loop,
            bulge_loop,
            hairpin_loop,
            triloops: read_loop_bonus(Pair::new(
                "triloop",
                config!("triloop.ds"),
                config!("triloop.dh"),
            ))?,
            tetraloops: read_loop_bonus(Pair::new(
                "tetraloop",
                config!("tetraloop.ds"),
                config!("tetraloop.dh"),
            ))?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_tables_parse() {
        let p = ThermoParams::load().unwrap();
        assert_eq!(p.triloops.len(), 16);
        assert_eq!(p.tetraloops.len(), 77);
        // Spot-check a value from each file against primer3_config.
        assert_eq!(p.stack(Base::A, Base::A, Base::T, Base::T), Thermo::new(-22.2, -7900.0));
        assert_eq!(p.hairpin_loop[2], Thermo::new(-11.28, 0.0));
        assert!(p.stack(Base::N, Base::A, Base::T, Base::T).h.is_infinite());
    }
}
