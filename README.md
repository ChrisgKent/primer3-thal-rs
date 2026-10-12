# primer3-thal-rs

A pure-Rust port of [primer3](https://github.com/primer3-org/primer3)'s
`thal` thermodynamic alignment. It calculates hairpin, homodimer and
heterodimer ΔG, ΔH, ΔS and Tm for DNA oligos, with results that match
[primer3-py](https://github.com/libnano/primer3-py).

- **Exact.** About 18,000 generated cases are checked against primer3-py
  2.3.0. All agree to within 1e-10 (the tests allow 1e-6), and ASCII structures
  match character for character.
- **Self-contained.** No C code and no runtime files: primer3's parameter
  tables are compiled in. With default features off, the crate has no
  dependencies.
- **Fast.** On an Apple M-series machine, a mean call on an 18–30 nt primer
  takes about 18 µs for a hairpin and about 94 µs for a heterodimer of two
  primers.

## Usage

```toml
[dependencies]
primer3-thal-rs = "0.1"
```

```rust
use primer3_thal_rs::{ThalArgs, ThalError, calc_hairpin, calc_heterodimer, calc_homodimer};

fn main() -> Result<(), ThalError> {
    // primer3-py's defaults: 50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM DNA, 37 °C.
    let args = ThalArgs::default();

    let hairpin = calc_hairpin("CCCCCATCCGATCAGGGGG", &args, false)?;
    println!("hairpin Tm {:.2} °C, ΔG {:.0} cal/mol", hairpin.tm, hairpin.dg);

    // Sequences are 5'->3' and can be &str, String, &[u8] or Vec<u8>.
    let primer: &[u8] = b"GTAAAACGACGGCCAGT";
    let dimer = calc_heterodimer(primer, "ACTGGCCGTCGTTTTAC", &args, true)?;
    println!("{}", dimer.ascii_structure.unwrap());

    // Other conditions, with the temperature in °C.
    let pcr = ThalArgs { mv: 100.0, dv: 2.0, dntp: 0.8, ..ThalArgs::default() }.with_temp_c(60.0);
    let homo = calc_homodimer(primer, &pcr, false)?;
    println!("homodimer ΔG at 60 °C: {:.0} cal/mol", homo.dg);
    Ok(())
}
```

| Function | primer3-py equivalent |
|---|---|
| `calc_hairpin(seq, &args, output_structure)` | `primer3.calc_hairpin` |
| `calc_homodimer(seq, &args, output_structure)` | `primer3.calc_homodimer` |
| `calc_heterodimer(seq1, seq2, &args, output_structure)` | `primer3.calc_heterodimer` |
| `calc_*_batch(...)` (feature `parallel`) | many inputs in parallel |

Each call returns a `ThermoResult`:

| Field | Meaning |
|---|---|
| `structure_found` | whether any structure was found |
| `tm` | melting temperature, °C |
| `dg` | ΔG at `args.temp`, cal/mol |
| `dh` | ΔH, cal/mol |
| `ds` | ΔS, cal/K/mol |
| `ascii_structure` | primer3's drawing, if requested |

Calls return `ThalError` for input that primer3 also rejects:
- an empty sequence
- a hairpin or homodimer longer than 60 nt
- a heterodimer where *both* sequences are longer than 60 nt
- `max_loop` outside 0–30

### Input handling

The same rules as primer3 apply:
- Sequences are case-insensitive.
- Any byte other than A/C/G/T, such as N or an ambiguity code, is a base that
  pairs with nothing.
- Pass primers exactly as written, 5′→3′. For a heterodimer the second
  sequence is reversed internally.

## Features

| Feature | Default | Purpose |
|---|---|---|
| `parallel` | yes | `calc_*_batch` functions, using rayon |
| `python` | no | PyO3 bindings, built with [maturin](https://www.maturin.rs/) as the `primer3_thal_rs` module |

To build the Python module (Python ≥ 3.8):

```sh
maturin build --release       # or: maturin develop
```

## Development

```sh
cargo test --release          # fixtures, error handling, no-panic property tests
cargo bench --bench batch     # mean time per call on 18-30 nt primers
```

The reference fixtures in `tests/data/` are generated from primer3-py by
`tests/gen_fixtures.py`. That script is seeded, so it produces the same
fixtures on every run. To regenerate them:

```sh
pip install "primer3-py==2.3.0"
python tests/gen_fixtures.py
```

The algorithm is a faithful port of `thal.c`. Every comparison and every
floating-point expression keeps primer3's order of operations, which is what
makes the results match exactly. Any change that alters the "worst |diff|"
reported by the fixture tests is a regression.

## Licence and credits

This crate is derived from primer3's `thal.c` and its parameter tables:

> Copyright (c) 1996–2016 Whitehead Institute for Biomedical Research,
> Steve Rozen and Helen Skaletsky.

They are licensed under the GNU General Public License, version 2 or later.
This crate is distributed under the same terms, **GPL-2.0-or-later**; see
[LICENSE](LICENSE).

The reference behaviour and default conditions follow
[primer3-py](https://github.com/libnano/primer3-py) by Ben Pruitt, Nick
Conway and contributors.
