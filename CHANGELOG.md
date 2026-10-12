# Changelog

## 0.1.0

First release.

- Pure-Rust port of primer3's `thal`: `calc_hairpin`, `calc_homodimer` and
  `calc_heterodimer`, matching primer3-py 2.3.0 to within 1e-10.
- Accepts sequences as any `AsRef<[u8]>` (`&str`, `String`, `&[u8]`,
  `Vec<u8>`).
- Returns `Result<ThermoResult, ThalError>` for input that primer3 rejects.
- Optional rayon batch functions (feature `parallel`, on by default).
- Optional PyO3 bindings (feature `python`).
