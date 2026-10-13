# Changelog

## 0.2.0

Additions beyond primer3. Nothing changes with the default `ThalArgs`, and
results still match primer3-py 2.3.0.

- `ThalArgs` gains `dmso_conc`, `dmso_fact` and `formamide_conc` (defaults 0,
  0.6, 0). When set, they correct the reported Tm as primer3's oligotm does.
- New `nn` module: thal's nearest-neighbour model for a perfectly matched
  duplex (constants, salt and DMSO/formamide corrections, `duplex_tm`,
  `duplex_annealing`, and the incremental `Duplex`). Its Tm equals
  `calc_heterodimer(seq, reverse_complement(seq)).tm`.
- `nn::bound_percent`: percentage bound at a temperature from a ΔH/ΔS, such
  as a `ThermoResult`'s, on the same terms as `duplex_annealing`.
- New `ThalError::InvalidBase` for the `nn` functions, which accept only
  A/C/G/T.
- Python bindings take `dmso_conc`, `dmso_fact` and `formamide_conc` keyword
  arguments, and `primer3_thal_rs.pyi` lists them.
- Degenerate duplexes (a single base, no salt at all) give thal's numbers,
  as primer3 does, rather than an error; tested against thal and primer3-py.
- Faster `nn`: a constant Watson-Crick stack table, and `bound_percent` in
  log space (one `exp`, no overflow). Hairpin and dimer calls skip the GC
  count unless DMSO or formamide is set.

Breaking: `ThalArgs` and `ThalError` have new fields and variants, so code
that builds `ThalArgs` with every field listed, or matches `ThalError`
exhaustively, needs updating.

## 0.1.0

First release.

- LLM transpilation of primer3-py v2.3.0's `thal` (bundled `thal.c`, parameter
  tables and `thermoanalysis.pyx` defaults) to Rust: `calc_hairpin`, `calc_homodimer` and
  `calc_heterodimer`, matching primer3-py 2.3.0 to within 1e-10.
- Accepts sequences as any `AsRef<[u8]>` (`&str`, `String`, `&[u8]`,
  `Vec<u8>`).
- Returns `Result<ThermoResult, ThalError>` for input that primer3 rejects.
- Optional rayon batch functions (feature `parallel`, on by default).
- Optional PyO3 bindings (feature `python`).
