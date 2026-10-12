#![allow(clippy::too_many_arguments)] // mirrors primer3-py's keyword arguments

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::*;

// Match primer3-py: bad max_loop is a ValueError, sequence problems a RuntimeError.
impl From<ThalError> for PyErr {
    fn from(e: ThalError) -> PyErr {
        match e {
            ThalError::InvalidMaxLoop(_) => PyValueError::new_err(e.to_string()),
            _ => PyRuntimeError::new_err(e.to_string()),
        }
    }
}

/// `ThalArgs` from primer3-py-style keyword arguments (temperature in °C).
fn thal_args(
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
) -> ThalArgs {
    ThalArgs {
        max_loop,
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        temp: 0.0,
    }
    .with_temp_c(temp_c)
}

/// A primer3-py `ThermoResult`-like dict.
fn to_dict(py: Python<'_>, res: ThermoResult) -> PyResult<Py<PyDict>> {
    let out = PyDict::new(py);
    out.set_item("structure_found", res.structure_found)?;
    out.set_item("tm", res.tm)?;
    out.set_item("ds", res.ds)?;
    out.set_item("dh", res.dh)?;
    out.set_item("dg", res.dg)?;
    out.set_item("ascii_structure", res.ascii_structure)?;
    Ok(out.unbind())
}

/// Dicts for a batch, failing on the first error.
fn to_dicts(
    py: Python<'_>,
    results: Vec<Result<ThermoResult, ThalError>>,
) -> PyResult<Vec<Py<PyDict>>> {
    results.into_iter().map(|r| to_dict(py, r?)).collect()
}

/// Hairpin thermodynamics of `seq`, as primer3.calc_hairpin.
///
/// Concentrations are mM (mv_conc, dv_conc, dntp_conc) and nM (dna_conc);
/// temp_c is in Celsius. Returns a dict with structure_found, tm (°C),
/// dg/dh (cal/mol), ds (cal/K/mol) and ascii_structure (str or None).
/// Raises RuntimeError for an empty or >60 nt sequence, ValueError for
/// max_loop outside 0..=30.
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
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    to_dict(py, calc_hairpin(seq, &args, output_structure)?)
}

/// calc_hairpin_py for a list of sequences, computed in parallel without
/// the GIL. Returns the dicts in input order; raises on the first error.
#[pyfunction]
fn calc_hairpin_batch_py(
    py: Python<'_>,
    seqs: Vec<String>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
) -> PyResult<Vec<Py<PyDict>>> {
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    let results = py.detach(|| calc_hairpin_batch(&seqs, &args, output_structure));
    to_dicts(py, results)
}

/// Homodimer thermodynamics of `seq`, as primer3.calc_homodimer.
/// Arguments and result as calc_hairpin_py.
#[pyfunction]
fn calc_homodimer_py(
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
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    to_dict(py, calc_homodimer(seq, &args, output_structure)?)
}

/// calc_homodimer_py for a list of sequences, computed in parallel without
/// the GIL. Returns the dicts in input order; raises on the first error.
#[pyfunction]
fn calc_homodimer_batch_py(
    py: Python<'_>,
    seqs: Vec<String>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
) -> PyResult<Vec<Py<PyDict>>> {
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    let results = py.detach(|| calc_homodimer_batch(&seqs, &args, output_structure));
    to_dicts(py, results)
}

/// Heterodimer thermodynamics of seq1 binding seq2 (both 5'->3'), as
/// primer3.calc_heterodimer. Arguments and result as calc_hairpin_py; only
/// raises for length if both sequences exceed 60 nt.
#[pyfunction]
fn calc_heterodimer_py(
    py: Python<'_>,
    seq1: &str,
    seq2: &str,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
) -> PyResult<Py<PyDict>> {
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    to_dict(py, calc_heterodimer(seq1, seq2, &args, output_structure)?)
}

/// calc_heterodimer_py for a list of (seq1, seq2) pairs, computed in
/// parallel without the GIL. Returns the dicts in input order; raises on
/// the first error.
#[pyfunction]
fn calc_heterodimer_batch_py(
    py: Python<'_>,
    pairs: Vec<(String, String)>,
    mv_conc: f64,
    dv_conc: f64,
    dntp_conc: f64,
    dna_conc: f64,
    temp_c: f64,
    max_loop: i32,
    output_structure: bool,
) -> PyResult<Vec<Py<PyDict>>> {
    let args = thal_args(mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop);
    let results = py.detach(|| calc_heterodimer_batch(&pairs, &args, output_structure));
    to_dicts(py, results)
}

/// The `primer3_thal_rs` Python module.
#[pymodule]
fn primer3_thal_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_hairpin_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_hairpin_batch_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_homodimer_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_homodimer_batch_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_heterodimer_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_heterodimer_batch_py, m)?)?;
    m.add("__doc__", "Rust port of primer3 thal (PyO3 bindings)")?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
