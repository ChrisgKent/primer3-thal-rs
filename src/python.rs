use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::*;

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
    let args = ThalArgs {
        max_loop,
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        temp: temp_c + 273.15,
    };
    let results = py.detach(|| calc_hairpin_batch(&seqs, &args, output_structure));
    let mut out = Vec::with_capacity(results.len());
    for res in results {
        let dict = PyDict::new(py);
        dict.set_item("structure_found", res.structure_found)?;
        dict.set_item("tm", res.tm)?;
        dict.set_item("ds", res.ds)?;
        dict.set_item("dh", res.dh)?;
        dict.set_item("dg", res.dg)?;
        match res.ascii_structure {
            Some(s) => dict.set_item("ascii_structure", s)?,
            None => dict.set_item("ascii_structure", py.None())?,
        };
        out.push(dict.unbind());
    }
    Ok(out)
}

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
    let args = ThalArgs {
        max_loop,
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        temp: temp_c + 273.15,
    };
    let res = calc_heterodimer(seq1, seq2, &args, output_structure);
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
    let args = ThalArgs {
        max_loop,
        mv: mv_conc,
        dv: dv_conc,
        dntp: dntp_conc,
        dna_conc,
        temp: temp_c + 273.15,
    };
    let results = py.detach(|| calc_heterodimer_batch(&pairs, &args, output_structure));
    let mut out = Vec::with_capacity(results.len());
    for res in results {
        let dict = PyDict::new(py);
        dict.set_item("structure_found", res.structure_found)?;
        dict.set_item("tm", res.tm)?;
        dict.set_item("ds", res.ds)?;
        dict.set_item("dh", res.dh)?;
        dict.set_item("dg", res.dg)?;
        match res.ascii_structure {
            Some(s) => dict.set_item("ascii_structure", s)?,
            None => dict.set_item("ascii_structure", py.None())?,
        };
        out.push(dict.unbind());
    }
    Ok(out)
}

#[pymodule]
fn primer3_py_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calc_hairpin_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_hairpin_batch_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_heterodimer_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_heterodimer_batch_py, m)?)?;
    m.add("__doc__", "Rust hairpin thermodynamics (PyO3 bindings)")?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
