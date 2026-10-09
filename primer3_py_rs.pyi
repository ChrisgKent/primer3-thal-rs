from typing import Optional, TypedDict, List, Tuple

class ThermoResult(TypedDict):
    structure_found: bool
    tm: float
    ds: float
    dh: float
    dg: float
    ascii_structure: Optional[str]

version: str


def calc_hairpin_py(
    seq: str,
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
    temp_c: float,
    max_loop: int,
    output_structure: bool,
) -> ThermoResult:
    """
    Compute hairpin thermodynamics for a single DNA sequence.

    Parameters
    - seq: DNA sequence (A/C/G/T).
    - mv_conc: monovalent cation concentration (mM).
    - dv_conc: divalent cation concentration (mM).
    - dntp_conc: dNTP concentration (mM).
    - dna_conc: DNA concentration (nM).
    - temp_c: temperature in Celsius.
    - max_loop: maximum loop size.
    - output_structure: whether to compute ASCII structure.

    Returns
    - ThermoResult dict with tm/ds/dh/dg and structure_found.
    - ascii_structure is provided when output_structure=True.
    """
    ...


def calc_hairpin_batch_py(
    seqs: List[str],
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
    temp_c: float,
    max_loop: int,
    output_structure: bool,
) -> List[ThermoResult]:
    """
    Compute hairpin thermodynamics for many sequences in parallel.

    Parameters match calc_hairpin_py, with seqs as a list of sequences.

    Returns
    - List of ThermoResult dicts in the same order as input.
    """
    ...


def calc_heterodimer_py(
    seq1: str,
    seq2: str,
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
    temp_c: float,
    max_loop: int,
    output_structure: bool,
) -> ThermoResult:
    """
    Compute heterodimer thermodynamics for two DNA sequences (thal_any alignment).

    Parameters match calc_hairpin_py; seq1/seq2 are the two sequences.

    Returns
    - ThermoResult dict with tm/ds/dh/dg and structure_found.
    """
    ...


def calc_heterodimer_batch_py(
    pairs: List[Tuple[str, str]],
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
    temp_c: float,
    max_loop: int,
    output_structure: bool,
) -> List[ThermoResult]:
    """
    Compute heterodimer thermodynamics for many pairs in parallel.

    Parameters match calc_heterodimer_py, with pairs as (seq1, seq2).

    Returns
    - List of ThermoResult dicts in the same order as input.
    """
    ...
