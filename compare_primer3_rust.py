import random

import primer3
import primer3_py_rs


def _fields(res):
    return res["tm"], res["ds"], res["dh"], res["dg"]


def main() -> None:
    random.seed(42)
    seqs = [
        "CCCCCATCCGATCAGGGGG",
        "GTGTCCAATCGACATTTAGGAGCTACGCGGCCAGGGGCAAA",
        "TGCATTGAATGCGTGTCACGTTATGCACGC",
        "GATTACAGATTACAGATTACA",
        "ACGTACGTACGTACGTACGT",
    ]
    for _ in range(10):
        seqs.append("".join(random.choice("ACGT") for _ in range(random.randint(20, 59))))

    for seq in seqs:
        py = primer3.calc_hairpin(seq)
        ru = primer3_py_rs.calc_hairpin(seq)
        py_tm, py_ds, py_dh, py_dg = py.tm, py.ds, py.dh, py.dg
        ru_tm, ru_ds, ru_dh, ru_dg = _fields(ru)
        assert abs(py_tm - ru_tm) < 0.05, (seq, py_tm, ru_tm)
        assert abs(py_ds - ru_ds) < 0.2, (seq, py_ds, ru_ds)
        assert abs(py_dh - ru_dh) < 1.0, (seq, py_dh, ru_dh)
        assert abs(py_dg - ru_dg) < 0.5, (seq, py_dg, ru_dg)


if __name__ == "__main__":
    main()
