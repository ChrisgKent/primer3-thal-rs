import os
import random

import primer3
from primer3_py_rs import calc_hairpin_py, calc_heterodimer_py


def _fields(res):
    if isinstance(res, dict):
        return res["tm"], res["ds"], res["dh"], res["dg"]
    return res.tm, res.ds, res.dh, res.dg


def main() -> int:
    seed = int(os.environ.get("SEED", "42"))
    count = int(os.environ.get("COUNT", "10"))
    min_len = int(os.environ.get("MIN_LEN", "20"))
    max_len = int(os.environ.get("MAX_LEN", "59"))

    random.seed(seed)
    seqs = [
        "CCCCCATCCGATCAGGGGG",
        "GTGTCCAATCGACATTTAGGAGCTACGCGGCCAGGGGCAAA",
        "TGCATTGAATGCGTGTCACGTTATGCACGC",
        "GATTACAGATTACAGATTACA",
        "ACGTACGTACGTACGTACGT",
    ]
    for _ in range(count):
        seqs.append(
            "".join(
                random.choice("ACGT") for _ in range(random.randint(min_len, max_len))
            )
        )

    mv_conc = 50.0
    dv_conc = 0.0
    dntp_conc = 0.8
    dna_conc = 50.0
    temp_c = 37.0
    max_loop = 30
    output_structure = False

    for seq in seqs:
        py = primer3.calc_hairpin(
            seq,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
            temp_c=temp_c,
            max_loop=max_loop,
            output_structure=output_structure,
        )
        ru_py = calc_hairpin_py(
            seq,
            mv_conc,
            dv_conc,
            dntp_conc,
            dna_conc,
            temp_c,
            max_loop,
            output_structure,
        )
        print(seq, py, ru_py, sep="\n")
        py_tm, py_ds, py_dh, py_dg = _fields(py)
        rpy_tm, rpy_ds, rpy_dh, rpy_dg = _fields(ru_py)

        assert abs(py_tm - rpy_tm) < 0.05, (seq, py_tm, rpy_tm)
        assert abs(py_ds - rpy_ds) < 0.2, (seq, py_ds, rpy_ds)
        assert abs(py_dh - rpy_dh) < 1.0, (seq, py_dh, rpy_dh)
        assert abs(py_dg - rpy_dg) < 0.5, (seq, py_dg, rpy_dg)

    for seq1 in seqs:
        for seq2 in seqs:
            py = primer3.calc_heterodimer(
                seq1,
                seq2,
                mv_conc=mv_conc,
                dv_conc=dv_conc,
                dntp_conc=dntp_conc,
                dna_conc=dna_conc,
                temp_c=temp_c,
                max_loop=max_loop,
                output_structure=output_structure,
            )
            ru = calc_heterodimer_py(
                seq1,
                seq2,
                mv_conc,
                dv_conc,
                dntp_conc,
                dna_conc,
                temp_c,
                max_loop,
                output_structure,
            )
            print(seq1, seq2, py, ru, sep="\n")
            py_tm, py_ds, py_dh, py_dg = _fields(py)
            ru_tm, ru_ds, ru_dh, ru_dg = _fields(ru)

            assert abs(py_tm - ru_tm) < 0.05, (seq1, seq2, py_tm, ru_tm)
            assert abs(py_ds - ru_ds) < 0.2, (seq1, seq2, py_ds, ru_ds)
            assert abs(py_dh - ru_dh) < 1.0, (seq1, seq2, py_dh, ru_dh)
            assert abs(py_dg - ru_dg) < 0.5, (seq1, seq2, py_dg, ru_dg)

    print(f"ok: {len(seqs)} sequences matched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
