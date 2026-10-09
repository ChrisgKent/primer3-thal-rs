import argparse
import json
from typing import Dict, Any, Iterable

from primer3_py_rs import calc_hairpin_py, calc_heterodimer_py


def _fields(res) -> Dict[str, Any]:
    if isinstance(res, dict):
        return {
            "structure_found": bool(res.get("structure_found", False)),
            "tm": float(res["tm"]),
            "ds": float(res["ds"]),
            "dh": float(res["dh"]),
            "dg": float(res["dg"]),
        }
    return {
        "structure_found": bool(getattr(res, "structure_found", False)),
        "tm": float(res.tm),
        "ds": float(res.ds),
        "dh": float(res.dh),
        "dg": float(res.dg),
    }


def iter_records(path: str) -> Iterable[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            yield json.loads(line)


def main() -> int:
    p = argparse.ArgumentParser(description="Run rust bindings against saved primer3 torture data.")
    p.add_argument("--hairpin", default="torture_hairpin.jsonl", help="Hairpin JSONL file")
    p.add_argument("--heterodimer", default="torture_heterodimer.jsonl", help="Heterodimer JSONL file")
    p.add_argument("--tm", type=float, default=0.05, help="TM tolerance")
    p.add_argument("--ds", type=float, default=0.2, help="dS tolerance")
    p.add_argument("--dh", type=float, default=1.0, help="dH tolerance")
    p.add_argument("--dg", type=float, default=0.5, help="dG tolerance")
    args = p.parse_args()

    count = 0
    for path in (args.hairpin, args.heterodimer):
        for rec in iter_records(path):
            func = rec["func"]
            a = rec["args"]
            expected = rec["result"]

            if func == "hairpin":
                seq = rec["seq"]
                res = calc_hairpin_py(
                    seq,
                    a["mv_conc"],
                    a["dv_conc"],
                    a["dntp_conc"],
                    a["dna_conc"],
                    a["temp_c"],
                    a["max_loop"],
                    a["output_structure"],
                )
            elif func == "heterodimer":
                seq1 = rec["seq1"]
                seq2 = rec["seq2"]
                res = calc_heterodimer_py(
                    seq1,
                    seq2,
                    a["mv_conc"],
                    a["dv_conc"],
                    a["dntp_conc"],
                    a["dna_conc"],
                    a["temp_c"],
                    a["max_loop"],
                    a["output_structure"],
                )
            else:
                raise ValueError(f"Unknown func: {func}")

            got = _fields(res)
            if got["structure_found"] != expected["structure_found"]:
                raise AssertionError((func, "structure_found", expected, got))
            if abs(got["tm"] - expected["tm"]) > args.tm:
                raise AssertionError((func, "tm", expected, got))
            if abs(got["ds"] - expected["ds"]) > args.ds:
                raise AssertionError((func, "ds", expected, got))
            if abs(got["dh"] - expected["dh"]) > args.dh:
                raise AssertionError((func, "dh", expected, got))
            if abs(got["dg"] - expected["dg"]) > args.dg:
                raise AssertionError((func, "dg", expected, got))

            count += 1

    print(f"ok: {count} cases matched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
