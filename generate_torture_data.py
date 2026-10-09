import argparse
import json
import random
import math
from typing import Dict, Any

import primer3


def _fields(res) -> Dict[str, Any]:
    return {
        "structure_found": bool(getattr(res, "structure_found", False)),
        "tm": float(res.tm),
        "ds": float(res.ds),
        "dh": float(res.dh),
        "dg": float(res.dg),
    }


def is_finite_result(r: Dict[str, Any]) -> bool:
    return all(math.isfinite(r[k]) for k in ("tm", "ds", "dh", "dg"))


def random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def sample_args(rng: random.Random) -> Dict[str, Any]:
    mv_conc = rng.choice([0.0, 10.0, 25.0, 50.0, 100.0])
    dv_conc = rng.choice([0.0, 0.5, 1.0, 2.0, 4.0])
    dntp_conc = rng.choice([0.0, 0.2, 0.8, 1.2, 2.0])
    dna_conc = rng.choice([10.0, 25.0, 50.0, 200.0])
    temp_c = rng.choice([20.0, 25.0, 37.0, 50.0, 60.0])
    max_loop = rng.choice([10, 20, 30])
    output_structure = False
    return {
        "mv_conc": mv_conc,
        "dv_conc": dv_conc,
        "dntp_conc": dntp_conc,
        "dna_conc": dna_conc,
        "temp_c": temp_c,
        "max_loop": max_loop,
        "output_structure": output_structure,
    }


def main() -> int:
    p = argparse.ArgumentParser(description="Generate primer3 torture data (hairpin + heterodimer).")
    p.add_argument("--out-hairpin", default="torture_hairpin.jsonl", help="Hairpin JSONL output")
    p.add_argument("--out-heterodimer", default="torture_heterodimer.jsonl", help="Heterodimer JSONL output")
    p.add_argument("--seed", type=int, default=42, help="RNG seed")
    p.add_argument("--hairpin", type=int, default=400, help="Number of hairpin cases")
    p.add_argument("--heterodimer", type=int, default=400, help="Number of heterodimer cases")
    p.add_argument("--min-len", type=int, default=10, help="Min sequence length")
    p.add_argument("--max-len", type=int, default=60, help="Max sequence length")
    args = p.parse_args()

    rng = random.Random(args.seed)

    hairpin_written = 0
    hetero_written = 0
    hairpin_skipped = 0
    hetero_skipped = 0

    with open(args.out_hairpin, "w", encoding="utf-8") as f:
        for _ in range(args.hairpin):
            seq = random_seq(rng, rng.randint(args.min_len, args.max_len))
            a = sample_args(rng)
            res = primer3.calc_hairpin(
                seq,
                mv_conc=a["mv_conc"],
                dv_conc=a["dv_conc"],
                dntp_conc=a["dntp_conc"],
                dna_conc=a["dna_conc"],
                temp_c=a["temp_c"],
                max_loop=a["max_loop"],
                output_structure=a["output_structure"],
            )
            result = _fields(res)
            if not is_finite_result(result):
                hairpin_skipped += 1
                continue
            rec = {
                "func": "hairpin",
                "seq": seq,
                "args": a,
                "result": result,
            }
            f.write(json.dumps(rec) + "\n")
            hairpin_written += 1

    with open(args.out_heterodimer, "w", encoding="utf-8") as f:
        for _ in range(args.heterodimer):
            seq1 = random_seq(rng, rng.randint(args.min_len, args.max_len))
            seq2 = random_seq(rng, rng.randint(args.min_len, args.max_len))
            a = sample_args(rng)
            res = primer3.calc_heterodimer(
                seq1,
                seq2,
                mv_conc=a["mv_conc"],
                dv_conc=a["dv_conc"],
                dntp_conc=a["dntp_conc"],
                dna_conc=a["dna_conc"],
                temp_c=a["temp_c"],
                max_loop=a["max_loop"],
                output_structure=a["output_structure"],
            )
            result = _fields(res)
            if not is_finite_result(result):
                hetero_skipped += 1
                continue
            rec = {
                "func": "heterodimer",
                "seq1": seq1,
                "seq2": seq2,
                "args": a,
                "result": result,
            }
            f.write(json.dumps(rec) + "\n")
            hetero_written += 1

    print(
        f"wrote {args.out_hairpin} ({hairpin_written} records, {hairpin_skipped} skipped) and "
        f"{args.out_heterodimer} ({hetero_written} records, {hetero_skipped} skipped)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
