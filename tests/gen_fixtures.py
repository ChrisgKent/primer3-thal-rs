"""Generate reference fixtures from primer3-py for the Rust thal port.

Usage: python tests/gen_fixtures.py [--seed 42] [--scale 1.0]

Writes tests/data/{hairpin,homodimer,heterodimer}.jsonl. The first line of each
file is a header recording the primer3-py version; every other line is a case:
  {"seq": ..., ["seq2": ...], "args": {...}, "result": {...}}  or
  {"seq": ..., "args": {...}, "error": "RuntimeError"}
"""

import argparse
import json
import math
import random
from pathlib import Path

import primer3

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "primer3_config"
OUT = ROOT / "tests" / "data"

COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")

MV = [0.0, 10.0, 25.0, 50.0, 100.0, 200.0]
DV = [0.0, 0.5, 1.5, 2.0, 4.0, 10.0]
DNTP = [0.0, 0.2, 0.6, 0.8, 2.0]
DNA = [10.0, 15.0, 50.0, 250.0, 1000.0]
TEMP = [20.0, 25.0, 37.0, 50.0, 60.0, 65.0, 72.0]


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n, gc=0.5):
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=w, k=n))


def rand_args(rng, output_structure=False):
    return {
        "mv_conc": rng.choice(MV),
        "dv_conc": rng.choice(DV),
        "dntp_conc": rng.choice(DNTP),
        "dna_conc": rng.choice(DNA),
        "temp_c": rng.choice(TEMP),
        "max_loop": rng.randint(0, 30),
        "output_structure": output_structure,
    }


def default_args(output_structure=False):
    return {
        "mv_conc": 50.0,
        "dv_conc": 1.5,
        "dntp_conc": 0.6,
        "dna_conc": 50.0,
        "temp_c": 37.0,
        "max_loop": 30,
        "output_structure": output_structure,
    }


def loop_keys(name):
    return [line.split()[0] for line in (CONFIG / name).read_text().splitlines() if line.strip()]


def mutate(rng, s, n):
    s = list(s)
    for _ in range(n):
        s[rng.randrange(len(s))] = rng.choice("ACGT")
    return "".join(s)


def stem_loop(rng, loop):
    """stem + loop + revcomp(stem), optionally with a mismatch or bulge."""
    stem = rand_seq(rng, rng.randint(3, 12), gc=rng.uniform(0.3, 0.8))
    right = revcomp(stem)
    r = rng.random()
    if r < 0.2:
        right = mutate(rng, right, 1)
    elif r < 0.3:
        k = rng.randrange(len(right) + 1)
        right = right[:k] + rng.choice("ACGT") + right[k:]
    flank5 = rand_seq(rng, rng.randint(0, 8))
    flank3 = rand_seq(rng, rng.randint(0, 8))
    return (flank5 + stem + loop + right + flank3)[:60]


def num(x):
    """JSON has no inf/nan; zero-salt conditions produce them, so keep as strings."""
    return x if math.isfinite(x) else str(x)


def run(func, seqs, args):
    kw = dict(args)
    try:
        res = func(*seqs, **kw)
    except Exception as e:  # primer3-py raises RuntimeError / ValueError
        return {"error": type(e).__name__}
    out = {
        "structure_found": bool(res.structure_found),
        "tm": num(res.tm),
        "dg": num(res.dg),
        "dh": num(res.dh),
        "ds": num(res.ds),
    }
    if args["output_structure"]:
        out["ascii_structure"] = res.ascii_structure
    return {"result": out}


def hairpin_cases(rng, scale):
    n = lambda k: max(1, int(k * scale))
    cases = []
    # Random sequences across the full allowed length range.
    for _ in range(n(3000)):
        cases.append((rand_seq(rng, rng.randint(1, 60), rng.uniform(0.2, 0.8)), rand_args(rng)))
    for length in [1, 2, 3, 4, 5, 6, 7, 8, 59, 60]:
        for _ in range(n(10)):
            cases.append((rand_seq(rng, length), rand_args(rng)))
    # Designed stem-loops hitting every triloop/tetraloop bonus entry.
    for key in loop_keys("triloop.ds") + loop_keys("tetraloop.ds"):
        for _ in range(n(8)):
            cases.append((stem_loop(rng, key), rand_args(rng)))
    for _ in range(n(1200)):
        loop = rand_seq(rng, rng.randint(3, 30))
        cases.append((stem_loop(rng, loop), rand_args(rng)))
    # Primer-like sequences at realistic and default conditions.
    for _ in range(n(1000)):
        s = rand_seq(rng, rng.randint(18, 32), rng.uniform(0.4, 0.6))
        cases.append((s, default_args() if rng.random() < 0.3 else rand_args(rng)))
    # Low complexity: homopolymers, dinucleotide and short tandem repeats.
    for _ in range(n(300)):
        unit = rand_seq(rng, rng.randint(1, 4))
        cases.append(((unit * 60)[: rng.randint(4, 60)], rand_args(rng)))
    # Lowercase, mixed case, N and other non-ACGT characters.
    for _ in range(n(300)):
        s = list(stem_loop(rng, rand_seq(rng, rng.randint(3, 8))))
        for i in range(len(s)):
            r = rng.random()
            if r < 0.3:
                s[i] = s[i].lower()
            elif r < 0.35:
                s[i] = rng.choice("NnRYX")
        cases.append(("".join(s), rand_args(rng)))
    # Structure strings.
    for _ in range(n(300)):
        cases.append((stem_loop(rng, rand_seq(rng, rng.randint(3, 12))), rand_args(rng, True)))
    # Inputs primer3 rejects.
    errs = [("", default_args()), (rand_seq(rng, 61), default_args()), (rand_seq(rng, 200), default_args())]
    for ml in [-1, 31, 100]:
        a = default_args()
        a["max_loop"] = ml
        errs.append((rand_seq(rng, 40), a))
    cases.extend(errs)
    return [{"seq": s, "args": a, **run(primer3.calc_hairpin, [s], a)} for s, a in cases]


def write(name, records):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.jsonl"
    with path.open("w") as f:
        f.write(json.dumps({"primer3_py_version": primer3.__version__, "cases": len(records)}) + "\n")
        for r in records:
            f.write(json.dumps(r, allow_nan=False) + "\n")
    print(f"{path.relative_to(ROOT)}: {len(records)} cases")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--scale", type=float, default=1.0)
    a = p.parse_args()
    write("hairpin", hairpin_cases(random.Random(a.seed), a.scale))


if __name__ == "__main__":
    main()
