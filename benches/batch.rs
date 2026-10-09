//! Mean time per call on primer-sized (18-30 nt) sequences, the main use in
//! primalschemers.
use criterion::{Criterion, black_box, criterion_group, criterion_main};
use pprof::criterion::{Output, PProfProfiler};

use primer3_py_rs::{ThalArgs, calc_hairpin, calc_heterodimer, calc_homodimer};

fn primers(n: usize) -> Vec<String> {
    // Small deterministic LCG so the workload is stable across runs.
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    (0..n)
        .map(|_| {
            let len = 18 + next() % 13;
            (0..len).map(|_| b"ACGT"[next() % 4] as char).collect()
        })
        .collect()
}

fn bench_per_call(c: &mut Criterion) {
    let args = ThalArgs::default();
    let seqs = primers(1000);
    let partner = "GTAAAACGACGGCCAGTGAATTC";

    // Each iteration is one independent call on the next primer, so the
    // reported time is the mean per call over realistic 18-30 nt primers.
    let mut group = c.benchmark_group("per_call_18_30nt");
    group.bench_function("hairpin", |b| {
        let mut k = 0;
        b.iter(|| {
            k = (k + 1) % seqs.len();
            calc_hairpin(black_box(&seqs[k]), &args, false)
        })
    });
    group.bench_function("homodimer", |b| {
        let mut k = 0;
        b.iter(|| {
            k = (k + 1) % seqs.len();
            calc_homodimer(black_box(&seqs[k]), &args, false)
        })
    });
    group.bench_function("heterodimer", |b| {
        let mut k = 0;
        b.iter(|| {
            k = (k + 1) % seqs.len();
            calc_heterodimer(black_box(&seqs[k]), partner, &args, false)
        })
    });
    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(
        1000,
        Output::Flamegraph(None),
    ));
    targets = bench_per_call
);
criterion_main!(benches);
