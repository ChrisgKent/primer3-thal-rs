//! Primer-sized workload, the main use in primalschemers: many 18-30 nt
//! sequences through hairpin, homodimer and heterodimer, single-threaded.
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

fn bench_batch(c: &mut Criterion) {
    let args = ThalArgs::default();
    let seqs = primers(1000);
    let partner = "GTAAAACGACGGCCAGTGAATTC";

    let mut group = c.benchmark_group("batch_1000_primers");
    group.sample_size(20);
    group.bench_function("hairpin", |b| {
        b.iter(|| {
            for s in &seqs {
                black_box(calc_hairpin(black_box(s), &args, false).ok());
            }
        })
    });
    group.bench_function("homodimer", |b| {
        b.iter(|| {
            for s in &seqs {
                black_box(calc_homodimer(black_box(s), &args, false).ok());
            }
        })
    });
    group.bench_function("heterodimer", |b| {
        b.iter(|| {
            for s in &seqs {
                black_box(calc_heterodimer(black_box(s), partner, &args, false).ok());
            }
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
    targets = bench_batch
);
criterion_main!(benches);
