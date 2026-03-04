use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pprof::criterion::{Output, PProfProfiler};

use primer3_py_rs::{calc_hairpin, ThalArgs};

fn bench_calc_hairpin(c: &mut Criterion) {
    let args = ThalArgs::hairpin_defaults();
    let seqs = [
        "CCCCCATCCGATCAGGGGG",
        "GTGTCCAATCGACATTTAGGAGCTACGCGGCCAGGGGCAAA",
        "TGCATTGAATGCGTGTCACGTTATGCACGC",
        "GATTACAGATTACAGATTACA",
        "ACGTACGTACGTACGTACGT",
        "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCG",
    ];

    let mut group = c.benchmark_group("calc_hairpin");
    for seq in seqs {
        group.bench_with_input(BenchmarkId::from_parameter(seq.len()), seq, |b, s| {
            b.iter(|| {
                let res = calc_hairpin(black_box(s), black_box(&args), black_box(false));
                black_box(res);
            })
        });
    }
    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(
        100,
        Output::Flamegraph(None),
    ));
    targets = bench_calc_hairpin
);
criterion_main!(benches);
