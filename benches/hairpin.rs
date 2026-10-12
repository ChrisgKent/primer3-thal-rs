use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use pprof::criterion::{Output, PProfProfiler};

use primer3_thal_rs::{ThalArgs, calc_hairpin};

fn bench_calc_hairpin(c: &mut Criterion) {
    let args = ThalArgs::default();
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
                let _ = black_box(res);
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
