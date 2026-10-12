use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use pprof::criterion::{Output, PProfProfiler};

use primer3_thal_rs::{ThalArgs, calc_heterodimer, calc_homodimer};

fn bench_calc_heterodimer(c: &mut Criterion) {
    let args = ThalArgs::default();
    let pairs = [
        ("CCCCCATCCGATCAGGGGG", "CCCCCATCCGATCAGGGGG"),
        (
            "GTGTCCAATCGACATTTAGGAGCTACGCGGCCAGGGGCAAA",
            "TGCATTGAATGCGTGTCACGTTATGCACGC",
        ),
        ("GATTACAGATTACAGATTACA", "ACGTACGTACGTACGTACGT"),
        (
            "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCG",
            "GCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTA",
        ),
    ];

    let mut group = c.benchmark_group("calc_heterodimer");
    for (seq1, seq2) in pairs {
        let id = format!("{}x{}", seq1.len(), seq2.len());
        group.bench_with_input(BenchmarkId::from_parameter(id), &(seq1, seq2), |b, s| {
            b.iter(|| {
                let res = calc_heterodimer(
                    black_box(s.0),
                    black_box(s.1),
                    black_box(&args),
                    black_box(false),
                );
                let _ = black_box(res);
            })
        });
    }
    group.finish();
}

fn bench_calc_homodimer(c: &mut Criterion) {
    let args = ThalArgs::default();
    let seqs = [
        "CCCCCATCCGATCAGGGGG",
        "TGCATTGAATGCGTGTCACGTTATGCACGC",
        "GTGTCCAATCGACATTTAGGAGCTACGCGGCCAGGGGCAAA",
        "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCGTCCGATCGTACGTAGCTAGC",
    ];

    let mut group = c.benchmark_group("calc_homodimer");
    for seq in seqs {
        group.bench_with_input(BenchmarkId::from_parameter(seq.len()), seq, |b, s| {
            b.iter(|| black_box(calc_homodimer(black_box(s), black_box(&args), false)))
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
    targets = bench_calc_heterodimer, bench_calc_homodimer
);
criterion_main!(benches);
