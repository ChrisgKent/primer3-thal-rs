use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
#[cfg(unix)]
use pprof::criterion::{Output, PProfProfiler};

use primer3_thal_rs::{ThalArgs, nn};

const SEQS: [&str; 3] = [
    "GATTACAGATTACAGATT",
    "ACCAACGATGGTGTGTCCAT",
    "TCATTGTATCCTCACATAACTCTCCCAAAT",
];

fn bench_nn(c: &mut Criterion) {
    let args = ThalArgs {
        mv: 100.0,
        dv: 2.0,
        dntp: 0.8,
        dna_conc: 15.0,
        ..ThalArgs::default()
    };

    // The digest walk: extend by one base and check Tm at every step
    let mut group = c.benchmark_group("duplex_walk");
    for seq in SEQS {
        group.bench_with_input(BenchmarkId::from_parameter(seq.len()), seq, |b, s| {
            b.iter(|| {
                let mut d = nn::Duplex::new();
                let mut tm = 0.0;
                for &base in black_box(s.as_bytes()) {
                    d.push(base).unwrap();
                    tm += d.tm(&args, false);
                }
                black_box(tm)
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("duplex_tm");
    for seq in SEQS {
        group.bench_with_input(BenchmarkId::from_parameter(seq.len()), seq, |b, s| {
            b.iter(|| black_box(nn::duplex_tm(black_box(s), &args).unwrap()))
        });
    }
    group.finish();

    c.bench_function("bound_percent", |b| {
        b.iter(|| {
            black_box(nn::bound_percent(
                black_box(-160_000.0),
                black_box(-430.0),
                &args,
                false,
                black_box(65.0),
            ))
        })
    });
}

/// Criterion with flamegraph profiling (`--profile-time`) where supported.
fn config() -> Criterion {
    let c = Criterion::default();
    #[cfg(unix)]
    let c = c.with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
    c
}

criterion_group!(
    name = benches;
    config = config();
    targets = bench_nn
);
criterion_main!(benches);
