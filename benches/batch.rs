//! Mean time per call on primer-sized (18-30 nt) sequences, the main use in
//! primalschemers. Each iteration is one independent call on the next input.
use criterion::{Criterion, black_box, criterion_group, criterion_main};
#[cfg(unix)]
use pprof::criterion::{Output, PProfProfiler};

use primer3_thal_rs::{ThalArgs, calc_hairpin, calc_heterodimer, calc_homodimer};

/// Small deterministic LCG so the workload is stable across runs.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as usize
    }

    fn seq(&mut self, len: usize) -> String {
        (0..len).map(|_| b"ACGT"[self.next() % 4] as char).collect()
    }

    fn primer(&mut self) -> String {
        let len = 18 + self.next() % 13;
        self.seq(len)
    }
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}

fn primers(n: usize) -> Vec<String> {
    let mut rng = Lcg(0x2545_f491_4f6c_dd1d);
    (0..n).map(|_| rng.primer()).collect()
}

/// Representative primer pairs: 80% two independent random primers, 20% a
/// partner whose 3' end is complementary to the primer's 3' end (4-10 nt),
/// the case dimer checks exist to catch.
fn primer_pairs(n: usize) -> Vec<(String, String)> {
    let mut rng = Lcg(0x9e37_79b9_7f4a_7c15);
    (0..n)
        .map(|_| {
            let a = rng.primer();
            let b = if rng.next().is_multiple_of(5) {
                let k = 4 + rng.next() % 7;
                let tail_len = 18 + rng.next() % 13 - k;
                let mut b = rng.seq(tail_len);
                b.push_str(&revcomp(&a[a.len() - k..]));
                b
            } else {
                rng.primer()
            };
            (a, b)
        })
        .collect()
}

fn bench_per_call(c: &mut Criterion) {
    let args = ThalArgs::default();
    let seqs = primers(1000);
    let pairs = primer_pairs(1000);

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
            k = (k + 1) % pairs.len();
            let (a, p) = &pairs[k];
            calc_heterodimer(black_box(a), black_box(p), &args, false)
        })
    });
    group.finish();
}

/// Criterion with flamegraph profiling (`--profile-time`) where supported.
fn config() -> Criterion {
    let c = Criterion::default();
    #[cfg(unix)]
    let c = c.with_profiler(PProfProfiler::new(1000, Output::Flamegraph(None)));
    c
}

criterion_group!(
    name = benches;
    config = config();
    targets = bench_per_call
);
criterion_main!(benches);
