// Copyright (c) 2022 Espresso Systems (espressosys.com)
// This file is part of the Jellyfish library.

// You should have received a copy of the MIT License
// along with the Jellyfish library. If not, see <https://mit-license.org/>.

use std::time::{Duration, Instant};

use ark_bls12_381::Bls12_381;
use ark_bn254::Bn254;
use ark_ec::pairing::Pairing;
use ark_ff::{UniformRand, PrimeField, BigInteger};
use ark_poly::DenseMultilinearExtension;
use ark_std::rand::Rng;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use jf_primitives::pcs::{
    prelude::{MultilinearKzgPCS, PolynomialCommitmentScheme, MLE},
    StructuredReferenceString,
    
};
use jf_utils::test_rng;

use ark_ec::scalar_mul::variable_base::SMALLNESS; // coefficient bit size that's considered small

const MIN_NUM_VARS: usize = 10;
const MAX_NUM_VARS: usize = 20;

/// Produce a random MLE with small coefficients
fn small_mle<F: PrimeField>(
    num_vars: usize,
    rng: &mut impl Rng,
) -> DenseMultilinearExtension<F> {
    let num_bits = F::MODULUS_BIT_SIZE as usize;
    
    let small_scalars = (0..(1 << num_vars))
        .map(|_| {
            let s = F::rand(rng).into_bigint();
            let mut bits = s.to_bits_le();
            bits.truncate(SMALLNESS);
            bits.resize(num_bits, false);
            let bigint = F::BigInt::from_bits_le(&bits);
            F::from_bigint(bigint).unwrap()
        })
        .collect::<Vec<_>>();

    DenseMultilinearExtension::from_evaluations_vec(num_vars, small_scalars)
}

/// Measure the time cost of {commit/open/verify} across a range of num_vars
pub fn bench_pcs_method<E: Pairing>(
    c: &mut Criterion,
    range: impl Iterator<Item = usize>,
    msg: &str,
    method: impl Fn(&<MultilinearKzgPCS<E> as PolynomialCommitmentScheme>::SRS, usize) -> Duration,
) {
    let mut group = c.benchmark_group(msg);

    let mut rng = &mut test_rng();

    for num_vars in range {
        let pp = MultilinearKzgPCS::<E>::gen_srs_for_testing(&mut rng, num_vars).unwrap();

        group.bench_with_input(
            BenchmarkId::from_parameter(num_vars),
            &num_vars,
            |b, num_vars| {
                b.iter_custom(|i| {
                    let mut time = Duration::from_nanos(0);
                    for _ in 0..i {
                        time += method(&pp, *num_vars);
                    }
                    time
                });
            },
        );
    }

    group.finish();
}

/// Report the time cost of a commitment
pub fn commit<E: Pairing>(
    pp: &<MultilinearKzgPCS<E> as PolynomialCommitmentScheme>::SRS,
    num_vars: usize,
) -> Duration {
    let rng = &mut test_rng();

    let (ml_ck, _ml_vk) = pp.0.trim(num_vars).unwrap();
    let (uni_ck, _uni_vk) = pp.1.trim(num_vars).unwrap();
    let ck = (ml_ck, uni_ck);

    let poly = MLE::from(small_mle(num_vars, rng));

    let start = Instant::now();
    let _ = MultilinearKzgPCS::commit(&ck, &poly).unwrap();
    start.elapsed()
}

/// Report the time cost of an opening
pub fn open<E: Pairing>(
    pp: &<MultilinearKzgPCS<E> as PolynomialCommitmentScheme>::SRS,
    num_vars: usize,
) -> Duration {
    let rng = &mut test_rng();

    let (ml_ck, _ml_vk) = pp.0.trim(num_vars).unwrap();
    let (uni_ck, _uni_vk) = pp.1.trim(num_vars).unwrap();
    let ck = (ml_ck, uni_ck);

    let poly = MLE::from(small_mle(num_vars, rng));
    let point: Vec<_> = (0..num_vars).map(|_| E::ScalarField::rand(rng)).collect();

    let start = Instant::now();
    let _ = MultilinearKzgPCS::open(&ck, &poly, &point).unwrap();
    start.elapsed()
}

/// Report the time cost of a verification
pub fn verify<E: Pairing>(
    pp: &<MultilinearKzgPCS<E> as PolynomialCommitmentScheme>::SRS,
    num_vars: usize,
) -> Duration {
    let rng = &mut test_rng();

    let (ml_ck, ml_vk) = pp.0.trim(num_vars).unwrap();
    let (uni_ck, uni_vk) = pp.1.trim(num_vars).unwrap();
    let ck = (ml_ck, uni_ck);
    let vk = (ml_vk, uni_vk);

    let poly = MLE::from(small_mle(num_vars, rng));
    let point: Vec<_> = (0..num_vars).map(|_| E::ScalarField::rand(rng)).collect();

    let commitment = MultilinearKzgPCS::commit(&ck, &poly).unwrap();

    let (proof, value) = MultilinearKzgPCS::open(&ck, &poly, &point).unwrap();

    let start = Instant::now();
    assert!(MultilinearKzgPCS::verify(&vk, &commitment, &point, &value, &proof).unwrap());
    start.elapsed()
}

fn kzg_254(c: &mut Criterion) {
    bench_pcs_method::<Bn254>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2),
        "commit_kzg_range_BN_254",
        commit::<Bn254>,
    );
    bench_pcs_method::<Bn254>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2),
        "open_kzg_range_BN_254",
        open::<Bn254>,
    );
    bench_pcs_method::<Bn254>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2),
        "verify_kzg_range_BN_254",
        verify::<Bn254>,
    );
}

fn kzg_381(c: &mut Criterion) {
    bench_pcs_method::<Bls12_381>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2),
        "commit_kzg_range_BLS_381",
        commit::<Bls12_381>,
    );
    bench_pcs_method::<Bls12_381>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2),
        "open_kzg_range_BLS_381",
        open::<Bls12_381>,
    );
    bench_pcs_method::<Bls12_381>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(2),
        "verify_kzg_range_BLS_381",
        verify::<Bls12_381>,
    );
}

criterion_group! {
    name = pcs_benches;
    config = Criterion::default();
    targets =
        kzg_254,
        kzg_381
}

criterion_main!(pcs_benches);
