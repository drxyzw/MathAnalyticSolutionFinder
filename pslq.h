#pragma once

#include <vector>

#if __has_include(<gmp_util.h>)
#include <gmp_util.h>
#else
#include <mpreal.h>
#endif

using mpfr::mpreal;

// Result of the PSLQ integer relation search.
struct PslqResult
{
	enum class Status
	{
		Found,              // integer relation found: sum relation[i] * x[i] ~ 0
		NormBoundReached,   // no relation with Euclidean norm < normBound exists
		PrecisionExhausted, // working precision is used up; normBound is still valid
		MaxIterReached,     // iteration limit; normBound is still valid
	};

	Status status = Status::MaxIterReached;
	std::vector<mpreal> relation; // integer coefficients (valid when status == Found)
	mpreal normBound;             // any integer relation m != 0 of x satisfies |m|_2 >= normBound
	mpreal residual;              // |sum relation[i] * x[i]| (valid when status == Found)
	int iterations = 0;
};

// PSLQ algorithm (Ferguson-Bailey-Arno 1999) in the current mpreal default precision.
// Searches for an integer vector m != 0 with sum m[i] * x[i] = 0.
// - maxNormBound: stop (NormBoundReached) once it is proven that no relation of norm < maxNormBound exists.
//                 Pass 0 to disable.
// - maxIter:      iteration limit.
PslqResult pslq(const std::vector<mpreal>& x, const mpreal& maxNormBound = 0, int maxIter = 1000000);

// Algorithm selection for pslq() / findIntegerPolynomial().
//   Standard    : original single-level PSLQ (default).
//   Multipair2L : Bailey's two-level multipair PSLQ (double-precision inner iterations).
// The initial value is taken once from the env var PSLQ_ALGO ("m2l" -> Multipair2L, anything else -> Standard);
// setPslqAlgorithm() overrides it. The setting is process-global and thread-safe (atomic).
enum class PslqAlgo { Standard, Multipair2L };
void setPslqAlgorithm(PslqAlgo a);
PslqAlgo getPslqAlgorithm();

// Explicit entry points (independent of the global switch).
PslqResult pslqStandard(const std::vector<mpreal>& x, const mpreal& maxNormBound = 0, int maxIter = 1000000);
PslqResult pslqMultipair2L(const std::vector<mpreal>& x, const mpreal& maxNormBound = 0, int maxIter = 1000000);
// Number of times the debug invariant check (env PSLQ_M2L_VERIFY) failed in this process.
long pslqMultipair2LInvariantViolations();

// Searches for an integer polynomial p of degree <= degree with p(x) = 0,
// by running PSLQ on (1, x, x^2, ..., x^degree).
// relation[k] is the coefficient of x^k.
PslqResult findIntegerPolynomial(const mpreal& x, int degree, const mpreal& maxNormBound = 0, int maxIter = 1000000);

const char* toString(PslqResult::Status s);
