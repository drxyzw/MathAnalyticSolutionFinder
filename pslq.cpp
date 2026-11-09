#include "pslq.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace
{
	using Matrix = std::vector<std::vector<mpreal>>;

	// Hermite-reduce H in place, applying the same integer operations to y, A and B.
	// Rows i in [iBegin, n), columns j from min(i - 1, jMax) down to 0.
	void reduce(Matrix& H, Matrix& A, Matrix& B, std::vector<mpreal>& y, size_t iBegin, size_t jMax)
	{
		const size_t n = y.size();
		for (size_t i = iBegin; i < n; ++i)
		{
			const size_t jStart = std::min(i - 1, jMax);
			for (size_t jj = jStart + 1; jj-- > 0;)
			{
				const mpreal t = mpfr::round(H[i][jj] / H[jj][jj]);
				if (t == 0)
					continue;
				y[jj] += t * y[i];
				for (size_t k = 0; k <= jj; ++k)
					H[i][k] -= t * H[jj][k];
				for (size_t k = 0; k < n; ++k)
				{
					A[i][k] -= t * A[jj][k];
					B[k][jj] += t * B[k][i];
				}
			}
		}
	}
}

const char* toString(PslqResult::Status s)
{
	switch (s)
	{
	case PslqResult::Status::Found: return "Found";
	case PslqResult::Status::NormBoundReached: return "NormBoundReached";
	case PslqResult::Status::PrecisionExhausted: return "PrecisionExhausted";
	case PslqResult::Status::MaxIterReached: return "MaxIterReached";
	}
	return "?";
}

PslqResult pslq(const std::vector<mpreal>& x, const mpreal& maxNormBound, int maxIter)
{
	const size_t n = x.size();
	if (n < 2)
		throw std::invalid_argument("pslq: need at least 2 numbers");

	const mp_prec_t prec = mpreal::get_default_prec();
	const double digits = mpfr::bits2digits(prec);
	// Fraction of the working digits regarded as reliable.
	const double reliable = 0.8 * digits;
	// |y_j| below this is treated as zero (relation detected).
	const mpreal detectEps = mpfr::pow(mpreal(10), -reliable);
	// Entries of A above this mean the precision is used up.
	const mpreal maxA = mpfr::pow(mpreal(10), reliable);
	// A generic (spurious) relation of norm N satisfies |m.x| ~ N^-(n-1), so relations with
	// N^n > 10^reliable cannot be told apart from numerical noise.
	const mpreal maxMeaningfulNorm = mpfr::pow(mpreal(10), reliable / n);
	const mpreal gamma = mpfr::sqrt(mpreal(4) / 3) * mpreal(1.0001); // gamma > sqrt(4/3)

	PslqResult result;

	// Step 1: normalize x and compute partial norms s_k = |(x_k, ..., x_{n-1})|.
	std::vector<mpreal> s(n);
	{
		mpreal acc = 0;
		for (size_t k = n; k-- > 0;)
		{
			acc += x[k] * x[k];
			s[k] = mpfr::sqrt(acc);
		}
	}
	if (s[0] == 0)
		throw std::invalid_argument("pslq: zero vector");
	std::vector<mpreal> y(n);
	const mpreal norm = s[0];
	for (size_t k = 0; k < n; ++k)
	{
		y[k] = x[k] / norm;
		s[k] /= norm;
	}
	for (size_t k = 0; k < n; ++k)
	{
		if (y[k] == 0)
		{
			// trivial relation: e_k
			result.status = PslqResult::Status::Found;
			result.relation.assign(n, mpreal(0));
			result.relation[k] = 1;
			result.residual = 0;
			return result;
		}
	}

	// Step 2: H is n x (n-1) lower trapezoidal.
	Matrix H(n, std::vector<mpreal>(n - 1, mpreal(0)));
	for (size_t j = 0; j + 1 < n; ++j)
	{
		H[j][j] = s[j + 1] / s[j];
		for (size_t i = j + 1; i < n; ++i)
			H[i][j] = -y[i] * y[j] / (s[j] * s[j + 1]);
	}

	// Step 3: A = B = I (A * B = I, B columns hold relation candidates).
	Matrix A(n, std::vector<mpreal>(n, mpreal(0)));
	Matrix B(n, std::vector<mpreal>(n, mpreal(0)));
	for (size_t i = 0; i < n; ++i)
		A[i][i] = B[i][i] = 1;

	// Step 4: full reduction.
	reduce(H, A, B, y, 1, n - 2);

	for (int iter = 1; iter <= maxIter; ++iter)
	{
		result.iterations = iter;

		// Step 1: select m maximizing gamma^(i+1) |H_ii|.
		size_t m = 0;
		{
			mpreal best = -1;
			mpreal g = gamma;
			for (size_t i = 0; i + 1 < n; ++i)
			{
				const mpreal v = g * mpfr::abs(H[i][i]);
				if (v > best)
				{
					best = v;
					m = i;
				}
				g *= gamma;
			}
		}

		// Step 2: exchange entries m and m+1.
		std::swap(y[m], y[m + 1]);
		std::swap(A[m], A[m + 1]);
		std::swap(H[m], H[m + 1]);
		for (size_t k = 0; k < n; ++k)
			std::swap(B[k][m], B[k][m + 1]);

		// Step 3: corner removal (restore lower trapezoidal shape).
		if (m + 2 < n)
		{
			const mpreal t0 = mpfr::sqrt(H[m][m] * H[m][m] + H[m][m + 1] * H[m][m + 1]);
			const mpreal t1 = H[m][m] / t0;
			const mpreal t2 = H[m][m + 1] / t0;
			for (size_t i = m; i < n; ++i)
			{
				const mpreal t3 = H[i][m];
				const mpreal t4 = H[i][m + 1];
				H[i][m] = t1 * t3 + t2 * t4;
				H[i][m + 1] = -t2 * t3 + t1 * t4;
			}
		}

		// Step 4: reduction.
		reduce(H, A, B, y, m + 1, m + 1);

		// Step 5: norm bound. Any relation has |m| >= 1 / max |H_jj|.
		mpreal maxH = 0;
		for (size_t j = 0; j + 1 < n; ++j)
			maxH = mpfr::max(maxH, mpfr::abs(H[j][j]));
		if (maxH > 0)
		{
			const mpreal bound = 1 / maxH;
			if (bound > result.normBound)
				result.normBound = bound;
		}

		// Step 6: termination tests.
		size_t jMin = 0;
		for (size_t j = 1; j < n; ++j)
			if (mpfr::abs(y[j]) < mpfr::abs(y[jMin]))
				jMin = j;
		if (mpfr::abs(y[jMin]) < detectEps || maxH == 0)
		{
			std::vector<mpreal> rel(n);
			mpreal norm2 = 0, r = 0;
			for (size_t k = 0; k < n; ++k)
			{
				rel[k] = B[k][jMin];
				norm2 += rel[k] * rel[k];
				r += rel[k] * x[k];
			}
			if (mpfr::sqrt(norm2) > maxMeaningfulNorm)
			{
				// tiny y_j caused by precision loss, not a genuine relation
				result.status = PslqResult::Status::PrecisionExhausted;
				return result;
			}
			result.status = PslqResult::Status::Found;
			result.relation = std::move(rel);
			result.residual = mpfr::abs(r);
			return result;
		}

		if (result.normBound >= maxMeaningfulNorm)
		{
			result.status = PslqResult::Status::PrecisionExhausted;
			return result;
		}

		if (maxNormBound > 0 && result.normBound >= maxNormBound)
		{
			result.status = PslqResult::Status::NormBoundReached;
			return result;
		}

		mpreal maxAbsA = 0;
		for (size_t i = 0; i < n; ++i)
			for (size_t k = 0; k < n; ++k)
				maxAbsA = mpfr::max(maxAbsA, mpfr::abs(A[i][k]));
		if (maxAbsA > maxA)
		{
			result.status = PslqResult::Status::PrecisionExhausted;
			return result;
		}
	}

	result.status = PslqResult::Status::MaxIterReached;
	return result;
}

PslqResult findIntegerPolynomial(const mpreal& x, int degree, const mpreal& maxNormBound, int maxIter)
{
	if (degree < 1)
		throw std::invalid_argument("findIntegerPolynomial: degree must be >= 1");
	std::vector<mpreal> powers(degree + 1);
	powers[0] = 1;
	for (int k = 1; k <= degree; ++k)
		powers[k] = powers[k - 1] * x;
	return pslq(powers, maxNormBound, maxIter);
}
