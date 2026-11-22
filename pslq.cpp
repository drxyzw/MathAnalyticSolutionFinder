#include "pslq.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <type_traits>
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

PslqResult pslqStandard(const std::vector<mpreal>& x, const mpreal& maxNormBound, int maxIter)
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

// ---------------------------------------------------------------------------------------------
// Algorithm switch
// ---------------------------------------------------------------------------------------------
namespace
{
	std::atomic<int>& algoFlag()
	{
		static std::atomic<int> flag([]
			{
				const char* e = std::getenv("PSLQ_ALGO");
				return (e && (std::strcmp(e, "m2l") == 0 || std::strcmp(e, "M2L") == 0)) ? 1 : 0;
			}());
		return flag;
	}
}

void setPslqAlgorithm(PslqAlgo a) { algoFlag().store(a == PslqAlgo::Multipair2L ? 1 : 0); }
PslqAlgo getPslqAlgorithm() { return algoFlag().load() == 1 ? PslqAlgo::Multipair2L : PslqAlgo::Standard; }

PslqResult pslq(const std::vector<mpreal>& x, const mpreal& maxNormBound, int maxIter)
{
	if (getPslqAlgorithm() == PslqAlgo::Multipair2L)
		return pslqMultipair2L(x, maxNormBound, maxIter);
	return pslqStandard(x, maxNormBound, maxIter);
}

// ---------------------------------------------------------------------------------------------
// Two-level multipair PSLQ (D. H. Bailey, "The two-level multipair PSLQ algorithm", 2025)
//
// Multiprecision state: y (normalized |y|=1 initially, y = B^T x / |x|), B (unimodular, columns are
// relation candidates), H (n x (n-1), kept lower trapezoidal; H = A H0 Q with A = B^-1 unimodular, Q orthogonal).
// Inner loop: many multipair iterations in double precision on copies dH, dy with integer matrices dA, dB
// (starting from the identity). Those integer matrices are exact as long as all entries stay below 2^50
// (checked after every elementary update; on violation the iteration is reverted). Afterwards
//   y := dB^T y,  B := B dB,  H := dA H,  followed by an LQ decomposition of H (multiprecision).
// The double-precision data only decides WHICH unimodular operations are applied; all invariants
// (y = B^T x, H = A H0 Q) are maintained exactly in multiprecision, so the relation and the norm bound
// 1/max|H_jj| (computed from the multiprecision H) do not depend on the double-precision arithmetic.
// ---------------------------------------------------------------------------------------------
namespace
{
	using MMat = std::vector<std::vector<mpreal>>;
	using DMat = std::vector<std::vector<double>>;

	constexpr double kHardLimit = 1125899906842624.0; // 2^50: any entry above this => precision failure
	constexpr double kSoftLimit = 1e13;               // stop the double-precision batch above this
	constexpr double kDepsRel = 1e-12;                // stop the batch if min|dy| < kDepsRel * max|dy|
	constexpr double kDrep = 1e-10;                   // below this dynamic range of y use multiprecision iterations
	constexpr int kMpBlock = 10;                      // multiprecision iterations per block
	constexpr int kRing = 8;                          // length of the y history for cycle detection
	constexpr double kBeta = 0.4;

	inline double rnd(double v) { return std::nearbyint(v); }
	inline mpreal rnd(const mpreal& v) { return mpfr::round(v); }
	inline double absv(double v) { return std::fabs(v); }
	inline mpreal absv(const mpreal& v) { return mpfr::abs(v); }
	inline double sqrtv(double v) { return std::sqrt(v); }
	inline mpreal sqrtv(const mpreal& v) { return mpfr::sqrt(v); }

	// One multipair PSLQ iteration (selection, exchange, corner removal, Hermite reduction).
	// A may be null (multiprecision level does not need it). For T = double returns false on a
	// precision failure (entry of A/B above kHardLimit or a non-finite quantity); the arrays are
	// then in an undefined state and the caller must restore them.
	template <class T>
	bool multipairIteration(std::vector<std::vector<T>>& H, std::vector<std::vector<T>>* A,
		std::vector<std::vector<T>>& B, std::vector<T>& y, int maxPairs, const T& gamma, double* maxAB)
	{
		constexpr bool isDouble = std::is_same<T, double>::value;
		const size_t n = y.size();

		// Selection: sort gamma^(i+1) |H_ii| in decreasing order, pick disjoint pairs (m, m+1).
		std::vector<std::pair<T, size_t>> keys;
		keys.reserve(n - 1);
		{
			T g = gamma;
			for (size_t i = 0; i + 1 < n; ++i)
			{
				keys.emplace_back(g * absv(H[i][i]), i);
				g *= gamma;
			}
		}
		std::stable_sort(keys.begin(), keys.end(),
			[](const std::pair<T, size_t>& a, const std::pair<T, size_t>& b) { return a.first > b.first; });
		std::vector<char> used(n, 0);
		std::vector<size_t> pairs;
		for (const auto& kv : keys)
		{
			if ((int)pairs.size() >= maxPairs)
				break;
			const size_t m = kv.second;
			if (used[m] || used[m + 1])
				continue;
			used[m] = used[m + 1] = 1;
			pairs.push_back(m);
		}
		if (pairs.empty())
			return false;

		// Exchange.
		size_t mMin = n, mMax = 0;
		for (size_t m : pairs)
		{
			mMin = std::min(mMin, m);
			mMax = std::max(mMax, m);
			std::swap(y[m], y[m + 1]);
			if (A) std::swap((*A)[m], (*A)[m + 1]);
			std::swap(H[m], H[m + 1]);
			for (size_t k = 0; k < n; ++k)
				std::swap(B[k][m], B[k][m + 1]);
		}

		// Corner removal.
		for (size_t m : pairs)
		{
			if (m + 2 >= n)
				continue;
			const T t0 = sqrtv(H[m][m] * H[m][m] + H[m][m + 1] * H[m][m + 1]);
			if (t0 == 0)
				continue;
			const T t1 = H[m][m] / t0;
			const T t2 = H[m][m + 1] / t0;
			for (size_t i = m; i < n; ++i)
			{
				const T t3 = H[i][m];
				const T t4 = H[i][m + 1];
				H[i][m] = t1 * t3 + t2 * t4;
				H[i][m + 1] = -t2 * t3 + t1 * t4;
			}
		}
		if constexpr (isDouble)
		{
			for (size_t m : pairs)
				if (!std::isfinite(H[m][m]) || !std::isfinite(H[m][m + 1]))
					return false;
		}

		// Hermite reduction of rows mMin+1 .. n-1, columns up to mMax+1.
		for (size_t i = mMin + 1; i < n; ++i)
		{
			const size_t jStart = std::min(i - 1, mMax + 1);
			for (size_t jj = jStart + 1; jj-- > 0;)
			{
				if (H[jj][jj] == 0)
					continue;
				const T t = rnd(H[i][jj] / H[jj][jj]);
				if constexpr (isDouble)
				{
					if (!std::isfinite(t) || std::fabs(t) > kHardLimit)
						return false;
				}
				if (t == 0)
					continue;
				y[jj] += t * y[i];
				for (size_t k = 0; k <= jj; ++k)
					H[i][k] -= t * H[jj][k];
				if constexpr (isDouble)
				{
					if (!std::isfinite(y[jj]))
						return false;
					for (size_t k = 0; k < n; ++k)
					{
						const double a = (*A)[i][k] - t * (*A)[jj][k];
						const double b = B[k][jj] + t * B[k][i];
						if (!(std::fabs(a) <= kHardLimit) || !(std::fabs(b) <= kHardLimit))
							return false;
						(*A)[i][k] = a;
						B[k][jj] = b;
						*maxAB = std::max(*maxAB, std::max(std::fabs(a), std::fabs(b)));
					}
				}
				else
				{
					for (size_t k = 0; k < n; ++k)
					{
						if (A) (*A)[i][k] -= t * (*A)[jj][k];
						B[k][jj] += t * B[k][i];
					}
				}
			}
		}
		return true;
	}

	// LQ decomposition (Householder, right-multiplication by orthogonal matrices):
	// M (n x (n-1)) -> lower trapezoidal. Row l is reduced by a reflection acting on columns l..m-1.
	void lqDecompose(MMat& M)
	{
		const size_t n = M.size();
		const size_t m = n - 1;
		std::vector<mpreal> v(m);
		mpreal N2, N, alpha, vtv, dot, f;
		for (size_t l = 0; l + 1 < m; ++l)
		{
			N2 = 0;
			for (size_t c = l; c < m; ++c)
				N2 += M[l][c] * M[l][c];
			if (N2 == 0)
				continue;
			N = mpfr::sqrt(N2);
			alpha = (M[l][l] >= 0) ? N : -N;
			for (size_t c = l; c < m; ++c)
				v[c] = M[l][c];
			v[l] += alpha;
			vtv = 2 * alpha * (alpha + M[l][l]);
			if (vtv == 0)
				continue;
			for (size_t j = l + 1; j < n; ++j)
			{
				dot = 0;
				for (size_t c = l; c < m; ++c)
					dot += M[j][c] * v[c];
				f = 2 * dot / vtv;
				for (size_t c = l; c < m; ++c)
					M[j][c] -= f * v[c];
			}
			M[l][l] = -alpha;
			for (size_t c = l + 1; c < m; ++c)
				M[l][c] = 0;
		}
	}

	using SparseCols = std::vector<std::vector<std::pair<int, double>>>;

	// y := dB^T y, B := B dB, H := dA H (all exact integer matrices), then LQ of H.
	void applyBatch(MMat& H, MMat& B, std::vector<mpreal>& y, const DMat& dA, const DMat& dB)
	{
		const size_t n = y.size();
		SparseCols colsB(n); // colsB[j] = nonzeros (l, dB[l][j])
		for (size_t l = 0; l < n; ++l)
			for (size_t j = 0; j < n; ++j)
				if (dB[l][j] != 0)
					colsB[j].emplace_back((int)l, dB[l][j]);

		mpreal tmp;
		mpreal acc;
		auto addMul = [&](mpreal& a, const mpreal& b, double c)
			{
				mpfr_mul_d(tmp.mpfr_ptr(), b.mpfr_srcptr(), c, MPFR_RNDN);
				mpfr_add(a.mpfr_ptr(), a.mpfr_srcptr(), tmp.mpfr_srcptr(), MPFR_RNDN);
			};

		// y
		{
			std::vector<mpreal> ny(n);
			for (size_t j = 0; j < n; ++j)
			{
				acc = 0;
				for (const auto& e : colsB[j])
					addMul(acc, y[e.first], e.second);
				ny[j] = acc;
			}
			y.swap(ny);
		}
		// B
		{
			std::vector<mpreal> row(n);
			for (size_t k = 0; k < n; ++k)
			{
				for (size_t j = 0; j < n; ++j)
				{
					acc = 0;
					for (const auto& e : colsB[j])
						addMul(acc, B[k][e.first], e.second);
					row[j] = acc;
				}
				B[k].swap(row);
			}
		}
		// H
		{
			MMat nH(n, std::vector<mpreal>(n - 1, mpreal(0)));
			for (size_t i = 0; i < n; ++i)
				for (size_t l = 0; l < n; ++l)
				{
					const double c = dA[i][l];
					if (c == 0)
						continue;
					const size_t jEnd = std::min(l, n - 2);
					for (size_t j = 0; j <= jEnd; ++j)
						addMul(nH[i][j], H[l][j], c);
				}
			H.swap(nH);
		}
		lqDecompose(H);
	}

	std::atomic<long> g_m2lInvariantViolations(0);
}

long pslqMultipair2LInvariantViolations() { return g_m2lInvariantViolations.load(); }

PslqResult pslqMultipair2L(const std::vector<mpreal>& x, const mpreal& maxNormBound, int maxIter)
{
	const size_t n = x.size();
	if (n < 2)
		throw std::invalid_argument("pslqMultipair2L: need at least 2 numbers");

	const mp_prec_t prec = mpreal::get_default_prec();
	const double digits = mpfr::bits2digits(prec);
	const double reliable = 0.8 * digits;
	const mpreal detectEps = mpfr::pow(mpreal(10), -reliable);
	const mpreal maxB = mpfr::pow(mpreal(10), reliable); // proxy for the "max |A|" guard of pslq()
	const mpreal maxMeaningfulNorm = mpfr::pow(mpreal(10), reliable / n);
	const mpreal gamma = mpfr::sqrt(mpreal(4) / 3) * mpreal(1.0001);
	const double gammaD = std::sqrt(4.0 / 3.0) * 1.0001;
	const int maxPairs = std::max(1, (int)(kBeta * n));
	const bool verify = std::getenv("PSLQ_M2L_VERIFY") != nullptr;

	PslqResult result;

	// Normalize x, partial norms.
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
		throw std::invalid_argument("pslqMultipair2L: zero vector");
	const mpreal xnorm = s[0];
	std::vector<mpreal> y(n);
	for (size_t k = 0; k < n; ++k)
	{
		y[k] = x[k] / xnorm;
		s[k] /= xnorm;
	}
	for (size_t k = 0; k < n; ++k)
	{
		if (y[k] == 0)
		{
			result.status = PslqResult::Status::Found;
			result.relation.assign(n, mpreal(0));
			result.relation[k] = 1;
			result.residual = 0;
			return result;
		}
	}

	MMat H(n, std::vector<mpreal>(n - 1, mpreal(0)));
	for (size_t j = 0; j + 1 < n; ++j)
	{
		H[j][j] = s[j + 1] / s[j];
		for (size_t i = j + 1; i < n; ++i)
			H[i][j] = -y[i] * y[j] / (s[j] * s[j + 1]);
	}
	MMat B(n, std::vector<mpreal>(n, mpreal(0)));
	for (size_t i = 0; i < n; ++i)
		B[i][i] = 1;

	// Initial full reduction (multiprecision, A not needed).
	{
		// rows 1..n-1, columns up to n-2: emulate with a "pair" reaching everything.
		for (size_t i = 1; i < n; ++i)
		{
			for (size_t jj = std::min(i - 1, n - 2) + 1; jj-- > 0;)
			{
				const mpreal t = mpfr::round(H[i][jj] / H[jj][jj]);
				if (t == 0)
					continue;
				y[jj] += t * y[i];
				for (size_t k = 0; k <= jj; ++k)
					H[i][k] -= t * H[jj][k];
				for (size_t k = 0; k < n; ++k)
					B[k][jj] += t * B[k][i];
			}
		}
	}

	int it = 0;
	long nBatches = 0, nFail = 0, nMpIter = 0, nImq = 0, dpIter = 0;
	bool imq = false;
	bool forceMp = false;
	std::deque<std::vector<mpreal>> mpRing;

	// Double-precision work arrays (persistent to avoid reallocation).
	DMat dH(n, std::vector<double>(n - 1)), dA(n, std::vector<double>(n)), dB(n, std::vector<double>(n));
	DMat sH = dH, sA = dA, sB = dB;
	std::vector<double> dy(n), sy(n);

	const double targetBoundD = [&]
		{
			mpreal t = maxMeaningfulNorm;
			if (maxNormBound > 0 && maxNormBound < t)
				t = maxNormBound;
			return t > 1e300 ? 1e300 : t.toDouble();
		}();

	// Verification of the multiprecision invariants (debug aid, env PSLQ_M2L_VERIFY).
	auto verifyInvariants = [&]()
		{
			// y == B^T x / |x| and y^T H == 0.
			mpreal worst = 0;
			for (size_t j = 0; j < n; ++j)
			{
				mpreal v = 0;
				for (size_t k = 0; k < n; ++k)
					v += x[k] * B[k][j];
				v /= xnorm;
				worst = mpfr::max(worst, mpfr::abs(v - y[j]));
			}
			mpreal worstH = 0;
			for (size_t j = 0; j + 1 < n; ++j)
			{
				mpreal v = 0;
				for (size_t i = j; i < n; ++i)
					v += y[i] * H[i][j];
				worstH = mpfr::max(worstH, mpfr::abs(v));
			}
			// y accuracy is limited by cancellation; thresholds are the same "reliable" level as for detection.
			const bool bad = worst > detectEps || worstH > detectEps;
			if (bad)
				g_m2lInvariantViolations++;
			if (bad || std::getenv("PSLQ_M2L_VERIFY")[0] == '2')
				fprintf(stderr, "[m2l verify] |y-B^Tx/|x||=%s  max|y^T H_col|=%s %s\n", worst.toString(3).c_str(), worstH.toString(3).c_str(), bad ? "VIOLATION" : "ok");
		};

	// Full check on the (LQ-consistent) multiprecision state. Returns true when result is final.
	auto fullCheck = [&]() -> bool
		{
			result.iterations = it;
			mpreal maxH = 0;
			for (size_t j = 0; j + 1 < n; ++j)
				maxH = mpfr::max(maxH, mpfr::abs(H[j][j]));
			if (maxH > 0)
			{
				const mpreal bound = 1 / maxH;
				if (bound > result.normBound)
					result.normBound = bound;
			}
			size_t jMin = 0;
			for (size_t j = 1; j < n; ++j)
				if (mpfr::abs(y[j]) < mpfr::abs(y[jMin]))
					jMin = j;
			bool done = false;
			if (std::getenv("PSLQ_M2L_TRACE"))
				fprintf(stderr, "[m2l trace] it=%d bound=%s ymin=%s\n", it, (1 / maxH).toString(4).c_str(), mpfr::abs(y[jMin]).toString(3).c_str());
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
				const mpreal resid = mpfr::abs(r);
				// Independent guards: norm below the meaningful level, and the residual of the integer vector
				// itself (computed directly from x, not from the tracked y) is tiny.
				if (mpfr::sqrt(norm2) > maxMeaningfulNorm || resid > 10 * detectEps * xnorm)
					result.status = PslqResult::Status::PrecisionExhausted;
				else
				{
					result.status = PslqResult::Status::Found;
					// The batch may overshoot the exact moment of detection. For data that satisfy the relation only
					// approximately the bound can then exceed |relation| (bound is for exact relations of the given x);
					// clamp so that the reported bound never exceeds the norm of the relation found.
					const mpreal relNorm = mpfr::sqrt(norm2);
					if (result.normBound > relNorm)
						result.normBound = relNorm;
					result.relation = std::move(rel);
					result.residual = resid;
				}
				done = true;
			}
			else if (result.normBound >= maxMeaningfulNorm)
			{
				result.status = PslqResult::Status::PrecisionExhausted;
				done = true;
			}
			else if (maxNormBound > 0 && result.normBound >= maxNormBound)
			{
				result.status = PslqResult::Status::NormBoundReached;
				done = true;
			}
			else
			{
				mpreal maxAbsB = 0;
				for (size_t i = 0; i < n; ++i)
					for (size_t k = 0; k < n; ++k)
						maxAbsB = mpfr::max(maxAbsB, mpfr::abs(B[i][k]));
				if (maxAbsB > maxB)
				{
					result.status = PslqResult::Status::PrecisionExhausted;
					done = true;
				}
			}
			if (done && verify)
				verifyInvariants();
			if (done && std::getenv("PSLQ_M2L_TRACE"))
				fprintf(stderr, "[m2l stats] n=%zu iters=%d dpIters=%ld batches=%ld dpFail=%ld mpIters=%ld imq=%ld\n", n, it, dpIter, nBatches, nFail, nMpIter, nImq);
			return done;
		};

	if (fullCheck())
		return result;

	while (it < maxIter)
	{
		// dynamic range of y
		double range;
		{
			mpreal mn = mpfr::abs(y[0]), mx = mn;
			for (size_t j = 1; j < n; ++j)
			{
				const mpreal a = mpfr::abs(y[j]);
				if (a < mn) mn = a;
				if (a > mx) mx = a;
			}
			range = (mn / mx).toDouble();
		}

		bool doMp = forceMp || range < kDrep;
		forceMp = false;

		if (!doMp)
		{
			// ---- double precision batch ----
			double ymax = 0;
			for (size_t j = 0; j < n; ++j)
				ymax = std::max(ymax, mpfr::abs(y[j]).toDouble());
			// y are O(1) here but may be tiny in absolute terms: scale by the max.
			mpreal hT = 0;
			for (size_t j = 0; j + 1 < n; ++j)
				hT = mpfr::max(hT, mpfr::abs(H[j][j]));
			if (!(hT > 0))
			{
				forceMp = true;
				continue;
			}
			{
				const mpreal ymaxm = [&] { mpreal m = 0; for (auto& v : y) m = mpfr::max(m, mpfr::abs(v)); return m; }();
				for (size_t j = 0; j < n; ++j)
					dy[j] = (y[j] / ymaxm).toDouble();
			}
			for (size_t i = 0; i < n; ++i)
			{
				for (size_t j = 0; j + 1 < n; ++j)
					dH[i][j] = (j <= i) ? (H[i][j] / hT).toDouble() : 0.0;
				for (size_t j = 0; j < n; ++j)
				{
					dA[i][j] = (i == j) ? 1.0 : 0.0;
					dB[i][j] = (i == j) ? 1.0 : 0.0;
				}
			}
			const double hTd = hT.toDouble(); // dH = H / hT
			double maxAB = 1.0;
			int count = 0;
			bool failed = false;
			std::deque<std::vector<double>> ring;
			const int allowed = maxIter - it;
			while (count < allowed)
			{
				sH = dH; sA = dA; sB = dB; sy = dy;
				const double savedMaxAB = maxAB;
				const bool ok = multipairIteration<double>(dH, &dA, dB, dy, imq ? 1 : maxPairs, gammaD, &maxAB);
				imq = false;
				if (!ok)
				{
					dH.swap(sH); dA.swap(sA); dB.swap(sB); dy.swap(sy);
					maxAB = savedMaxAB;
					failed = true;
					++nFail;
					break;
				}
				++count;
				// cycle detection on dy
				for (const auto& old : ring)
					if (old == dy)
					{
						imq = true;
						break;
					}
				ring.push_back(dy);
				if ((int)ring.size() > kRing)
					ring.pop_front();
				// batch stop criteria
				double mn = std::fabs(dy[0]), mx = mn;
				for (size_t j = 1; j < n; ++j)
				{
					const double a = std::fabs(dy[j]);
					mn = std::min(mn, a);
					mx = std::max(mx, a);
				}
				if (!(mx > 0) || mn < kDepsRel * mx)
					break;
				if (maxAB > kSoftLimit)
					break;
				double hd = 0;
				for (size_t j = 0; j + 1 < n; ++j)
					hd = std::max(hd, std::fabs(dH[j][j]));
				if (hd > 0 && 1.0 / (hTd * hd) >= targetBoundD)
					break;
			}
			(void)failed;
			if (count == 0)
			{
				forceMp = true;
				continue;
			}
			it += count;
			++nBatches;
			dpIter += count;
			applyBatch(H, B, y, dA, dB);
		}
		else
		{
			// ---- multiprecision block ----
			for (int b = 0; b < kMpBlock && it < maxIter; ++b)
			{
				multipairIteration<mpreal>(H, nullptr, B, y, imq ? 1 : maxPairs, gamma, nullptr);
				imq = false;
				++it;
				++nMpIter;
				for (const auto& old : mpRing)
				{
					bool same = true;
					for (size_t j = 0; j < n && same; ++j)
						same = (old[j] == y[j]);
					if (same)
					{
						imq = true;
						++nImq;
						break;
					}
				}
				mpRing.push_back(y);
				if ((int)mpRing.size() > kRing)
					mpRing.pop_front();
				mpreal mn = mpfr::abs(y[0]);
				for (size_t j = 1; j < n; ++j)
					mn = mpfr::min(mn, mpfr::abs(y[j]));
				if (mn < detectEps)
					break;
			}
		}

		if (fullCheck())
			return result;
	}

	result.iterations = it;
	result.status = PslqResult::Status::MaxIterReached;
	return result;
}
