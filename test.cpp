#include <filesystem>
#include <string>
#include <numeric>
#include <atomic>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <chrono>
#include <functional>
#include <array>
#include <map>
#include <memory>
#include <set>
#include <cstdio>
#include <sstream>
#include "check_recurring_decimal.h"
#include "pslq.h"
#include "test.h"

using mutex_t = std::mutex;
mutex_t mutex;

void test_mpreal()
{
	//const int digits = 500;

	mpreal::set_default_prec(mpfr::digits2bits(digits));

	const mpreal pi = mpfr::const_pi();

	std::cout.precision(digits);

	std::cout << "pi = " << pi << std::endl;

	std::fstream fs;

	fs.open("pi.txt", std::ios::out);
	fs.precision(digits);
	fs << pi << std::endl;
	fs.close();
	std::cout << std::filesystem::current_path() << std::endl;

	mpreal apery;
	std::fstream fi;
	fi.open("apery.txt", std::ios::in);
	fi.precision(digits);
	fi >> apery;

	std::cout << "apery: " << apery << std::endl;
}

void test_multithreading_trigonometric_functions()
{
	const mpreal pi = mpfr::const_pi();

	std::cout.precision(digits);

	const mpreal pi_over_180 = pi / 180.;

	// conver degree to radian
	auto rad = [&pi_over_180](const auto& x_deg) -> mpreal
		{
			return x_deg * pi_over_180; // in rad
		};

	//std::cout << "sin (45 deg) = " << sin(rad(mpreal(45.))) << std::endl;

	std::thread t1([&rad](auto x)
		{
			std::unique_lock lock(mutex);
			std::cout << "sin (" << x << "deg) = " << sin(rad(mpreal(x))) << std::endl;
		}, mpreal(45.));

	std::thread t2([&rad](auto x)
		{
			std::unique_lock lock(mutex);
			std::cout << "cos (" << x << "deg) = " << sin(rad(mpreal(x))) << std::endl;
		}, mpreal(45.));

	t1.join();
	t2.join();
}

void test_check_recurring_decimal_performance()
{
	mpreal::set_default_prec(mpfr::digits2bits(digits));
	std::cout.precision(digits);

	// const mpreal x = (mpreal)(2.) / (mpreal)(3.);
	// const mpreal x = (mpreal)(7777768117.11) + (mpreal)(2.) / (mpreal)(71.) / (mpreal)(100.);
	// const mpreal x = mpfr::const_pi();
	const mpreal x = mpfr::sqrt(6);
	int p, n;
	auto t0 = std::chrono::steady_clock::now();
	bool r = checkRecurringDecimal(x, p, n);
	auto t1 = std::chrono::steady_clock::now();
	int p_fast, n_fast;
	bool r_fast = checkRecurringDecimalFast(x, p_fast, n_fast);
	auto t2 = std::chrono::steady_clock::now();
	std::cout << "x: " << x << std::endl;
	std::cout << "p: " << p << std::endl;
	std::cout << "n: " << n << std::endl;
	std::cout << "recurring decimal?: " << (r ? "true" : "false") << ", period: " << p << std::endl;
	std::cout << "fast: recurring decimal?: " << (r_fast ? "true" : "false") << ", period: " << p_fast << ", n: " << n_fast
		<< ((r == r_fast && p == p_fast && n == n_fast) ? " (same)" : " (DIFFERENT)") << std::endl;
	std::cout << "elapsed: original " << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms, fast "
		<< std::chrono::duration<double, std::milli>(t2 - t1).count() << " ms" << std::endl;
}

namespace
{
	std::string polynomialToString(const std::vector<mpreal>& c)
	{
		std::string out;
		for (size_t k = c.size(); k-- > 0;)
		{
			if (c[k] == 0)
				continue;
			const bool neg = c[k] < 0;
			const std::string a = mpfr::abs(c[k]).toString("%.0Rf");
			if (!out.empty())
				out += neg ? " - " : " + ";
			else if (neg)
				out += "-";
			if (k == 0 || a != "1")
				out += a;
			if (k >= 1)
				out += "x";
			if (k >= 2)
				out += "^" + std::to_string(k);
		}
		return out + " = 0";
	}

	void reportPolynomialSearch(const char* name, const mpreal& x, int maxDegree)
	{
		std::cout << "=== " << name << " ===" << std::endl;
		for (int d = 1; d <= maxDegree; ++d)
		{
			const auto t0 = std::chrono::steady_clock::now();
			const PslqResult r = findIntegerPolynomial(x, d);
			const auto t1 = std::chrono::steady_clock::now();
			std::cout << "degree " << d << ": " << toString(r.status) << ", iter " << r.iterations << ", "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms";
			if (r.status == PslqResult::Status::Found)
			{
				std::cout << "\n    " << polynomialToString(r.relation)
					<< "   (residual " << r.residual.toString("%.3Re") << ")" << std::endl;
				break;
			}
			std::cout << ", no polynomial with |coeffs|_2 < " << r.normBound.toString("%.3Re") << std::endl;
		}
	}
}

void test_pslq_apery_over_pi_cubed()
{
	const int pslqDigits = 1000;
	const int maxDegree = 20;
	mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
	std::cout.precision(50);

	// sanity checks with known algebraic numbers
	reportPolynomialSearch("sqrt(2) + sqrt(3)", mpfr::sqrt(mpreal(2)) + mpfr::sqrt(mpreal(3)), 8);
	reportPolynomialSearch("2^(1/5) + 1", mpfr::pow(mpreal(2), mpreal(1) / 5) + 1, 8);
	reportPolynomialSearch("zeta(2) / pi^2", mpfr::zeta(mpreal(2)) / mpfr::pow(mpfr::const_pi(), 2), 8);

	const mpreal zeta3 = mpfr::zeta(mpreal(3));
	const mpreal pi = mpfr::const_pi();
	const mpreal x = zeta3 / (pi * pi * pi);
	std::cout << "zeta(3)       = " << zeta3 << std::endl;
	std::cout << "zeta(3)/pi^3  = " << x << std::endl;
	reportPolynomialSearch("zeta(3) / pi^3", x, maxDegree);
}

namespace
{
	struct LogBasis
	{
		std::string name; // symbol printed for the base, e.g. "pi", "2", "(1+sqrt2)"
		mpreal log;       // log of the base; for "e" this is 1
	};

	std::string rationalToString(long long num, long long den)
	{
		if (den < 0)
		{
			num = -num;
			den = -den;
		}
		const long long g = std::gcd(num < 0 ? -num : num, den);
		num /= g;
		den /= g;
		return den == 1 ? std::to_string(num) : std::to_string(num) + "/" + std::to_string(den);
	}

	// Searches x = prod_i base_i^(rational_i) via PSLQ on (log x, log base_1, ..., log base_k).
	void reportMultiplicativeSearch(const char* name, const mpreal& x, const std::vector<LogBasis>& basis)
	{
		std::cout << "=== " << name << " ===" << std::endl;
		std::vector<mpreal> v;
		v.push_back(mpfr::log(x));
		for (const auto& b : basis)
			v.push_back(b.log);

		const auto t0 = std::chrono::steady_clock::now();
		const PslqResult r = pslq(v);
		const auto t1 = std::chrono::steady_clock::now();
		std::cout << toString(r.status) << ", iter " << r.iterations << ", "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms" << std::endl;

		if (r.status != PslqResult::Status::Found)
		{
			std::cout << "    no relation q*log(x) = sum c_i*log(base_i) with |(q, c)|_2 < "
				<< r.normBound.toString("%.3Re") << std::endl;
			return;
		}
		const long long q = r.relation[0].toLLong();
		if (q == 0)
		{
			std::cout << "    relation among basis constants only (x not involved):";
			for (size_t i = 0; i < basis.size(); ++i)
				std::cout << " " << r.relation[i + 1].toLLong() << "*log(" << basis[i].name << ")";
			std::cout << std::endl;
			return;
		}
		// q*log x + sum c_i log b_i = 0  =>  x = prod b_i^(-c_i/q)
		std::cout << "    x =";
		bool first = true;
		for (size_t i = 0; i < basis.size(); ++i)
		{
			const long long c = r.relation[i + 1].toLLong();
			if (c == 0)
				continue;
			std::cout << (first ? " " : " * ") << basis[i].name << "^(" << rationalToString(-c, q) << ")";
			first = false;
		}
		if (first)
			std::cout << " 1";
		std::cout << "   (residual " << r.residual.toString("%.3Re") << ")" << std::endl;
	}
}

// Tests whether zeta(3)/pi^3 = alpha * exp(r) * pi^s with r, s rational and alpha in the
// multiplicative group generated by rational powers of primes and of a few quadratic units.
// (A general algebraic alpha cannot be handled this way: log(alpha) has no finite basis.)
void test_pslq_apery_multiplicative()
{
	const int pslqDigits = 1000;
	mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
	std::cout.precision(50);

	const mpreal pi = mpfr::const_pi();
	const mpreal sqrt2 = mpfr::sqrt(mpreal(2)), sqrt3 = mpfr::sqrt(mpreal(3)), sqrt5 = mpfr::sqrt(mpreal(5));

	std::vector<LogBasis> basis;
	basis.push_back({ "e", mpreal(1) });
	basis.push_back({ "pi", mpfr::log(pi) });
	for (int p : { 2, 3, 5, 7, 11, 13, 17, 19, 23 })
		basis.push_back({ std::to_string(p), mpfr::log(mpreal(p)) });
	basis.push_back({ "((1+sqrt5)/2)", mpfr::log((1 + sqrt5) / 2) });
	basis.push_back({ "(1+sqrt2)", mpfr::log(1 + sqrt2) });
	basis.push_back({ "(2+sqrt3)", mpfr::log(2 + sqrt3) });

	// sanity check: 3^(1/3) * e^(2/5) * pi^(-7/4) * (1+sqrt2)^(1/2) / 7^(2/3)
	const mpreal known = mpfr::cbrt(mpreal(3)) * mpfr::exp(mpreal(2) / 5) * mpfr::pow(pi, mpreal(-7) / 4)
		* mpfr::sqrt(1 + sqrt2) / mpfr::pow(mpreal(7), mpreal(2) / 3);
	reportMultiplicativeSearch("3^(1/3) * e^(2/5) * pi^(-7/4) * (1+sqrt2)^(1/2) * 7^(-2/3)", known, basis);

	const mpreal x = mpfr::zeta(mpreal(3)) / (pi * pi * pi);
	reportMultiplicativeSearch("zeta(3) / pi^3", x, basis);
}

namespace
{
	struct Rational
	{
		long long num, den;
		std::string str() const { return rationalToString(num, den); }
	};

	// All reduced fractions p/q with 1 <= q <= maxDen and |p/q| <= maxAbs, sorted by value.
	std::vector<Rational> enumerateRationals(int maxDen, int maxAbs)
	{
		std::vector<Rational> out;
		for (long long q = 1; q <= maxDen; ++q)
			for (long long p = -maxAbs * q; p <= maxAbs * q; ++p)
				if (std::gcd(p < 0 ? -p : p, q) == 1)
					out.push_back({ p, q });
		std::sort(out.begin(), out.end(), [](const Rational& a, const Rational& b) { return a.num * b.den < b.num * a.den; });
		return out;
	}

	// Searches x = alpha * e^r * pi^s with alpha algebraic of degree <= degree and
	// r, s rational (denominator <= maxDen, |r|, |s| <= maxAbs), by running PSLQ on
	// (1, y, ..., y^degree) for y = x * e^(-r) * pi^(-s) over the whole (r, s) grid.
	// xAtPrec(digits) must return x computed at the given number of digits (used to verify hits).
	template <class XAtPrec>
	void reportAlgebraicExpPiSearch(const char* name, XAtPrec xAtPrec, int pslqDigits, int degree, int maxDen, int maxAbs, bool piGrid = true)
	{
		std::cout << "=== " << name << " ===" << std::endl;
		mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
		const mpreal x = xAtPrec(pslqDigits);
		const mpreal pi = mpfr::const_pi();
		const std::vector<Rational> rats = enumerateRationals(maxDen, maxAbs);
		// without the pi grid, s is fixed to 0
		const std::vector<Rational> ratsS = piGrid ? rats : std::vector<Rational>{ { 0, 1 } };
		const size_t nr = rats.size(), ns = ratsS.size();
		std::vector<mpreal> expMinus(nr), piMinus(ns);
		for (size_t i = 0; i < nr; ++i)
			expMinus[i] = mpfr::exp(-mpreal(rats[i].num) / rats[i].den);
		for (size_t i = 0; i < ns; ++i)
			piMinus[i] = mpfr::pow(pi, -mpreal(ratsS[i].num) / ratsS[i].den);
		const size_t total = nr * ns;
		std::cout << "degree <= " << degree << ", " << pslqDigits << " digits, denominators <= " << maxDen
			<< (piGrid ? ", |r|,|s| <= " : ", |r| <= ") << maxAbs << ": " << total << (piGrid ? " (r, s) pairs" : " values of r (s = 0)") << std::endl;

		std::vector<PslqResult> results(total);
		std::atomic<size_t> next{ 0 }, done{ 0 };
		const unsigned nThreads = std::max(1u, std::thread::hardware_concurrency());
		const auto t0 = std::chrono::steady_clock::now();
		std::vector<std::thread> threads;
		for (unsigned t = 0; t < nThreads; ++t)
		{
			threads.emplace_back([&]()
				{
					mpreal::set_default_prec(mpfr::digits2bits(pslqDigits)); // default precision is thread-local
					for (size_t k; (k = next++) < total;)
					{
						const mpreal y = x * expMinus[k / ns] * piMinus[k % ns];
						results[k] = findIntegerPolynomial(y, degree);
						const size_t d = ++done;
						if (d % 1000 == 0)
						{
							std::unique_lock lock(mutex);
							std::cout << "  " << d << " / " << total << std::endl;
						}
					}
				});
		}
		for (auto& th : threads)
			th.join();
		const auto t1 = std::chrono::steady_clock::now();
		std::cout << "elapsed " << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count() << " s" << std::endl;

		// Verify hits at doubled precision; collect the weakest norm bound among the non-hits.
		mpreal minBound = -1;
		size_t minBoundIdx = 0;
		int hits = 0;
		for (size_t k = 0; k < total; ++k)
		{
			PslqResult& r = results[k];
			const Rational& er = rats[k / ns];
			const Rational& ps = ratsS[k % ns];
			if (r.status != PslqResult::Status::Found)
			{
				if (minBound < 0 || r.normBound < minBound)
				{
					minBound = r.normBound;
					minBoundIdx = k;
				}
				continue;
			}
			// PSLQ may return a multiple of the minimal polynomial; retry with increasing degree.
			{
				const mpreal y = x * expMinus[k / ns] * piMinus[k % ns];
				for (int d = 1; d < degree; ++d)
				{
					PslqResult rd = findIntegerPolynomial(y, d);
					if (rd.status == PslqResult::Status::Found)
					{
						r = std::move(rd);
						break;
					}
				}
			}
			mpreal::set_default_prec(mpfr::digits2bits(2 * pslqDigits));
			const mpreal er2 = mpreal(er.num) / er.den, ps2 = mpreal(ps.num) / ps.den;
			const mpreal y2 = xAtPrec(2 * pslqDigits) * mpfr::exp(-er2) * mpfr::pow(mpfr::const_pi(), -ps2);
			mpreal v = 0;
			for (size_t i = r.relation.size(); i-- > 0;)
				v = v * y2 + r.relation[i];
			mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
			const bool genuine = mpfr::abs(v) < mpfr::pow(mpreal(10), -1.5 * pslqDigits);
			hits += genuine;
			std::cout << (genuine ? "  HIT" : "  spurious (failed verification)") << ": r = " << er.str() << ", s = " << ps.str()
				<< ", alpha = x*e^(-r)*pi^(-s) is a root of " << polynomialToString(r.relation)
				<< "   (residual at " << 2 * pslqDigits << " digits: " << mpfr::abs(v).toString("%.3Re") << ")" << std::endl;
		}
		std::cout << hits << " verified hit(s)." << std::endl;
		if (minBound >= 0)
			std::cout << "For every other (r, s): no integer polynomial of degree <= " << degree
				<< " with |coeffs|_2 < " << minBound.toString("%.3Re") << " (weakest at r = " << rats[minBoundIdx / ns].str()
				<< ", s = " << ratsS[minBoundIdx % ns].str() << ")" << std::endl;
	}
}

void test_pslq_apery_algebraic_exp_pi_grid()
{
	std::cout.precision(50);

	// sanity check: (2^(1/3) + 1) * e^(1/2) * pi^(-2/3)
	reportAlgebraicExpPiSearch("(2^(1/3) + 1) * e^(1/2) * pi^(-2/3)  [small grid]", [](int d)
		{
			mpreal::set_default_prec(mpfr::digits2bits(d));
			return (mpfr::cbrt(mpreal(2)) + 1) * mpfr::exp(mpreal(1) / 2) * mpfr::pow(mpfr::const_pi(), mpreal(-2) / 3);
		}, 500, 8, 3, 1);

	reportAlgebraicExpPiSearch("zeta(3) / pi^3", [](int d)
		{
			mpreal::set_default_prec(mpfr::digits2bits(d));
			const mpreal pi = mpfr::const_pi();
			return mpfr::zeta(mpreal(3)) / (pi * pi * pi);
		}, 500, 8, 6, 6);
}

namespace
{
	// Runs findIntegerPolynomial for every degree 1..maxDegree on each target, in parallel over (target, degree).
	// xAtPrec(digits) must return x at the given number of digits (Found results are re-checked at 2x digits).
	struct PolySearchTarget
	{
		std::string name;
		std::function<mpreal(int)> xAtPrec;
	};

	struct PolySearchSummary
		{
			std::vector<std::vector<PslqResult>> results; // [target][degree]
			std::vector<int> verified;
		};

		PolySearchSummary reportPolynomialSearchParallel(const std::vector<PolySearchTarget>& targets, int pslqDigits, int maxDegree,
			bool onlyMaxDegree = false, unsigned maxThreads = 0)
	{
		mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
		// onlyMaxDegree: run only degree == maxDegree (the vector (1,x,..,x^d) contains every lower-degree polynomial).
		const int minDeg = onlyMaxDegree ? maxDegree : 1;
		const size_t nt = targets.size();
		std::vector<mpreal> xs(nt);
		for (size_t t = 0; t < nt; ++t)
			xs[t] = targets[t].xAtPrec(pslqDigits);

		// jobs: (target, degree), biggest degree first for load balancing
		std::vector<std::pair<size_t, int>> jobs;
		for (int d = maxDegree; d >= minDeg; --d)
			for (size_t t = 0; t < nt; ++t)
				jobs.push_back({ t, d });
		std::vector<std::vector<PslqResult>> results(nt, std::vector<PslqResult>(maxDegree + 1));
		std::vector<std::vector<double>> seconds(nt, std::vector<double>(maxDegree + 1));

		std::atomic<size_t> next{ 0 }, done{ 0 };
		const unsigned nThreads = maxThreads ? maxThreads : std::max(1u, std::thread::hardware_concurrency());
		const auto t0 = std::chrono::steady_clock::now();
		std::vector<std::thread> threads;
		for (unsigned th = 0; th < nThreads; ++th)
		{
			threads.emplace_back([&]()
				{
					mpreal::set_default_prec(mpfr::digits2bits(pslqDigits)); // default precision is thread-local
					for (size_t k; (k = next++) < jobs.size();)
					{
						const auto [t, d] = jobs[k];
						const auto s0 = std::chrono::steady_clock::now();
						results[t][d] = findIntegerPolynomial(xs[t], d);
						seconds[t][d] = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
						const size_t c = ++done;
						std::unique_lock lock(mutex);
						std::cout << "  [" << c << " / " << jobs.size() << "] " << targets[t].name << ", degree " << d << ": "
							<< toString(results[t][d].status) << " (" << static_cast<long>(seconds[t][d] * 1000) << " ms, normBound "
								<< results[t][d].normBound.toString("%.2Re") << ", elapsed "
								<< static_cast<long>(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()) << " s)" << std::endl;
					}
				});
		}
		for (auto& th : threads)
			th.join();
		std::cout << "elapsed " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count()
			<< " s" << std::endl;

		std::vector<int> verifiedCounts(nt, 0);
		for (size_t t = 0; t < nt; ++t)
		{
			std::cout << "=== " << targets[t].name << " (" << pslqDigits << " digits) ===" << std::endl;
			int verified = 0;
			for (int d = minDeg; d <= maxDegree; ++d)
			{
				const PslqResult& r = results[t][d];
				std::cout << "degree " << d << ": " << toString(r.status) << ", iter " << r.iterations << ", "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration<double>(seconds[t][d])).count() << " ms";
				if (r.status != PslqResult::Status::Found)
				{
					std::cout << ", no polynomial with |coeffs|_2 < " << r.normBound.toString("%.3Re") << std::endl;
					continue;
				}
				mpreal::set_default_prec(mpfr::digits2bits(2 * pslqDigits));
				const mpreal x2 = targets[t].xAtPrec(2 * pslqDigits);
				mpreal v = 0;
				for (size_t i = r.relation.size(); i-- > 0;)
					v = v * x2 + r.relation[i];
				const bool genuine = mpfr::abs(v) < mpfr::pow(mpreal(10), -1.5 * pslqDigits);
				mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
				verified += genuine;
				int trueDeg = 0;
				for (size_t i = 0; i < r.relation.size(); ++i)
					if (r.relation[i] != 0) trueDeg = static_cast<int>(i);
				std::cout << (genuine ? "  VERIFIED" : "  SPURIOUS") << " (true degree " << trueDeg << ")\n    " << polynomialToString(r.relation)
					<< "   (residual at " << 2 * pslqDigits << " digits: " << mpfr::abs(v).toString("%.3Re") << ")" << std::endl;
			}
			std::cout << verified << " verified relation(s) for " << targets[t].name << std::endl;
			verifiedCounts[t] = verified;
		}
		return { results, verifiedCounts };
	}
}

void test_pslq_zeta3_zeta5_ratio_parallel(int maxDegree)
{
	const int pslqDigits = 500;
	std::cout.precision(60);
	auto zetaRatio = [](int p, int q)
		{
			return [p, q](int digits)
				{
					mpreal::set_default_prec(mpfr::digits2bits(digits));
					return mpfr::pow(mpfr::zeta(mpreal(3)), mpreal(1) / p) / mpfr::pow(mpfr::zeta(mpreal(5)), mpreal(1) / q);
				};
		};
	reportPolynomialSearchParallel({
		{ "zeta(3)^(1/3) / zeta(5)^(1/5)", zetaRatio(3, 5) },
		{ "zeta(3)^(1/2) / zeta(5)^(1/4)", zetaRatio(2, 4) },
		}, pslqDigits, maxDegree);
}

void test_pslq_zeta3_zeta7_over_zeta5_squared(int pslqDigits, int maxDegree)
{
	std::cout.precision(60);
	reportPolynomialSearchParallel({
		{ "zeta(3) * zeta(7) / zeta(5)^2", [](int digits)
			{
				mpreal::set_default_prec(mpfr::digits2bits(digits));
				const mpreal z5 = mpfr::zeta(mpreal(5));
				return mpfr::zeta(mpreal(3)) * mpfr::zeta(mpreal(7)) / (z5 * z5);
			} },
		}, pslqDigits, maxDegree);
}

void test_pslq_zeta3_zeta7_over_zeta5_squared_exp_grid(int pslqDigits, int degree, int maxDen, int maxAbs)
{
	std::cout.precision(50);

	// sanity check: (2^(1/3) + 1) * e^(5/7)
	reportAlgebraicExpPiSearch("(2^(1/3) + 1) * e^(5/7)  [sanity]", [](int d)
		{
			mpreal::set_default_prec(mpfr::digits2bits(d));
			return (mpfr::cbrt(mpreal(2)) + 1) * mpfr::exp(mpreal(5) / 7);
		}, pslqDigits, degree, maxDen, maxAbs, false);

	reportAlgebraicExpPiSearch("zeta(3) * zeta(7) / zeta(5)^2", [](int d)
		{
			mpreal::set_default_prec(mpfr::digits2bits(d));
			const mpreal z5 = mpfr::zeta(mpreal(5));
			return mpfr::zeta(mpreal(3)) * mpfr::zeta(mpreal(7)) / (z5 * z5);
		}, pslqDigits, degree, maxDen, maxAbs, false);
}

// Tests whether zeta(3)*zeta(7)/zeta(5)^2 = alpha * e^r with r rational (any denominator) and alpha in the
// multiplicative group generated by rational powers of primes and of a few quadratic units,
// via PSLQ on (log x, 1, log 2, log 3, ..., log units).
void test_pslq_zeta3_zeta7_over_zeta5_squared_multiplicative(int pslqDigits, int maxPrime)
{
	mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
	std::cout.precision(50);

	const mpreal sqrt2 = mpfr::sqrt(mpreal(2)), sqrt3 = mpfr::sqrt(mpreal(3)), sqrt5 = mpfr::sqrt(mpreal(5));
	std::vector<LogBasis> basis;
	basis.push_back({ "e", mpreal(1) });
	for (int p = 2; p <= maxPrime; ++p)
	{
		bool isPrime = true;
		for (int d = 2; d * d <= p && isPrime; ++d)
			isPrime = p % d != 0;
		if (isPrime)
			basis.push_back({ std::to_string(p), mpfr::log(mpreal(p)) });
	}
	basis.push_back({ "((1+sqrt5)/2)", mpfr::log((1 + sqrt5) / 2) });
	basis.push_back({ "(1+sqrt2)", mpfr::log(1 + sqrt2) });
	basis.push_back({ "(2+sqrt3)", mpfr::log(2 + sqrt3) });
	std::cout << "primes <= " << maxPrime << ", " << pslqDigits << " digits, PSLQ vector length " << basis.size() + 1 << std::endl;

	// sanity check with a large-denominator r: 2^(1/3) * e^(357/1000) * (1+sqrt2)^(-1/2) * 11^(3/7)
	const mpreal known = mpfr::cbrt(mpreal(2)) * mpfr::exp(mpreal(357) / 1000) / mpfr::sqrt(1 + sqrt2)
		* mpfr::pow(mpreal(11), mpreal(3) / 7);
	reportMultiplicativeSearch("2^(1/3) * e^(357/1000) * (1+sqrt2)^(-1/2) * 11^(3/7)  [sanity]", known, basis);

	const mpreal z5 = mpfr::zeta(mpreal(5));
	const mpreal x = mpfr::zeta(mpreal(3)) * mpfr::zeta(mpreal(7)) / (z5 * z5);
	reportMultiplicativeSearch("zeta(3) * zeta(7) / zeta(5)^2", x, basis);
}

namespace
{
	// PSLQ on (log b_1, ..., log b_n): searches prod b_i^(c_i) = 1 with integers c_i.
	void reportLogRelationSearch(const char* name, const std::vector<LogBasis>& basis)
	{
		std::cout << "=== " << name << " ===" << std::endl;
		std::vector<mpreal> v;
		for (const auto& b : basis)
			v.push_back(b.log);
		const auto t0 = std::chrono::steady_clock::now();
		const PslqResult r = pslq(v);
		const auto t1 = std::chrono::steady_clock::now();
		std::cout << toString(r.status) << ", iter " << r.iterations << ", "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms" << std::endl;
		if (r.status != PslqResult::Status::Found)
		{
			std::cout << "    no relation prod base_i^(c_i) = 1 with integer c, |c|_2 < " << r.normBound.toString("%.3Re") << std::endl;
			return;
		}
		std::cout << "    ";
		bool first = true;
		for (size_t i = 0; i < basis.size(); ++i)
		{
			const std::string c = r.relation[i].toString("%.0Rf");
			if (c == "0" || c == "-0")
				continue;
			std::cout << (first ? "" : " * ") << basis[i].name << "^(" << c << ")";
			first = false;
		}
		std::cout << " = 1   (residual " << r.residual.toString("%.3Re") << ")" << std::endl;
	}
}

// PSLQ on (log zeta(3), log zeta(5), log zeta(7), 1, log 2, log 3, ..., log p, log units):
// searches zeta(3)^a * zeta(5)^b * zeta(7)^c = e^r * prod p^(a_p) * units with rational exponents.
void test_pslq_zeta357_multiplicative(int pslqDigits, int maxPrime)
{
	mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
	std::cout.precision(50);

	const mpreal sqrt2 = mpfr::sqrt(mpreal(2)), sqrt3 = mpfr::sqrt(mpreal(3)), sqrt5 = mpfr::sqrt(mpreal(5));
	const mpreal z3 = mpfr::zeta(mpreal(3)), z5 = mpfr::zeta(mpreal(5)), z7 = mpfr::zeta(mpreal(7));
	std::vector<LogBasis> rest;
	rest.push_back({ "e", mpreal(1) });
	for (int p = 2; p <= maxPrime; ++p)
	{
		bool isPrime = true;
		for (int d = 2; d * d <= p && isPrime; ++d)
			isPrime = p % d != 0;
		if (isPrime)
			rest.push_back({ std::to_string(p), mpfr::log(mpreal(p)) });
	}
	rest.push_back({ "((1+sqrt5)/2)", mpfr::log((1 + sqrt5) / 2) });
	rest.push_back({ "(1+sqrt2)", mpfr::log(1 + sqrt2) });
	rest.push_back({ "(2+sqrt3)", mpfr::log(2 + sqrt3) });
	std::cout << "primes <= " << maxPrime << ", " << pslqDigits << " digits, PSLQ vector length " << rest.size() + 3 << std::endl;

	// sanity check: replace zeta(7) by y = zeta(3)^2 * zeta(5)^(-1/3) * 2^(1/3) * e^(357/1000)
	{
		const mpreal y = z3 * z3 / mpfr::cbrt(z5) * mpfr::cbrt(mpreal(2)) * mpfr::exp(mpreal(357) / 1000);
		std::vector<LogBasis> basis{ { "zeta(3)", mpfr::log(z3) }, { "zeta(5)", mpfr::log(z5) }, { "y", mpfr::log(y) } };
		basis.insert(basis.end(), rest.begin(), rest.end());
		reportLogRelationSearch("y = zeta(3)^2 * zeta(5)^(-1/3) * 2^(1/3) * e^(357/1000)  [sanity]", basis);
	}

	std::vector<LogBasis> basis{ { "zeta(3)", mpfr::log(z3) }, { "zeta(5)", mpfr::log(z5) }, { "zeta(7)", mpfr::log(z7) } };
	basis.insert(basis.end(), rest.begin(), rest.end());
	reportLogRelationSearch("zeta(3), zeta(5), zeta(7)", basis);
}

// PSLQ on (log zeta(s_1), ..., log zeta(s_k), [log pi], 1, log 2, ..., log p, log units):
// searches prod zeta(s_i)^(a_i) * pi^b = e^r * prod p^(a_p) * units with rational exponents.
void test_pslq_zeta_multiplicative(int pslqDigits, int maxPrime, const std::vector<int>& zetaArgs, bool includePi, bool includeEuler)
{
	mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
	std::cout.precision(50);

	const mpreal sqrt2 = mpfr::sqrt(mpreal(2)), sqrt3 = mpfr::sqrt(mpreal(3)), sqrt5 = mpfr::sqrt(mpreal(5));
	std::vector<LogBasis> zetas;
	for (int s : zetaArgs)
		zetas.push_back({ "zeta(" + std::to_string(s) + ")", mpfr::log(mpfr::zeta(mpreal(s))) });
	std::vector<LogBasis> rest;
	if (includePi)
		rest.push_back({ "pi", mpfr::log(mpfr::const_pi()) });
	if (includeEuler)
		rest.push_back({ "gamma", mpfr::log(mpfr::const_euler()) });
	rest.push_back({ "e", mpreal(1) });
	for (int p = 2; p <= maxPrime; ++p)
	{
		bool isPrime = true;
		for (int d = 2; d * d <= p && isPrime; ++d)
			isPrime = p % d != 0;
		if (isPrime)
			rest.push_back({ std::to_string(p), mpfr::log(mpreal(p)) });
	}
	rest.push_back({ "((1+sqrt5)/2)", mpfr::log((1 + sqrt5) / 2) });
	rest.push_back({ "(1+sqrt2)", mpfr::log(1 + sqrt2) });
	rest.push_back({ "(2+sqrt3)", mpfr::log(2 + sqrt3) });
	std::cout << "primes <= " << maxPrime << ", " << pslqDigits << " digits, PSLQ vector length "
		<< zetas.size() + rest.size() << std::endl;

	// sanity check: replace the last zeta by y = zeta(s_1)^2 * zeta(s_2)^(-1/3) * pi^(5/4) * 2^(1/3) * e^(357/1000)
	if (zetas.size() >= 3)
	{
		const mpreal z1 = mpfr::zeta(mpreal(zetaArgs[0])), z2 = mpfr::zeta(mpreal(zetaArgs[1]));
		mpreal y = z1 * z1 / mpfr::cbrt(z2) * mpfr::cbrt(mpreal(2)) * mpfr::exp(mpreal(357) / 1000);
		if (includePi)
			y *= mpfr::pow(mpfr::const_pi(), mpreal(5) / 4);
		if (includeEuler)
			y *= mpfr::pow(mpfr::const_euler(), mpreal(-2) / 7);
		std::vector<LogBasis> basis(zetas.begin(), zetas.end() - 1);
		basis.push_back({ "y", mpfr::log(y) });
		basis.insert(basis.end(), rest.begin(), rest.end());
		const std::string name = "y = " + zetas[0].name + "^2 * " + zetas[1].name + "^(-1/3)"
			+ (includePi ? " * pi^(5/4)" : "") + (includeEuler ? " * gamma^(-2/7)" : "") + " * 2^(1/3) * e^(357/1000)  [sanity]";
		reportLogRelationSearch(name.c_str(), basis);
	}

	std::vector<LogBasis> basis = zetas;
	basis.insert(basis.end(), rest.begin(), rest.end());
	std::string name;
	for (const auto& z : zetas)
		name += (name.empty() ? "" : ", ") + z.name;
	if (includePi)
		name += ", pi";
	if (includeEuler)
		name += ", gamma";
	reportLogRelationSearch(name.c_str(), basis);
}

void test_pslq_zeta3_zeta5_ratio(int p, int q)
{
	// x = zeta(3)^(1/p) / zeta(5)^(1/q): search for an integer polynomial P with P(x) = 0
	// via PSLQ on (1, x, ..., x^d) at 500 digits.
	const int pslqDigits = 500;
	const int maxDegree = 24;
	mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
	std::cout.precision(60);

	const std::string expo = "^(1/" + std::to_string(p) + ") / ";
	const std::string expq = "^(1/" + std::to_string(q) + ")";

	// sanity check with a known algebraic number of the same shape
	const std::string sanityName = "2" + expo + "3" + expq;
	reportPolynomialSearch(sanityName.c_str(),
		mpfr::pow(mpreal(2), mpreal(1) / p) / mpfr::pow(mpreal(3), mpreal(1) / q), maxDegree);

	const mpreal zeta3 = mpfr::zeta(mpreal(3));
	const mpreal zeta5 = mpfr::zeta(mpreal(5));
	const mpreal x = mpfr::pow(zeta3, mpreal(1) / p) / mpfr::pow(zeta5, mpreal(1) / q);
	const std::string name = "zeta(3)" + expo + "zeta(5)" + expq;
	std::cout << name << " = " << x << std::endl;
	reportPolynomialSearch(name.c_str(), x, maxDegree);
}

// Li_3(z) = sum z^k / k^3 for real |z| <= 1/2, at the current default precision.
mpreal polylog3(const mpreal& z)
{
	const mpreal eps = mpfr::pow(mpreal(2), -(long)mpreal::get_default_prec() - 40);
	mpreal sum = 0, zk = 1;
	for (long k = 1;; ++k)
	{
		zk *= z;
		const mpreal term = zk / mpreal(k * k * k);
		sum += term;
		if (mpfr::abs(term) < eps) break;
	}
	return sum;
}

// PSLQ on the vector computeVector(digits); prints the relation with labels and verifies it at 2*digits.
struct LinearSearchOutcome { bool found = false; bool verified = false; std::string normBound; std::vector<std::string> coeffs; };
LinearSearchOutcome reportLinearRelationSearch(std::ostream& out, const std::string& name, const std::function<std::vector<mpreal>(int digits)>& computeVector,
	const std::vector<std::string>& labels, int digits)
{
	LinearSearchOutcome outcome;
	const auto t0 = std::chrono::steady_clock::now();
	mpreal::set_default_prec(mpfr::digits2bits(digits));
	const std::vector<mpreal> c = computeVector(digits);
	const int n = (int)c.size();
	out << "=== " << name << " (n = " << n << ", digits = " << digits << ")" << std::endl;
	PslqResult r = pslq(c);
	outcome.normBound = r.normBound.toString("%.2Re");
	out << "status = " << toString(r.status) << ", iterations = " << r.iterations
		<< ", normBound = " << r.normBound.toString("%.3Re") << std::endl;
	if (r.status == PslqResult::Status::Found)
	{
		std::string s;
		for (int i = 0; i < n; ++i)
		{
			if (r.relation[i] == 0) continue;
			const std::string v = r.relation[i].toString("%.0Rf");
			s += (v[0] == '-' ? " " : " +") + v + "*" + labels[i];
		}
		out << "relation:" << s << std::endl;
		for (int i = 0; i < n; ++i) outcome.coeffs.push_back(r.relation[i].toString("%.0Rf"));
		mpreal::set_default_prec(mpfr::digits2bits(2 * digits));
		const std::vector<mpreal> c2 = computeVector(2 * digits);
		mpreal v = 0;
		for (int i = 0; i < n; ++i) v += r.relation[i] * c2[i];
		const bool ok = mpfr::abs(v) < mpfr::pow(mpreal(10), -1.6 * digits);
		out << "residual at " << 2 * digits << " digits = " << mpfr::abs(v).toString("%.3Re")
			<< (ok ? "  VERIFIED" : "  FAILED verification") << std::endl;
		outcome.found = true; outcome.verified = ok;
		mpreal::set_default_prec(mpfr::digits2bits(digits));
	}
	out << "power estimate 10^(0.8*D/n) = 1e" << (int)(0.8 * digits / n) << std::endl;
	out << "elapsed " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s" << std::endl;
	return outcome;
}

void reportLinearRelationSearch(const std::string& name, const std::function<std::vector<mpreal>(int digits)>& computeVector,
	const std::vector<std::string>& labels, int digits)
{
	reportLinearRelationSearch(std::cout, name, computeVector, labels, digits);
}

void test_pslq_zeta3_additive_elementary(int digits, bool includeLog5)
{
	std::cout.precision(30);
	// letters
	std::vector<std::string> let = { "pi", "log2", "log3" };
	if (includeLog5) let.push_back("log5");
	const int nl = (int)let.size();
	auto letters = [nl](int) {
		std::vector<mpreal> v = { mpfr::const_pi(), mpfr::log(mpreal(2)), mpfr::log(mpreal(3)) };
		if (nl == 4) v.push_back(mpfr::log(mpreal(5)));
		return v;
	};
	std::vector<std::string> baseLabels;
	std::vector<std::function<mpreal(const std::vector<mpreal>&, const mpreal&)>> unused;
	std::vector<std::array<int, 3>> mons;
	for (int a = 0; a < nl; ++a)
		for (int b = a; b < nl; ++b)
			for (int c = b; c < nl; ++c)
			{
				mons.push_back({ a, b, c });
				baseLabels.push_back(let[a] + "*" + let[b] + "*" + let[c]);
			}
	for (int a = 0; a < nl; ++a) baseLabels.push_back("G*" + let[a]);
	auto basis = [=](int d) {
		const std::vector<mpreal> L = letters(d);
		std::vector<mpreal> v;
		for (auto& m : mons) v.push_back(L[m[0]] * L[m[1]] * L[m[2]]);
		const mpreal G = mpfr::const_catalan();
		for (int a = 0; a < nl; ++a) v.push_back(G * L[a]);
		return v;
	};
	auto withPrefix = [&](std::vector<std::string> pre, std::function<std::vector<mpreal>(int)> pf) {
		std::vector<std::string> lab = pre;
		lab.insert(lab.end(), baseLabels.begin(), baseLabels.end());
		return std::make_pair(lab, [=](int d) {
			std::vector<mpreal> v = pf(d);
			const std::vector<mpreal> b = basis(d);
			v.insert(v.end(), b.begin(), b.end());
			return v;
		});
	};
	auto z3 = [](int) { return mpfr::zeta(mpreal(3)); };

	// V1
	reportLinearRelationSearch("V1 [sanity] Li3(1/2), zeta3, pi^2 log2, log2^3  expect (24,-21,2,-4)",
		[&](int d) {
			const mpreal pi = mpfr::const_pi(), l2 = mpfr::log(mpreal(2));
			return std::vector<mpreal>{ polylog3(mpreal(1) / 2), z3(d), pi * pi * l2, l2 * l2 * l2 };
		},
		{ "Li3(1/2)", "zeta(3)", "pi^2*log2", "log2^3" }, digits);
	// V2
	{
		auto p = withPrefix({ "Li3(1/2)", "zeta(3)" }, [&](int d) { return std::vector<mpreal>{ polylog3(mpreal(1) / 2), z3(d) }; });
		reportLinearRelationSearch("V2 [sanity] Li3(1/2), zeta3, full basis  expect same as V1", p.second, p.first, digits);
	}
	// V3
	{
		auto p = withPrefix({ "Li3(1/3)", "Li3(-1/3)", "zeta(3)" },
			[&](int d) { return std::vector<mpreal>{ polylog3(mpreal(1) / 3), polylog3(mpreal(-1) / 3), z3(d) }; });
		reportLinearRelationSearch("V3 [sanity] Li3(1/3), Li3(-1/3), zeta3, full basis  expect (-12, 6, 13) & pi^2 log3 -1, log3^3 +1", p.second, p.first, digits);
	}
	// M
	{
		auto p = withPrefix({ "zeta(3)" }, [&](int d) { return std::vector<mpreal>{ z3(d) }; });
		reportLinearRelationSearch("M zeta3, full basis", p.second, p.first, digits);
	}
}

// Homogeneous weight-w monomials pi^(2a) * prod zeta(k)^(b_k) (k odd, 3 <= k <= maxZetaArg), one PSLQ job per weight,
// each paired with a check job that appends a known-dependent element.
void test_pslq_zeta_weight_homogeneous(int digits, int minWeight, int maxWeight, int maxZetaArg)
{
	std::cout.precision(30);
	struct Job { int w; bool check; };
	struct Mono { int a; std::vector<int> b; }; // pi^(2a) * prod zeta(3+2i)^b[i]
	const int nz = (maxZetaArg >= 3) ? (maxZetaArg - 3) / 2 + 1 : 0;
	auto enumerate = [&](int w) {
		std::vector<Mono> res;
		std::vector<int> b(nz, 0);
		std::function<void(int, int)> rec = [&](int i, int rem) {
			if (i == nz) { if (rem % 2 == 0) res.push_back({ rem / 2, b }); return; }
			const int k = 3 + 2 * i;
			for (int e = 0; e * k <= rem; ++e) { b[i] = e; rec(i + 1, rem - e * k); }
			b[i] = 0;
		};
		rec(0, w);
		return res;
	};
	auto label = [&](const Mono& m) {
		std::string s;
		if (m.a) s += "pi^" + std::to_string(2 * m.a);
		for (int i = 0; i < nz; ++i) if (m.b[i])
			s += (s.empty() ? "" : "*") + std::string("z") + std::to_string(3 + 2 * i) + (m.b[i] > 1 ? "^" + std::to_string(m.b[i]) : "");
		return s;
	};
	std::vector<Job> jobs;
	for (int w = maxWeight; w >= minWeight; --w) { jobs.push_back({ w, false }); jobs.push_back({ w, true }); }
	const int nW = maxWeight + 1;
	std::vector<LinearSearchOutcome> mainOut(nW), checkOut(nW);
	std::vector<int> sizes(nW, 0);
	std::mutex outMutex;
	std::atomic<size_t> next{ 0 };
	const auto t0 = std::chrono::steady_clock::now();
	std::vector<std::thread> threads;
	for (int th = 0; th < 10; ++th)
		threads.emplace_back([&]() {
			for (size_t k; (k = next++) < jobs.size();)
			{
				const Job job = jobs[k];
				const int w = job.w;
				mpreal::set_default_prec(mpfr::digits2bits(digits)); // thread-local
				const std::vector<Mono> monos = enumerate(w);
				sizes[w] = (int)monos.size();
				std::vector<std::string> labels;
				for (auto& m : monos) labels.push_back(label(m));
				const bool even = (w % 2 == 0);
				std::string checkName;
				if (job.check)
				{
					if (!even && monos.size() < 3) continue;
					if (even) { labels.push_back("zeta(" + std::to_string(w) + ")"); checkName = "adds zeta(w); expect rational*pi^w"; }
					else { labels.push_back("y=3m1-7m2+2m_last"); checkName = "adds y; expect (3,-7,..,2,-1)"; }
				}
				auto compute = [&](int) {
					const mpreal pi = mpfr::const_pi();
					std::vector<mpreal> z(nz);
					for (int i = 0; i < nz; ++i) z[i] = mpfr::zeta(mpreal(3 + 2 * i));
					std::vector<mpreal> v;
					for (auto& m : monos)
					{
						mpreal t = mpfr::pow(pi, mpreal(2 * m.a));
						for (int i = 0; i < nz; ++i) if (m.b[i]) t *= mpfr::pow(z[i], mpreal(m.b[i]));
						v.push_back(t);
					}
					if (job.check)
					{
						if (even) v.push_back(mpfr::zeta(mpreal(w)));
						else v.push_back(3 * v[0] - 7 * v[1] + 2 * v[monos.size() - 1]);
					}
					return v;
				};
				std::ostringstream os;
				os.precision(30);
				const std::string name = std::string(job.check ? "CHECK" : "MAIN") + " w=" + std::to_string(w) + " " + checkName;
				const LinearSearchOutcome o = reportLinearRelationSearch(os, name, compute, labels, digits);
				std::lock_guard<std::mutex> lk(outMutex);
				std::cout << os.str() << std::flush;
				(job.check ? checkOut : mainOut)[w] = o;
			}
		});
	for (auto& t : threads) t.join();
	std::cout << "=== SUMMARY (digits = " << digits << ", maxZetaArg = " << maxZetaArg << ")" << std::endl;
	for (int w = minWeight; w <= maxWeight; ++w)
	{
		const bool canCheck = (w % 2 == 0) || sizes[w] >= 3;
		const bool checkOk = canCheck && checkOut[w].found && checkOut[w].verified;
		std::cout << "w=" << w << " n=" << sizes[w] << " power~1e" << (int)(0.8 * digits / std::max(1, sizes[w]))
			<< " check=" << (canCheck ? (checkOk ? "PASS" : "FAIL") : "N/A") << " main=";
		if (canCheck && !checkOk) std::cout << "INVALID";
		else if (mainOut[w].found) std::cout << (mainOut[w].verified ? "RELATION FOUND (VERIFIED)" : "found but spurious");
		else std::cout << "no relation";
		std::cout << std::endl;
	}
	std::cout << "elapsed " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s" << std::endl;
}

// S(s,j,k) = sum_{n>=1} s^n / (n^j C(2n,n)^k); if h2, the extra factor H_{n-1}^{(2)} is inserted.
static mpreal centralBinomialSeries(long s, int j, int k, bool h2 = false)
{
	const mpreal eps = mpfr::pow(mpreal(2), -(long)mpreal::get_default_prec() - 40);
	mpreal sum = 0, q = 1, H = 0;
	for (long n = 1;; ++n)
	{
		// q_n = s^n / C(2n,n)^k
		q *= mpreal(s);
		for (int i = 0; i < k; ++i) q = q * mpreal(n) / mpreal(2 * (2 * n - 1));
		mpreal nj = 1;
		for (int i = 0; i < j; ++i) nj *= mpreal(n);
		if (n > 1) H += mpreal(1) / (mpreal(n - 1) * mpreal(n - 1));
		const mpreal term = h2 ? q * H / nj : q / nj;
		sum += term;
		if (n > 10 && mpfr::abs(term) < eps) break;
	}
	return sum;
}

// zeta(3) (or zeta(5)) against central-binomial series S(s,w,k) and weight-w elementary constants.
// Sanity: V1/V2 (Koecher zeta(3), pi^2/18), V1b (Koecher zeta(5) with H^(2) series), and the main basis once with and once without the Koecher series.
void test_pslq_zeta_central_binomial(int digits, int targetZeta)
{
	std::cout.precision(30);
	const int w = targetZeta;
	struct Ser { long s; int j, k; bool h2; };
	typedef std::vector<std::string> Labels;
	auto lab = [](const Ser& x) {
		return std::string(x.h2 ? "H2*" : "") + "S(" + std::to_string(x.s) + "," + std::to_string(x.j) + "," + std::to_string(x.k) + ")";
	};
	auto zt = [w](int) { return mpfr::zeta(mpreal(w)); };
	// elementary constants of weight w
	Labels eLab;
	std::function<std::vector<mpreal>()> elem;
	if (w == 3)
	{
		eLab = { "pi^3", "pi^2*log2", "log2^3" };
		elem = []() { const mpreal pi = mpfr::const_pi(), l = mpfr::log(mpreal(2)); return std::vector<mpreal>{ pi * pi * pi, pi * pi * l, l * l * l }; };
	}
	else
	{
		eLab = { "pi^2*zeta3", "zeta3*log2^2", "pi^4*log2", "pi^2*log2^3", "log2^5" };
		elem = []() {
			const mpreal pi = mpfr::const_pi(), l = mpfr::log(mpreal(2)), z3 = mpfr::zeta(mpreal(3));
			const mpreal pi2 = pi * pi;
			return std::vector<mpreal>{ pi2 * z3, z3 * l * l, pi2 * pi2 * l, pi2 * l * l * l, l * l * l * l * l };
		};
	}
	auto build = [&](const std::vector<Ser>& sers, bool useElem, bool useZeta) {
		Labels L;
		if (useZeta) L.push_back("zeta(" + std::to_string(w) + ")");
		for (auto& x : sers) L.push_back(lab(x));
		if (useElem) L.insert(L.end(), eLab.begin(), eLab.end());
		auto f = [=](int d) {
			std::vector<mpreal> v;
			if (useZeta) v.push_back(mpfr::zeta(mpreal(w)));
			for (auto& x : sers) v.push_back(centralBinomialSeries(x.s, x.j, x.k, x.h2));
			if (useElem) { auto e = elem(); v.insert(v.end(), e.begin(), e.end()); }
			return v;
		};
		return std::make_pair(L, std::function<std::vector<mpreal>(int)>(f));
	};
	auto run = [&](const std::string& name, const std::vector<Ser>& sers, bool useElem, bool useZeta) {
		auto p = build(sers, useElem, useZeta);
		reportLinearRelationSearch(name, p.second, p.first, digits);
	};
	// sanity
	if (w == 3)
	{
		run("V1 [sanity] zeta3, S(-1,3,1)  expect 2*zeta3 + 5*S(-1,3,1) = 0", { {-1,3,1,false} }, false, true);
		{
			auto L = Labels{ "S(1,2,1)", "pi^2" };
			reportLinearRelationSearch("V2 [sanity] S(1,2,1), pi^2  expect 18*S - pi^2 = 0",
				[](int) { return std::vector<mpreal>{ centralBinomialSeries(1, 2, 1), mpfr::const_pi() * mpfr::const_pi() }; }, L, digits);
		}
	}
	else
	{
		run("V1b [sanity] zeta5, S(-1,5,1), H2*S(-1,3,1)  expect 2*zeta5 + 4*S(-1,5,1) - 5*H2*S(-1,3,1) = 0", { {-1,5,1,false}, {-1,3,1,true} }, false, true);
		{
			auto L = Labels{ "S(1,2,1)", "pi^2" };
			reportLinearRelationSearch("V2 [sanity] S(1,2,1), pi^2  expect 18*S - pi^2 = 0",
				[](int) { return std::vector<mpreal>{ centralBinomialSeries(1, 2, 1), mpfr::const_pi() * mpfr::const_pi() }; }, L, digits);
		}
	}
	// main bases
	std::vector<Ser> base;
	for (long s : { 1, 2, 3, -2, -3 }) base.push_back({ s, w, 1, false });
	for (long s : { 1, -1, 4, -4 }) base.push_back({ s, w, 2, false });
	std::vector<Ser> withK = base;
	withK.push_back({ -1, w, 1, false });
	if (w == 5) withK.push_back({ -1, 3, 1, true });
	run("V3 [sanity] main basis WITH Koecher series: relation must appear", withK, true, true);
	run("M main basis WITHOUT Koecher series", base, true, true);
}

void test_pslq_zeta3_elementary_functions(int digits, int maxDegree)
{
	std::cout.precision(60);
	struct Fn
	{
		std::string name;
		std::function<mpreal(const mpreal&)> f;
		std::function<mpreal(int)> checkC; // known x0 with f(x0) algebraic
		int checkK;
		std::string checkNote;
	};
	auto pi = []() { return mpfr::const_pi(); };
	const std::vector<Fn> fns = {
		{ "exp", [](const mpreal& t) { return mpfr::exp(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "exp(x0)=1+2^(1/3): deg 3" },
		{ "log", [](const mpreal& t) { return mpfr::log(t); }, [](int) { return mpfr::exp(mpfr::sqrt(mpreal(2))); }, 0, "x0=exp(sqrt2): log=sqrt2" },
		{ "sin", [](const mpreal& t) { return mpfr::sin(t); }, [](int) { return mpfr::const_pi() / 5; }, 0, "c=pi/5: deg 4" },
		{ "cos", [](const mpreal& t) { return mpfr::cos(t); }, [](int) { return mpfr::const_pi() / 7; }, 0, "c=pi/7: deg 3" },
		{ "tan", [](const mpreal& t) { return mpfr::tan(t); }, [](int) { return mpfr::const_pi() / 5; }, 0, "c=pi/5: deg 4" },
		{ "sinh", [](const mpreal& t) { return mpfr::sinh(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "c=log(1+2^(1/3)): deg 3" },
		{ "cosh", [](const mpreal& t) { return mpfr::cosh(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "c=log(1+2^(1/3)): deg 3" },
		{ "tanh", [](const mpreal& t) { return mpfr::tanh(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "c=log(1+2^(1/3)): deg 3" },
		{ "arctan", [](const mpreal& t) { return mpfr::atan(t); }, [](int) { return mpfr::tan(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=tan(sqrt2/3): sqrt2/3" },
		{ "arcsinh", [](const mpreal& t) { return mpfr::asinh(t); }, [](int) { return mpfr::sinh(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=sinh(sqrt2/3): sqrt2/3" },
		{ "arccosh", [](const mpreal& t) { return mpfr::acosh(t); }, [](int) { return mpfr::cosh(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=cosh(sqrt2/3): sqrt2/3" },
		{ "arcsin(1/t)", [](const mpreal& t) { return mpfr::asin(1 / t); }, [](int) { return 1 / mpfr::sin(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=1/sin(sqrt2/3): sqrt2/3" },
		{ "arccos(1/t)", [](const mpreal& t) { return mpfr::acos(1 / t); }, [](int) { return 1 / mpfr::cos(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=1/cos(sqrt2/3): sqrt2/3" },
		{ "log(t-1)", [](const mpreal& t) { return mpfr::log(t - 1); }, [](int) { return 1 + mpfr::exp(mpfr::sqrt(mpreal(2))); }, 0, "x0=1+exp(sqrt2): sqrt2" },
		{ "t^t", [](const mpreal& t) { return mpfr::pow(t, t); }, [](int) {
				// solve c^c = 2 by Newton on g(c) = c log c - log 2
				const mpreal l2 = mpfr::log(mpreal(2));
				mpreal c = 1.5;
				for (int i = 0; i < 200; ++i)
					c -= (c * mpfr::log(c) - l2) / (mpfr::log(c) + 1);
				return c; }, 0, "c^c=2: rational 2" },
		{ "Gamma", [](const mpreal& t) { return mpfr::gamma(t); }, [](int) { return mpreal(4); }, 0, "c=4: Gamma=6 (weak check)" },
		{ "exp(1/t)", [](const mpreal& t) { return mpfr::exp(1 / t); }, [](int) { return 1 / mpfr::log(mpreal(2)); }, 0, "c=1/log2: 2" },
	};
	(void)pi;
	const std::vector<int> ks = { -2, -1, 0, 1, 2, 3 }; // targets f(zeta(3) / pi^k)
	std::vector<PolySearchTarget> targets;
	for (const auto& fn : fns)
		for (int k : ks)
		{
			mpreal::set_default_prec(mpfr::digits2bits(50));
			if (mpfr::isnan(fn.f(mpfr::zeta(mpreal(3)) / mpfr::pow(mpfr::const_pi(), mpreal(k)))))
			{
				std::cout << "skip " << fn.name << "(zeta3 / pi^" << k << "): outside the domain of the function" << std::endl;
				continue;
			}
			targets.push_back({ fn.name + "(zeta3 / pi^" + std::to_string(k) + ")", [fn, k](int d)
				{
					mpreal::set_default_prec(mpfr::digits2bits(d));
					return fn.f(mpfr::zeta(mpreal(3)) / mpfr::pow(mpfr::const_pi(), mpreal(k)));
				} });
		}
	const size_t nMain = targets.size();
	for (const auto& fn : fns)
		targets.push_back({ "[CHECK] " + fn.name + "(c / pi^2), c = x0 pi^2  (" + fn.checkNote + ")", [fn](int d)
			{
				mpreal::set_default_prec(mpfr::digits2bits(d));
				const mpreal pi2 = mpfr::const_pi() * mpfr::const_pi();
				const mpreal c = fn.checkC(d) * pi2; // known x0 times pi^2, so that c / pi^2 = x0
				return fn.f(c / pi2);
			} });
	std::cout << "NOTE: targets are f(zeta(3)/pi^k). [CHECK] targets replace zeta(3)/pi^k by a known x0 (built as c/pi^2 with c = x0 pi^2) with f(x0) algebraic. "
		"t^t uses x0^x0=2, Gamma uses x0=4 (weak), exp(1/t) uses x0=1/log2." << std::endl;

	const auto sm = reportPolynomialSearchParallel(targets, digits, maxDegree);

	std::cout << "\n===== SUMMARY (" << digits << " digits, degree <= " << maxDegree << ") =====" << std::endl;
	std::cout << "target | verified relations (first degree) | norm bound at max degree / degree 4 / degree 8" << std::endl;
	for (size_t t = 0; t < targets.size(); ++t)
	{
		if (t == nMain)
			std::cout << "--- validation cases ---" << std::endl;
		int firstDeg = 0;
		for (int d = 1; d <= maxDegree && !firstDeg; ++d)
			if (sm.results[t][d].status == PslqResult::Status::Found)
				firstDeg = d;
		std::cout << targets[t].name << " | ";
		if (sm.verified[t] > 0)
			std::cout << "FOUND (verified), first at degree " << firstDeg;
		else if (firstDeg)
			std::cout << "SPURIOUS only";
		else
			std::cout << "none";
		auto nb = [&](int d) { return d <= maxDegree ? sm.results[t][d].normBound.toString("%.2Re") : std::string("n/a"); };
		if (firstDeg == 0)
			std::cout << " | " << nb(maxDegree) << " / " << nb(4) << " / " << nb(8);
		std::cout << std::endl;
	}
}

void test_pslq_zeta3_elementary_functions_relations(int digits)
{
	std::cout.precision(30);
	struct Fn
	{
		std::string name;
		std::function<mpreal(const mpreal&)> f;
		std::function<mpreal(const mpreal&)> inv; // x0 with f(x0) = v
		bool small; // true: value range is bounded/small, use v/4
	};
	auto newton = [](const std::function<mpreal(const mpreal&)>& g, const std::function<mpreal(const mpreal&)>& dg, mpreal x)
	{
		const mpreal eps = mpfr::pow(mpreal(2), -(long)mpreal::get_default_prec() + 20);
		for (int i = 0; i < 200; ++i)
		{
			const mpreal dx = g(x) / dg(x);
			x -= dx;
			if (mpfr::abs(dx) < eps) break;
		}
		return x;
	};
	const std::vector<Fn> fns = {
		{ "exp", [](const mpreal& t) { return mpfr::exp(t); }, [](const mpreal& v) { return mpfr::log(v); }, false },
		{ "log", [](const mpreal& t) { return mpfr::log(t); }, [](const mpreal& v) { return mpfr::exp(v); }, false },
		{ "sin", [](const mpreal& t) { return mpfr::sin(t); }, [](const mpreal& v) { return mpfr::asin(v); }, true },
		{ "cos", [](const mpreal& t) { return mpfr::cos(t); }, [](const mpreal& v) { return mpfr::acos(v); }, true },
		{ "tan", [](const mpreal& t) { return mpfr::tan(t); }, [](const mpreal& v) { return mpfr::atan(v); }, true },
		{ "sinh", [](const mpreal& t) { return mpfr::sinh(t); }, [](const mpreal& v) { return mpfr::asinh(v); }, false },
		{ "cosh", [](const mpreal& t) { return mpfr::cosh(t); }, [](const mpreal& v) { return mpfr::acosh(v); }, false },
		{ "tanh", [](const mpreal& t) { return mpfr::tanh(t); }, [](const mpreal& v) { return mpfr::atanh(v); }, true },
		{ "arctan", [](const mpreal& t) { return mpfr::atan(t); }, [](const mpreal& v) { return mpfr::tan(v); }, true },
		{ "arcsinh", [](const mpreal& t) { return mpfr::asinh(t); }, [](const mpreal& v) { return mpfr::sinh(v); }, false },
		{ "arccosh", [](const mpreal& t) { return mpfr::acosh(t); }, [](const mpreal& v) { return mpfr::cosh(v); }, false },
		{ "arcsin(1/t)", [](const mpreal& t) { return mpfr::asin(1 / t); }, [](const mpreal& v) { return 1 / mpfr::sin(v); }, true },
		{ "arccos(1/t)", [](const mpreal& t) { return mpfr::acos(1 / t); }, [](const mpreal& v) { return 1 / mpfr::cos(v); }, true },
		{ "log(t-1)", [](const mpreal& t) { return mpfr::log(t - 1); }, [](const mpreal& v) { return 1 + mpfr::exp(v); }, false },
		{ "t^t", [](const mpreal& t) { return mpfr::pow(t, t); }, [newton](const mpreal& v) {
				const mpreal lv = mpfr::log(v);
				return newton([&](const mpreal& c) { return c * mpfr::log(c) - lv; }, [](const mpreal& c) { return mpfr::log(c) + 1; }, mpreal(2)); }, false },
		{ "Gamma", [](const mpreal& t) { return mpfr::gamma(t); }, [newton](const mpreal& v) {
				const mpreal lv = mpfr::log(v);
				return newton([&](const mpreal& c) { return mpfr::lngamma(c) - lv; }, [](const mpreal& c) { return mpfr::digamma(c); }, mpreal(3)); }, false },
		{ "exp(1/t)", [](const mpreal& t) { return mpfr::exp(1 / t); }, [](const mpreal& v) { return 1 / mpfr::log(v); }, false },
	};
	const std::vector<std::string> labelsA = { "log|y|", "1", "log(pi)", "log2", "log3", "log5", "log7", "log11", "log13", "log(gamma)",
		"log(1+sqrt2)", "log(phi)", "log(2+sqrt3)" };
	const std::vector<std::string> labelsB = { "y", "1", "pi", "pi^2", "pi^3", "log2", "log3", "e", "gamma", "G", "sqrt2", "sqrt3", "sqrt5", "zeta3" };
	auto vecA = [](const mpreal& y)
	{
		return std::vector<mpreal>{ mpfr::log(mpfr::abs(y)), mpreal(1), mpfr::log(mpfr::const_pi()), mpfr::log(mpreal(2)), mpfr::log(mpreal(3)),
			mpfr::log(mpreal(5)), mpfr::log(mpreal(7)), mpfr::log(mpreal(11)), mpfr::log(mpreal(13)), mpfr::log(mpfr::const_euler()),
			mpfr::log(1 + mpfr::sqrt(mpreal(2))), mpfr::log((1 + mpfr::sqrt(mpreal(5))) / 2), mpfr::log(2 + mpfr::sqrt(mpreal(3))) };
	};
	auto vecB = [](const mpreal& y)
	{
		const mpreal pi = mpfr::const_pi();
		return std::vector<mpreal>{ y, mpreal(1), pi, pi * pi, pi * pi * pi, mpfr::log(mpreal(2)), mpfr::log(mpreal(3)), mpfr::exp(mpreal(1)),
			mpfr::const_euler(), mpfr::const_catalan(), mpfr::sqrt(mpreal(2)), mpfr::sqrt(mpreal(3)), mpfr::sqrt(mpreal(5)), mpfr::zeta(mpreal(3)) };
	};
	auto valA = [](bool small) { return mpfr::cbrt(mpreal(2)) * mpfr::sqrt(mpfr::const_pi()) * mpfr::exp(mpreal(1) / 5) / (small ? 4 : 1); };
	auto valB = [](bool small)
	{
		return (2 * mpfr::const_pi() - 3 - mpfr::log(mpreal(2)) / 5 - mpfr::sqrt(mpreal(3)) / 7) / (small ? 4 : 1);
	};

	struct Job
	{
		std::string name;
		bool isCheck, isA;
		std::function<std::vector<mpreal>(int)> compute;
		const std::vector<std::string>* labels;
		std::string text;
		LinearSearchOutcome oc;
		bool skipped = false;
	};
	std::vector<Job> jobs;
	const std::vector<int> ks = { -3, -2, -1, 0 };
	std::cout << "NOTE: targets y = f(zeta3 * pi^k), k in {-3,-2,-1,0}. (A) multiplicative basis n=13, (B) additive basis n=14. "
		"[CHECK] jobs use x0 = f^-1(v) with v = 2^(1/3) pi^(1/2) e^(1/5) [/4 for bounded-range f] (A), v = 2pi-3-log2/5-sqrt3/7 [/4] (B)." << std::endl;
	for (const auto& fn : fns)
	{
		for (int k : ks)
		{
			mpreal::set_default_prec(mpfr::digits2bits(50));
			const mpreal y50 = fn.f(mpfr::zeta(mpreal(3)) * mpfr::pow(mpfr::const_pi(), mpreal(k)));
			if (mpfr::isnan(y50) || mpfr::isinf(y50))
			{
				std::cout << "skip " << fn.name << "(zeta3 * pi^" << k << "): outside the domain" << std::endl;
				continue;
			}
			auto yOf = [fn, k](int d)
			{
				mpreal::set_default_prec(mpfr::digits2bits(d));
				return fn.f(mpfr::zeta(mpreal(3)) * mpfr::pow(mpfr::const_pi(), mpreal(k)));
			};
			const std::string nm = fn.name + "(zeta3*pi^" + std::to_string(k) + ")";
			jobs.push_back({ "(A) " + nm, false, true, [yOf, vecA](int d) { return vecA(yOf(d)); }, &labelsA });
			jobs.push_back({ "(B) " + nm, false, false, [yOf, vecB](int d) { return vecB(yOf(d)); }, &labelsB });
		}
		const bool sm = fn.small;
		auto chk = [fn, sm](int d, bool isA, const std::function<mpreal(bool)>& val)
		{
			mpreal::set_default_prec(mpfr::digits2bits(d));
			(void)isA;
			return fn.f(fn.inv(val(sm)));
		};
		jobs.push_back({ "(A) [CHECK] " + fn.name + "(x0), f(x0) = " + std::string(sm ? "2^(1/3)pi^(1/2)e^(1/5)/4" : "2^(1/3)pi^(1/2)e^(1/5)"), true, true,
			[chk, vecA, valA](int d) { return vecA(chk(d, true, valA)); }, &labelsA });
		jobs.push_back({ "(B) [CHECK] " + fn.name + "(x0), f(x0) = " + std::string(sm ? "(2pi-3-log2/5-sqrt3/7)/4" : "2pi-3-log2/5-sqrt3/7"), true, false,
			[chk, vecB, valB](int d) { return vecB(chk(d, false, valB)); }, &labelsB });
	}

	std::mutex mtx;
	std::atomic<size_t> next(0);
	const auto t0 = std::chrono::steady_clock::now();
	auto worker = [&]()
	{
		for (;;)
		{
			const size_t i = next++;
			if (i >= jobs.size()) return;
			std::ostringstream os;
			os.precision(30);
			Job& j = jobs[i];
			j.oc = reportLinearRelationSearch(os, j.name, j.compute, *j.labels, digits);
			j.text = os.str();
			std::lock_guard<std::mutex> lk(mtx);
			std::cout << j.text << std::endl;
		}
	};
	std::vector<std::thread> th;
	for (int t = 0; t < 10; ++t) th.emplace_back(worker);
	for (auto& t : th) t.join();

	std::cout << "\n===== SUMMARY (" << digits << " digits) =====" << std::endl;
	// validity per (function, type)
	int nFail = 0;
	for (const auto& j : jobs)
		if (j.isCheck && !(j.oc.found && j.oc.verified)) { ++nFail; std::cout << "CHECK FAIL: " << j.name << std::endl; }
	std::cout << "validation: " << (jobs.size() ? "" : "") << (nFail == 0 ? "all PASS" : "some FAIL") << std::endl;
	for (const auto& j : jobs)
	{
		if (j.isCheck)
		{
			std::cout << (j.oc.found && j.oc.verified ? "PASS " : "FAIL ") << j.name << std::endl;
			continue;
		}
		// find corresponding check
		bool valid = true;
		const std::string fname = j.name.substr(4, j.name.find('(', 4) - 4);
		for (const auto& c : jobs)
			if (c.isCheck && c.isA == j.isA && c.name.compare(4 + 8, fname.size() + 1, fname + "(") == 0 && !(c.oc.found && c.oc.verified)) valid = false;
		std::cout << j.name << " | ";
		if (!valid) std::cout << "INVALID (check failed)";
		else if (j.oc.found) std::cout << (j.oc.verified ? "FOUND (verified)" : "found but NOT verified");
		else std::cout << "none, normBound " << j.oc.normBound;
		std::cout << std::endl;
	}
	std::cout << "total elapsed " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s" << std::endl;
}

// Compositions y = f(g(zeta3 * pi^k)), k in {-3..0}, f,g in 12 elementary functions (inverse pairs and f = g excluded).
// (1) doAlgebraic: y algebraic of degree <= maxDegree?  (2) doRelations: multiplicative (A) / additive (B) PSLQ vectors, as in
// test_pslq_zeta3_elementary_functions_relations.
// Env ZCOMP_STRIDE=N (smoke test only): keep every N-th main target.
void test_pslq_zeta3_composite_functions(int digits, int maxDegree, bool doAlgebraic, bool doRelations)
{
	std::cout.precision(30);
	struct Fn
	{
		std::string name;
		std::function<mpreal(const mpreal&)> f;
		std::function<mpreal(const mpreal&)> inv;
		bool small;
	};
	const std::vector<Fn> fns = {
		{ "exp", [](const mpreal& t) { return mpfr::exp(t); }, [](const mpreal& v) { return mpfr::log(v); }, false },
		{ "log", [](const mpreal& t) { return mpfr::log(t); }, [](const mpreal& v) { return mpfr::exp(v); }, false },
		{ "sin", [](const mpreal& t) { return mpfr::sin(t); }, [](const mpreal& v) { return mpfr::asin(v); }, true },
		{ "cos", [](const mpreal& t) { return mpfr::cos(t); }, [](const mpreal& v) { return mpfr::acos(v); }, true },
		{ "tan", [](const mpreal& t) { return mpfr::tan(t); }, [](const mpreal& v) { return mpfr::atan(v); }, true },
		{ "sinh", [](const mpreal& t) { return mpfr::sinh(t); }, [](const mpreal& v) { return mpfr::asinh(v); }, false },
		{ "cosh", [](const mpreal& t) { return mpfr::cosh(t); }, [](const mpreal& v) { return mpfr::acosh(v); }, false },
		{ "tanh", [](const mpreal& t) { return mpfr::tanh(t); }, [](const mpreal& v) { return mpfr::atanh(v); }, true },
		{ "arctan", [](const mpreal& t) { return mpfr::atan(t); }, [](const mpreal& v) { return mpfr::tan(v); }, true },
		{ "arcsinh", [](const mpreal& t) { return mpfr::asinh(t); }, [](const mpreal& v) { return mpfr::sinh(v); }, false },
		{ "arcsin", [](const mpreal& t) { return mpfr::asin(t); }, [](const mpreal& v) { return mpfr::sin(v); }, true },
		{ "arccos", [](const mpreal& t) { return mpfr::acos(t); }, [](const mpreal& v) { return mpfr::cos(v); }, true },
	};
	auto isInversePair = [](const std::string& a, const std::string& b)
	{
		static const std::vector<std::pair<std::string, std::string>> inv = { { "exp", "log" }, { "sin", "arcsin" }, { "tan", "arctan" }, { "sinh", "arcsinh" }, { "cos", "arccos" } };
		for (auto& p : inv) if ((p.first == a && p.second == b) || (p.first == b && p.second == a)) return true;
		return false;
	};
	const std::vector<int> ks = { -3, -2, -1, 0 };
	auto tOf = [](int k) { return mpfr::zeta(mpreal(3)) * mpfr::pow(mpfr::const_pi(), mpreal(k)); };
	struct Tgt { std::string name; std::function<mpreal(int)> y; };
	std::vector<Tgt> main;
	int nDomain = 0, nRange = 0;
	for (const auto& f : fns)
		for (const auto& g : fns)
		{
			if (f.name == g.name || isInversePair(f.name, g.name)) continue;
			for (int k : ks)
			{
				mpreal::set_default_prec(mpfr::digits2bits(50));
				const mpreal y50 = f.f(g.f(tOf(k)));
				if (mpfr::isnan(y50) || mpfr::isinf(y50)) { ++nDomain; continue; }
				const mpreal ay = mpfr::abs(y50);
				if (ay < mpreal("1e-30") || ay > mpreal("1e30")) { ++nRange; continue; }
				main.push_back({ f.name + "(" + g.name + "(zeta3*pi^" + std::to_string(k) + "))", [f, g, k, tOf](int d)
					{
						mpreal::set_default_prec(mpfr::digits2bits(d));
						return f.f(g.f(tOf(k)));
					} });
			}
		}
	if (const char* s = std::getenv("ZCOMP_STRIDE"))
	{
		const size_t st = std::max(1, std::atoi(s));
		std::vector<Tgt> kept;
		for (size_t i = 0; i < main.size(); i += st) kept.push_back(main[i]);
		main = kept;
		std::cout << "SMOKE: ZCOMP_STRIDE=" << st << std::endl;
	}
	std::cout << "NOTE: targets y = f(g(zeta3*pi^k)), k in {-3..0}; dropped " << nDomain << " (domain error), " << nRange << " (|y|<1e-30 or >1e30). "
		"FINAL TARGET COUNT = " << main.size() << std::endl;

	const auto t0 = std::chrono::steady_clock::now();
	if (doAlgebraic)
	{
		std::vector<PolySearchTarget> targets;
		for (auto& m : main) targets.push_back({ m.name, m.y });
		const size_t nMain = targets.size();
		auto addChk = [&](const std::string& nm, std::function<mpreal()> v)
		{
			targets.push_back({ "[CHECK] " + nm, [v](int d) { mpreal::set_default_prec(mpfr::digits2bits(d)); return v(); } });
		};
		addChk("cos(arcsin(1/3)) = sqrt8/3 (deg 2)", []() { return mpfr::cos(mpfr::asin(mpreal(1) / 3)); });
		addChk("sin(arctan(2)) = 2/sqrt5 (deg 2)", []() { return mpfr::sin(mpfr::atan(mpreal(2))); });
		addChk("cosh(arcsinh(sqrt3)) = 2 (deg 1)", []() { return mpfr::cosh(mpfr::asinh(mpfr::sqrt(mpreal(3)))); });
		addChk("cos(arctan(3^(1/4))) = 1/sqrt(1+sqrt3) (deg 4)", []() { return mpfr::cos(mpfr::atan(mpfr::root(mpreal(3), 4))); });
		addChk("tanh(log 5) = 12/13 (deg 1)", []() { return mpfr::tanh(mpfr::log(mpreal(5))); });
		addChk("sin(arccos(2^(1/3)/2)) (deg 6)", []() { return mpfr::sin(mpfr::acos(mpfr::cbrt(mpreal(2)) / 2)); });
		const std::vector<int> expectDeg = { 2, 2, 1, 4, 1, 6 };

		const auto sm = reportPolynomialSearchParallel(targets, digits, maxDegree);

		std::cout << "\n===== ALGEBRAIC SUMMARY (" << digits << " digits, degree <= " << maxDegree << ", " << nMain << " targets) =====" << std::endl;
		int nFound = 0;
		for (size_t t = 0; t < nMain; ++t)
		{
			int firstDeg = 0;
			for (int d = 1; d <= maxDegree && !firstDeg; ++d)
				if (sm.results[t][d].status == PslqResult::Status::Found) firstDeg = d;
			if (firstDeg) { ++nFound; std::cout << targets[t].name << " | " << (sm.verified[t] > 0 ? "FOUND (verified)" : "SPURIOUS only") << ", first at degree " << firstDeg << std::endl; }
		}
		std::cout << "main targets with a relation: " << nFound << std::endl;
		std::cout << "min normBound per degree over all main targets without relation:" << std::endl;
		for (int d = 1; d <= maxDegree; ++d)
		{
			bool have = false; mpreal mn; std::string who;
			for (size_t t = 0; t < nMain; ++t)
			{
				bool anyFound = false;
				for (int e = 1; e <= maxDegree; ++e) if (sm.results[t][e].status == PslqResult::Status::Found) anyFound = true;
				if (anyFound) continue;
				const mpreal nb = sm.results[t][d].normBound;
				if (!have || nb < mn) { mn = nb; who = targets[t].name; have = true; }
			}
			if (have) std::cout << "  degree " << d << ": " << mn.toString("%.3Re") << "  (" << who << ")" << std::endl;
		}
		int nFail = 0;
		for (size_t c = 0; c < expectDeg.size(); ++c)
		{
			const size_t t = nMain + c;
			int firstDeg = 0;
			for (int d = 1; d <= maxDegree && !firstDeg; ++d)
				if (sm.results[t][d].status == PslqResult::Status::Found) firstDeg = d;
			if (expectDeg[c] > maxDegree) { std::cout << "SKIP (needs maxDegree >= " << expectDeg[c] << ") " << targets[t].name << std::endl; continue; }
			const bool ok = firstDeg == expectDeg[c] && sm.verified[t] > 0;
			if (!ok) ++nFail;
			std::cout << (ok ? "PASS " : "FAIL ") << targets[t].name << " : first degree " << firstDeg << ", verified " << sm.verified[t] << std::endl;
		}
		std::cout << "algebraic validation: " << (nFail == 0 ? "all PASS" : "some FAIL") << std::endl;
		std::cout << "elapsed (algebraic) " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s" << std::endl;
	}

	if (doRelations)
	{
		const auto t1 = std::chrono::steady_clock::now();
		const std::vector<std::string> labelsA = { "log|y|", "1", "log(pi)", "log2", "log3", "log5", "log7", "log11", "log13", "log(gamma)",
			"log(1+sqrt2)", "log(phi)", "log(2+sqrt3)" };
		const std::vector<std::string> labelsB = { "y", "1", "pi", "pi^2", "pi^3", "log2", "log3", "e", "gamma", "G", "sqrt2", "sqrt3", "sqrt5", "zeta3" };
		auto vecA = [](const mpreal& y)
		{
			return std::vector<mpreal>{ mpfr::log(mpfr::abs(y)), mpreal(1), mpfr::log(mpfr::const_pi()), mpfr::log(mpreal(2)), mpfr::log(mpreal(3)),
				mpfr::log(mpreal(5)), mpfr::log(mpreal(7)), mpfr::log(mpreal(11)), mpfr::log(mpreal(13)), mpfr::log(mpfr::const_euler()),
				mpfr::log(1 + mpfr::sqrt(mpreal(2))), mpfr::log((1 + mpfr::sqrt(mpreal(5))) / 2), mpfr::log(2 + mpfr::sqrt(mpreal(3))) };
		};
		auto vecB = [](const mpreal& y)
		{
			const mpreal pi = mpfr::const_pi();
			return std::vector<mpreal>{ y, mpreal(1), pi, pi * pi, pi * pi * pi, mpfr::log(mpreal(2)), mpfr::log(mpreal(3)), mpfr::exp(mpreal(1)),
				mpfr::const_euler(), mpfr::const_catalan(), mpfr::sqrt(mpreal(2)), mpfr::sqrt(mpreal(3)), mpfr::sqrt(mpreal(5)), mpfr::zeta(mpreal(3)) };
		};
		auto valA = [](bool small) { return mpfr::cbrt(mpreal(2)) * mpfr::sqrt(mpfr::const_pi()) * mpfr::exp(mpreal(1) / 5) / (small ? 4 : 1); };
		auto valB = [](bool small) { return (2 * mpfr::const_pi() - 3 - mpfr::log(mpreal(2)) / 5 - mpfr::sqrt(mpreal(3)) / 7) / (small ? 4 : 1); };

		struct Job
		{
			std::string name;
			bool isCheck, isA;
			std::function<std::vector<mpreal>(int)> compute;
			const std::vector<std::string>* labels;
			std::string text;
			LinearSearchOutcome oc;
		};
		std::vector<Job> jobs;
		for (auto& m : main)
		{
			auto yOf = m.y;
			jobs.push_back({ "(A) " + m.name, false, true, [yOf, vecA](int d) { return vecA(yOf(d)); }, &labelsA });
			jobs.push_back({ "(B) " + m.name, false, false, [yOf, vecB](int d) { return vecB(yOf(d)); }, &labelsB });
		}
		// single-function checks (same as the single-function test)
		for (const auto& fn : fns)
		{
			const bool sm = fn.small;
			auto chk = [fn, sm](int d, const std::function<mpreal(bool)>& val)
			{
				mpreal::set_default_prec(mpfr::digits2bits(d));
				return fn.f(fn.inv(val(sm)));
			};
			jobs.push_back({ "(A) [CHECK] " + fn.name + "(x0), f(x0) = 2^(1/3)pi^(1/2)e^(1/5)" + (sm ? "/4" : ""), true, true,
				[chk, vecA, valA](int d) { return vecA(chk(d, valA)); }, &labelsA });
			jobs.push_back({ "(B) [CHECK] " + fn.name + "(x0), f(x0) = 2pi-3-log2/5-sqrt3/7" + (sm ? " over 4" : ""), true, false,
				[chk, vecB, valB](int d) { return vecB(chk(d, valB)); }, &labelsB });
		}
		// composite checks: y = f(g(x0)) with x0 = g^-1(f^-1(v))
		for (const auto& f : fns)
			for (const auto& g : fns)
			{
				const bool pickA = (f.name == "cosh" && g.name == "arcsinh") || (f.name == "sin" && g.name == "arctan") || (f.name == "exp" && g.name == "arcsinh");
				if (!pickA) continue;
				const bool sm = f.small;
				auto chk = [f, g, sm](int d, const std::function<mpreal(bool)>& val)
				{
					mpreal::set_default_prec(mpfr::digits2bits(d));
					return f.f(g.f(g.inv(f.inv(val(sm)))));
				};
				jobs.push_back({ "(A) [CHECK] composite " + f.name + "(" + g.name + "(x0))", true, true, [chk, vecA, valA](int d) { return vecA(chk(d, valA)); }, &labelsA });
				jobs.push_back({ "(B) [CHECK] composite " + f.name + "(" + g.name + "(x0))", true, false, [chk, vecB, valB](int d) { return vecB(chk(d, valB)); }, &labelsB });
			}
		// fixed known composites (B): sinh(log 2) = 3/4 ; log(sinh(log 2)) = log3 - 2 log2 ; log(cosh(log 3)) = log5 - log3 (A: y = log(5/3) not in basis; tested as B-like via explicit relation below)
		jobs.push_back({ "(B) [CHECK] sinh(log 2) = 3/4", true, false, [vecB](int d) { mpreal::set_default_prec(mpfr::digits2bits(d)); return vecB(mpfr::sinh(mpfr::log(mpreal(2)))); }, &labelsB });
		jobs.push_back({ "(B) [CHECK] log(sinh(log 2)) = log3 - 2log2", true, false, [vecB](int d) { mpreal::set_default_prec(mpfr::digits2bits(d)); return vecB(mpfr::log(mpfr::sinh(mpfr::log(mpreal(2))))); }, &labelsB });
		jobs.push_back({ "(B') [CHECK] log(cosh(log 3)) = log5 - log3 (vector y, 1, log3, log5)", true, false, [](int d) {
			mpreal::set_default_prec(mpfr::digits2bits(d));
			return std::vector<mpreal>{ mpfr::log(mpfr::cosh(mpfr::log(mpreal(3)))), mpreal(1), mpfr::log(mpreal(3)), mpfr::log(mpreal(5)), mpfr::const_pi() }; },
			&labelsB });
		jobs.back().labels = nullptr;
		static const std::vector<std::string> labelsBp = { "y", "1", "log3", "log5", "pi" };
		jobs.back().labels = &labelsBp;

		std::cout << "RELATIONS: " << jobs.size() << " jobs (" << main.size() << " main targets x 2 + checks); (A) n=13 multiplicative, (B) n=14 additive." << std::endl;
		std::mutex mtx;
		std::atomic<size_t> next(0);
		auto worker = [&]()
		{
			for (;;)
			{
				const size_t i = next++;
				if (i >= jobs.size()) return;
				std::ostringstream os;
				os.precision(30);
				Job& j = jobs[i];
				const auto tj = std::chrono::steady_clock::now();
				j.oc = reportLinearRelationSearch(os, j.name, j.compute, *j.labels, digits);
				os << "elapsed " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tj).count() << " s" << std::endl;
				j.text = os.str();
				std::lock_guard<std::mutex> lk(mtx);
				std::cout << j.text << std::endl;
			}
		};
		std::vector<std::thread> th;
		for (int t = 0; t < 10; ++t) th.emplace_back(worker);
		for (auto& t : th) t.join();

		std::cout << "\n===== RELATIONS SUMMARY (" << digits << " digits) =====" << std::endl;
		int nFail = 0, nFound = 0;
		for (const auto& j : jobs)
			if (j.isCheck)
			{
				const bool ok = j.oc.found && j.oc.verified;
				if (!ok) ++nFail;
				std::cout << (ok ? "PASS " : "FAIL ") << j.name << std::endl;
			}
		std::cout << "relations validation: " << (nFail == 0 ? "all PASS" : "some FAIL") << std::endl;
		for (const auto& j : jobs)
			if (!j.isCheck && j.oc.found)
			{
				++nFound;
				std::cout << j.name << " | " << (j.oc.verified ? "FOUND (verified)" : "found but NOT verified") << std::endl;
			}
		std::cout << "main jobs with a relation: " << nFound << std::endl;
		// min normBound over A and B separately (normBound is a string)
		for (int isA = 1; isA >= 0; --isA)
		{
			bool have = false; mpreal mn; std::string who;
			for (const auto& j : jobs)
				if (!j.isCheck && !j.oc.found && (j.isA == (bool)isA) && !j.oc.normBound.empty())
				{
					mpreal nb(j.oc.normBound);
					if (!have || nb < mn) { mn = nb; who = j.name; have = true; }
				}
			if (have) std::cout << "min normBound " << (isA ? "(A)" : "(B)") << ": " << mn.toString("%.3Re") << "  (" << who << ")" << std::endl;
		}
		std::cout << "elapsed (relations) " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count() << " s" << std::endl;
	}
	std::cout << "total elapsed " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s" << std::endl;
}


namespace
{
	// Shared by the *_highdeg wrappers. Prints per-target status, summary, min normBound, and validation results.
	struct HdCheck
	{
		size_t index;      // index into targets
		int minTrueDeg;    // required true degree of the found relation (>= ; 0 = any)
	};

	void reportHighDegSummary(const std::vector<PolySearchTarget>& targets, size_t nMain, const PolySearchSummary& sm,
		const std::vector<int>& checkMinDeg, int digits, int degree)
	{
		auto trueDegOf = [](const PslqResult& r)
		{
			int td = 0;
			for (size_t i = 0; i < r.relation.size(); ++i) if (r.relation[i] != 0) td = static_cast<int>(i);
			return td;
		};
		std::cout << "\n===== HIGHDEG SUMMARY (" << digits << " digits, only degree " << degree << " [covers all degrees <= " << degree << "], "
			<< nMain << " main targets) =====" << std::endl;
		int nFound = 0;
		bool have = false; mpreal mn; std::string who;
		for (size_t t = 0; t < nMain; ++t)
		{
			const PslqResult& r = sm.results[t][degree];
			const bool found = r.status == PslqResult::Status::Found;
			std::cout << targets[t].name << " | " << toString(r.status) << " | ";
			if (found)
			{
				++nFound;
				std::cout << (sm.verified[t] > 0 ? "FOUND (verified)" : "SPURIOUS") << ", true degree " << trueDegOf(r) << std::endl;
			}
			else
			{
				std::cout << "normBound " << r.normBound.toString("%.3Re") << std::endl;
				if (!have || r.normBound < mn) { mn = r.normBound; who = targets[t].name; have = true; }
			}
		}
		std::cout << "main targets with a relation: " << nFound << " / " << nMain << std::endl;
		if (have) std::cout << "MIN normBound over targets without relation: " << mn.toString("%.3Re") << "  (" << who << ")" << std::endl;
		int nFail = 0;
		for (size_t c = 0; c < checkMinDeg.size(); ++c)
		{
			const size_t t = nMain + c;
			const PslqResult& r = sm.results[t][degree];
			const bool found = r.status == PslqResult::Status::Found;
			const int td = found ? trueDegOf(r) : -1;
			const bool ok = found && sm.verified[t] > 0 && td >= checkMinDeg[c];
			if (!ok) ++nFail;
			std::cout << (ok ? "PASS " : "FAIL ") << targets[t].name << " : status " << toString(r.status) << ", true degree " << td
				<< " (required >= " << checkMinDeg[c] << "), verified " << sm.verified[t] << std::endl;
		}
		std::cout << "validation: " << (nFail == 0 ? "all PASS" : "some FAIL") << " (" << checkMinDeg.size() << " cases)" << std::endl;
	}

	void applyStrideEnv(std::vector<PolySearchTarget>& main)
	{
		if (const char* s = std::getenv("ZHD_STRIDE"))
		{
			const size_t st = std::max(1, std::atoi(s));
			std::vector<PolySearchTarget> kept;
			for (size_t i = 0; i < main.size(); i += st) kept.push_back(main[i]);
			main = kept;
			std::cout << "SMOKE: ZHD_STRIDE=" << st << std::endl;
		}
	}
}

// Single functions f(t), t = zeta3 * pi^k (k in ks; NOTE: opposite sign convention to test_pslq_zeta3_elementary_functions),
// algebraic of degree <= `degree`?  Only degree == `degree` is run (the vector (1,..,x^d) contains all lower degrees).
void test_pslq_zeta3_elementary_functions_highdeg(int digits, int degree, std::vector<int> ks)
{
	std::cout.precision(60);
	struct Fn
	{
		std::string name;
		std::function<mpreal(const mpreal&)> f;
		std::function<mpreal(int)> checkC; // known x0 with f(x0) algebraic
		int checkK;
		std::string checkNote;
	};
	const std::vector<Fn> fns = {
		{ "exp", [](const mpreal& t) { return mpfr::exp(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "exp(x0)=1+2^(1/3): deg 3" },
		{ "log", [](const mpreal& t) { return mpfr::log(t); }, [](int) { return mpfr::exp(mpfr::sqrt(mpreal(2))); }, 0, "x0=exp(sqrt2): log=sqrt2" },
		{ "sin", [](const mpreal& t) { return mpfr::sin(t); }, [](int) { return mpfr::const_pi() / 5; }, 0, "c=pi/5: deg 4" },
		{ "cos", [](const mpreal& t) { return mpfr::cos(t); }, [](int) { return mpfr::const_pi() / 7; }, 0, "c=pi/7: deg 3" },
		{ "tan", [](const mpreal& t) { return mpfr::tan(t); }, [](int) { return mpfr::const_pi() / 5; }, 0, "c=pi/5: deg 4" },
		{ "sinh", [](const mpreal& t) { return mpfr::sinh(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "c=log(1+2^(1/3)): deg 3" },
		{ "cosh", [](const mpreal& t) { return mpfr::cosh(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "c=log(1+2^(1/3)): deg 3" },
		{ "tanh", [](const mpreal& t) { return mpfr::tanh(t); }, [](int) { return mpfr::log(1 + mpfr::cbrt(mpreal(2))); }, 0, "c=log(1+2^(1/3)): deg 3" },
		{ "arctan", [](const mpreal& t) { return mpfr::atan(t); }, [](int) { return mpfr::tan(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=tan(sqrt2/3): sqrt2/3" },
		{ "arcsinh", [](const mpreal& t) { return mpfr::asinh(t); }, [](int) { return mpfr::sinh(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=sinh(sqrt2/3): sqrt2/3" },
		{ "arccosh", [](const mpreal& t) { return mpfr::acosh(t); }, [](int) { return mpfr::cosh(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=cosh(sqrt2/3): sqrt2/3" },
		{ "arcsin(1/t)", [](const mpreal& t) { return mpfr::asin(1 / t); }, [](int) { return 1 / mpfr::sin(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=1/sin(sqrt2/3): sqrt2/3" },
		{ "arccos(1/t)", [](const mpreal& t) { return mpfr::acos(1 / t); }, [](int) { return 1 / mpfr::cos(mpfr::sqrt(mpreal(2)) / 3); }, 0, "x0=1/cos(sqrt2/3): sqrt2/3" },
		{ "log(t-1)", [](const mpreal& t) { return mpfr::log(t - 1); }, [](int) { return 1 + mpfr::exp(mpfr::sqrt(mpreal(2))); }, 0, "x0=1+exp(sqrt2): sqrt2" },
		{ "t^t", [](const mpreal& t) { return mpfr::pow(t, t); }, [](int) {
				// solve c^c = 2 by Newton on g(c) = c log c - log 2
				const mpreal l2 = mpfr::log(mpreal(2));
				mpreal c = 1.5;
				for (int i = 0; i < 200; ++i)
					c -= (c * mpfr::log(c) - l2) / (mpfr::log(c) + 1);
				return c; }, 0, "c^c=2: rational 2" },
		{ "Gamma", [](const mpreal& t) { return mpfr::gamma(t); }, [](int) { return mpreal(4); }, 0, "c=4: Gamma=6 (weak check)" },
		{ "exp(1/t)", [](const mpreal& t) { return mpfr::exp(1 / t); }, [](int) { return 1 / mpfr::log(mpreal(2)); }, 0, "c=1/log2: 2" },
	};
	std::vector<PolySearchTarget> main;
	for (const auto& fn : fns)
		for (int k : ks)
		{
			mpreal::set_default_prec(mpfr::digits2bits(50));
			if (mpfr::isnan(fn.f(mpfr::zeta(mpreal(3)) * mpfr::pow(mpfr::const_pi(), mpreal(k)))))
			{
				std::cout << "skip " << fn.name << "(zeta3*pi^" << k << "): outside the domain of the function" << std::endl;
				continue;
			}
			main.push_back({ fn.name + "(zeta3*pi^" + std::to_string(k) + ")", [fn, k](int d)
				{
					mpreal::set_default_prec(mpfr::digits2bits(d));
					return fn.f(mpfr::zeta(mpreal(3)) * mpfr::pow(mpfr::const_pi(), mpreal(k)));
				} });
		}
	applyStrideEnv(main);
	std::vector<PolySearchTarget> targets = main;
	const size_t nMain = targets.size();
	std::vector<int> checkMinDeg;
	for (const auto& fn : fns)
	{
		targets.push_back({ "[CHECK] " + fn.name + "(c / pi^2), c = x0 pi^2  (" + fn.checkNote + ")", [fn](int d)
			{
				mpreal::set_default_prec(mpfr::digits2bits(d));
				const mpreal pi2 = mpfr::const_pi() * mpfr::const_pi();
				const mpreal c = fn.checkC(d) * pi2;
				return fn.f(c / pi2);
			} });
		checkMinDeg.push_back(1);
	}
	// high-degree checks: cos(2 pi/65) has degree phi(65)/2 = 24; log(exp(2^(1/3)+3^(1/8))) = 2^(1/3)+3^(1/8) has degree 24.
	const int hd = 24;
	if (degree >= hd)
	{
		targets.push_back({ "[CHECK] cos(c / pi^2), x0 = 2 pi/65 (algebraic degree 24)", [](int d)
			{
				mpreal::set_default_prec(mpfr::digits2bits(d));
				const mpreal pi2 = mpfr::const_pi() * mpfr::const_pi();
				const mpreal c = (2 * mpfr::const_pi() / 65) * pi2;
				return mpfr::cos(c / pi2);
			} });
		checkMinDeg.push_back(hd);
		targets.push_back({ "[CHECK] log(c / pi^2), x0 = exp(2^(1/3)+3^(1/8)) (algebraic degree 24)", [](int d)
			{
				mpreal::set_default_prec(mpfr::digits2bits(d));
				const mpreal pi2 = mpfr::const_pi() * mpfr::const_pi();
				const mpreal c = mpfr::exp(mpfr::cbrt(mpreal(2)) + mpfr::root(mpreal(3), 8)) * pi2;
				return mpfr::log(c / pi2);
			} });
		checkMinDeg.push_back(hd);
	}
	std::cout << "NOTE: targets f(zeta3*pi^k), k in {";
	for (int k : ks) std::cout << " " << k;
	std::cout << " }. TARGET COUNT = " << nMain << " main + " << checkMinDeg.size() << " validation" << std::endl;
	const auto sm = reportPolynomialSearchParallel(targets, digits, degree, true, 10);
	reportHighDegSummary(targets, nMain, sm, checkMinDeg, digits, degree);
}

// Compositions y = f(g(zeta3 * pi^k)), k in ks, algebraicity only (degree == `degree` only).
void test_pslq_zeta3_composite_functions_highdeg(int digits, int degree, std::vector<int> ks)
{
	std::cout.precision(30);
	struct Fn
	{
		std::string name;
		std::function<mpreal(const mpreal&)> f;
		std::function<mpreal(const mpreal&)> inv;
		bool small;
	};
	const std::vector<Fn> fns = {
		{ "exp", [](const mpreal& t) { return mpfr::exp(t); }, [](const mpreal& v) { return mpfr::log(v); }, false },
		{ "log", [](const mpreal& t) { return mpfr::log(t); }, [](const mpreal& v) { return mpfr::exp(v); }, false },
		{ "sin", [](const mpreal& t) { return mpfr::sin(t); }, [](const mpreal& v) { return mpfr::asin(v); }, true },
		{ "cos", [](const mpreal& t) { return mpfr::cos(t); }, [](const mpreal& v) { return mpfr::acos(v); }, true },
		{ "tan", [](const mpreal& t) { return mpfr::tan(t); }, [](const mpreal& v) { return mpfr::atan(v); }, true },
		{ "sinh", [](const mpreal& t) { return mpfr::sinh(t); }, [](const mpreal& v) { return mpfr::asinh(v); }, false },
		{ "cosh", [](const mpreal& t) { return mpfr::cosh(t); }, [](const mpreal& v) { return mpfr::acosh(v); }, false },
		{ "tanh", [](const mpreal& t) { return mpfr::tanh(t); }, [](const mpreal& v) { return mpfr::atanh(v); }, true },
		{ "arctan", [](const mpreal& t) { return mpfr::atan(t); }, [](const mpreal& v) { return mpfr::tan(v); }, true },
		{ "arcsinh", [](const mpreal& t) { return mpfr::asinh(t); }, [](const mpreal& v) { return mpfr::sinh(v); }, false },
		{ "arcsin", [](const mpreal& t) { return mpfr::asin(t); }, [](const mpreal& v) { return mpfr::sin(v); }, true },
		{ "arccos", [](const mpreal& t) { return mpfr::acos(t); }, [](const mpreal& v) { return mpfr::cos(v); }, true },
	};
	auto isInversePair = [](const std::string& a, const std::string& b)
	{
		static const std::vector<std::pair<std::string, std::string>> inv = { { "exp", "log" }, { "sin", "arcsin" }, { "tan", "arctan" }, { "sinh", "arcsinh" }, { "cos", "arccos" } };
		for (auto& p : inv) if ((p.first == a && p.second == b) || (p.first == b && p.second == a)) return true;
		return false;
	};
	auto tOf = [](int k) { return mpfr::zeta(mpreal(3)) * mpfr::pow(mpfr::const_pi(), mpreal(k)); };
	std::vector<PolySearchTarget> main;
	int nDomain = 0, nRange = 0;
	for (const auto& f : fns)
		for (const auto& g : fns)
		{
			if (f.name == g.name || isInversePair(f.name, g.name)) continue;
			for (int k : ks)
			{
				mpreal::set_default_prec(mpfr::digits2bits(50));
				const mpreal y50 = f.f(g.f(tOf(k)));
				if (mpfr::isnan(y50) || mpfr::isinf(y50)) { ++nDomain; continue; }
				const mpreal ay = mpfr::abs(y50);
				if (ay < mpreal("1e-30") || ay > mpreal("1e30")) { ++nRange; continue; }
				main.push_back({ f.name + "(" + g.name + "(zeta3*pi^" + std::to_string(k) + "))", [f, g, k, tOf](int d)
					{
						mpreal::set_default_prec(mpfr::digits2bits(d));
						return f.f(g.f(tOf(k)));
					} });
			}
		}
	applyStrideEnv(main);
	std::vector<PolySearchTarget> targets = main;
	const size_t nMain = targets.size();
	std::vector<int> checkMinDeg;
	auto addChk = [&](const std::string& nm, std::function<mpreal()> v, int minDeg)
	{
		targets.push_back({ "[CHECK] " + nm, [v](int d) { mpreal::set_default_prec(mpfr::digits2bits(d)); return v(); } });
		checkMinDeg.push_back(minDeg);
	};
	addChk("cos(arcsin(1/3)) = sqrt8/3 (deg 2)", []() { return mpfr::cos(mpfr::asin(mpreal(1) / 3)); }, 2);
	addChk("sin(arctan(2)) = 2/sqrt5 (deg 2)", []() { return mpfr::sin(mpfr::atan(mpreal(2))); }, 2);
	addChk("cosh(arcsinh(sqrt3)) = 2 (deg 1)", []() { return mpfr::cosh(mpfr::asinh(mpfr::sqrt(mpreal(3)))); }, 1);
	addChk("cos(arctan(3^(1/4))) = 1/sqrt(1+sqrt3) (deg 4)", []() { return mpfr::cos(mpfr::atan(mpfr::root(mpreal(3), 4))); }, 4);
	addChk("tanh(log 5) = 12/13 (deg 1)", []() { return mpfr::tanh(mpfr::log(mpreal(5))); }, 1);
	addChk("sin(arccos(2^(1/3)/2)) (deg 6)", []() { return mpfr::sin(mpfr::acos(mpfr::cbrt(mpreal(2)) / 2)); }, 6);
	// high-degree: cosh(arccosh(a)) with a = 2^(1/4)+3^(1/4) (degree 16); cos(arctan(tan(2 pi/51))) = cos(2 pi/51) (degree phi(51)/2 = 16)
	const int hd = 16;
	if (degree >= hd)
	{
		addChk("cosh(arccosh(2^(1/4)+3^(1/4))) (deg 16)", []() { return mpfr::cosh(mpfr::acosh(mpfr::root(mpreal(2), 4) + mpfr::root(mpreal(3), 4))); }, hd);
		addChk("cos(arctan(tan(2 pi/51))) = cos(2 pi/51) (deg 16)", []() { return mpfr::cos(mpfr::atan(mpfr::tan(2 * mpfr::const_pi() / 51))); }, hd);
	}
	std::cout << "NOTE: targets y = f(g(zeta3*pi^k)), k in {";
	for (int k : ks) std::cout << " " << k;
	std::cout << " }; dropped " << nDomain << " (domain error), " << nRange << " (|y|<1e-30 or >1e30). TARGET COUNT = " << nMain
		<< " main + " << checkMinDeg.size() << " validation" << std::endl;
	const auto sm = reportPolynomialSearchParallel(targets, digits, degree, true, 10);
	reportHighDegSummary(targets, nMain, sm, checkMinDeg, digits, degree);
}


// Multi-level compositions y = f1(f2(...(fL(t)))), t = zeta3 * pi^k, L = minLevel..maxLevel; algebraicity only (degree == `degree` only).
// Pruning: adjacent equal or mutually inverse functions are dropped; domain/range errors, numerical duplicates (also against levels < minLevel)
// and targets that lose precision are dropped. Env ZHD_STRIDE (smoke test only) subsamples the final targets.
void test_pslq_zeta3_multilevel_composite(int digits, int degree, std::vector<int> ks, int minLevel, int maxLevel)
{
	std::cout.precision(30);
	enum Dom { None, Pos, Unit, Tan };
	struct Fn
	{
		std::string name;
		std::function<mpreal(const mpreal&)> f;
		Dom dom;
	};
	const std::vector<Fn> fns = {
		{ "exp", [](const mpreal& t) { return mpfr::exp(t); }, None },
		{ "log", [](const mpreal& t) { return mpfr::log(t); }, Pos },
		{ "sin", [](const mpreal& t) { return mpfr::sin(t); }, None },
		{ "cos", [](const mpreal& t) { return mpfr::cos(t); }, None },
		{ "tan", [](const mpreal& t) { return mpfr::tan(t); }, Tan },
		{ "sinh", [](const mpreal& t) { return mpfr::sinh(t); }, None },
		{ "cosh", [](const mpreal& t) { return mpfr::cosh(t); }, None },
		{ "tanh", [](const mpreal& t) { return mpfr::tanh(t); }, None },
		{ "arctan", [](const mpreal& t) { return mpfr::atan(t); }, None },
		{ "arcsinh", [](const mpreal& t) { return mpfr::asinh(t); }, None },
		{ "arcsin", [](const mpreal& t) { return mpfr::asin(t); }, Unit },
		{ "arccos", [](const mpreal& t) { return mpfr::acos(t); }, Unit },
	};
	const int NF = static_cast<int>(fns.size());
	auto idx = [&](const std::string& n) { for (int i = 0; i < NF; ++i) if (fns[i].name == n) return i; return -1; };
	std::vector<std::vector<char>> forbidden(NF, std::vector<char>(NF, 0)); // forbidden[a][b]: a adjacent to b is dropped
	for (int i = 0; i < NF; ++i) forbidden[i][i] = 1;
	const char* invs[][2] = { { "exp", "log" }, { "sin", "arcsin" }, { "cos", "arccos" }, { "tan", "arctan" }, { "sinh", "arcsinh" } };
	for (auto& p : invs) { forbidden[idx(p[0])][idx(p[1])] = 1; forbidden[idx(p[1])][idx(p[0])] = 1; }
	// algebraic pairs {outer, inner}: outer(inner(x)) is an algebraic function A(x), so A(h) is algebraic iff h is.
	// Outermost-pair case is exactly covered by the shorter composite; for inner pairs the reduction is NOT valid
	// (f1(A(t)) is a new target), so by default only the outermost pair is pruned. ZML_ALGPAIR_ALL=1 prunes every adjacent pair.
	const char* algs[][2] = { { "cos", "arcsin" }, { "tan", "arcsin" }, { "sin", "arccos" }, { "tan", "arccos" }, { "sin", "arctan" }, { "cos", "arctan" },
		{ "cosh", "arcsinh" }, { "tanh", "arcsinh" }, { "exp", "arcsinh" }, { "sinh", "log" }, { "cosh", "log" }, { "tanh", "log" } };
	std::vector<std::vector<char>> algOI(NF, std::vector<char>(NF, 0)); // algOI[outer][inner]
	for (auto& p : algs) algOI[idx(p[0])][idx(p[1])] = 1;
	const bool algAll = std::getenv("ZML_ALGPAIR_ALL") && std::atoi(std::getenv("ZML_ALGPAIR_ALL")) != 0;
	if (algAll)
		for (int o = 0; o < NF; ++o) for (int in = 0; in < NF; ++in) if (algOI[o][in]) forbidden[in][o] = 1;

	// apply(i, x): reason 0 = ok, 1 = domain/nan/inf, 2 = range (|y|>1e30 or <1e-30)
	const mpreal big("1e30"), small("1e-30"), poleEps("1e-10");
	auto apply = [&](int i, const mpreal& x, mpreal& out) -> int
	{
		const Fn& fn = fns[i];
		if (fn.dom == Pos && !(x > 0)) return 1;
		if (fn.dom == Unit && mpfr::abs(x) > 1) return 1;
		if (fn.dom == Tan && mpfr::abs(mpfr::cos(x)) < poleEps) return 1;
		out = fn.f(x);
		if (mpfr::isnan(out) || mpfr::isinf(out)) return 1;
		const mpreal a = mpfr::abs(out);
		if (a > big || a < small) return 2;
		return 0;
	};
	// zeta(3) and pi are cached per precision (thread-local): recomputing mpfr::zeta for each of ~35k targets
	// made the serial target-evaluation phase take hours.
	auto tOf = [](int k)
	{
		thread_local mp_prec_t cachedPrec = 0;
		thread_local mpreal z3, piV;
		const mp_prec_t p = mpreal::get_default_prec();
		if (cachedPrec != p) { z3 = mpfr::zeta(mpreal(3)); piV = mpfr::const_pi(); cachedPrec = p; }
		return z3 * mpfr::pow(piV, mpreal(k));
	};

	// number of rule-valid continuations (outward) of length n after function `last`
	std::vector<std::vector<double>> comp(maxLevel + 1, std::vector<double>(NF, 1.0));
	for (int n = 1; n <= maxLevel; ++n)
		for (int a = 0; a < NF; ++a)
		{
			double s = 0;
			for (int b = 0; b < NF; ++b) if (!forbidden[a][b]) s += comp[n - 1][b];
			comp[n][a] = s;
		}

	struct Entry
	{
		mpreal v;
		int level, k;
		std::vector<int> app; // application order: app[0] innermost
		int fate = 0;         // 0 keep; 3 dup of lower level; 4 dup within run; 5 precision loss; 6 algebraic outermost pair
		bool ref = false;     // reference entry (not a target)
	};
	std::vector<Entry> all;
	const int nk = static_cast<int>(ks.size());
	std::vector<double> algTop(maxLevel + 1, 0), dropDomain(maxLevel + 1, 0), dropRange(maxLevel + 1, 0), valid(maxLevel + 1, 0), total(maxLevel + 1, 0);
	mpreal::set_default_prec(256);
	for (int L = 1; L <= maxLevel; ++L)
	{
		double tot = nk;
		for (int i = 0; i < L; ++i) tot *= NF;
		total[L] = tot;
		double v = 0;
		for (int a = 0; a < NF; ++a) v += comp[L - 1][a];
		valid[L] = v * nk;
	}
	for (int ki = 0; ki < nk; ++ki)
	{
		const mpreal t = tOf(ks[ki]);
		std::vector<int> app;
		std::function<void(const mpreal&)> dfs = [&](const mpreal& x)
		{
			const int depth = static_cast<int>(app.size());
			if (depth >= maxLevel) return;
			for (int i = 0; i < NF; ++i)
			{
				if (depth > 0 && forbidden[app.back()][i]) continue;
				mpreal y;
				const int r = apply(i, x, y);
				if (r)
				{
					for (int l = depth + 1; l <= maxLevel; ++l)
						(r == 1 ? dropDomain : dropRange)[l] += comp[l - depth - 1][i];
					continue;
				}
				app.push_back(i);
				all.push_back({ y, depth + 1, ks[ki], app, 0, false });
				if (!algAll && depth >= 1 && depth + 1 >= minLevel && algOI[i][app[app.size() - 2]]) { all.back().fate = 6; ++algTop[depth + 1]; }
				dfs(y);
				app.pop_back();
			}
		};
		dfs(t);
	}
	std::cout << "enumerated " << all.size() << " raw composites" << std::endl;
	// generation order: level, then outer-first lexicographic sequence, then k
	std::vector<size_t> ord(all.size());
	std::iota(ord.begin(), ord.end(), 0);
	auto outerFirst = [&](const Entry& e) { std::vector<int> s(e.app.rbegin(), e.app.rend()); return s; };
	std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b)
		{
			if (all[a].level != all[b].level) return all[a].level < all[b].level;
			const auto sa = outerFirst(all[a]), sb = outerFirst(all[b]);
			if (sa != sb) return sa < sb;
			return all[a].k < all[b].k;
		});
	std::vector<size_t> rank(all.size());
	for (size_t r = 0; r < ord.size(); ++r) rank[ord[r]] = r;

	// numerical dedup, per k (same t): key = max(|v|, 1/|v|) so that y, -y, 1/y, -1/y coincide.
	// Reference entries: t itself (level 0) and every composite of level < minLevel (already tested).
	const mpreal tol("1e-40");
	std::vector<double> dupLower(maxLevel + 1, 0), dupRun(maxLevel + 1, 0), precLoss(maxLevel + 1, 0);
	std::vector<double> dupLowerInv(maxLevel + 1, 0), dupRunInv(maxLevel + 1, 0);
	for (int k : ks)
	{
		all.push_back({ tOf(k), 0, k, {}, 0, true });
		all.back().ref = true;
		rank.push_back(0);
	}
	auto keyOf = [&](const Entry& e) { const mpreal a = mpfr::abs(e.v); return a < 1 ? 1 / a : a; };
	for (int k : ks)
	{
		std::vector<size_t> ids;
		for (size_t i = 0; i < all.size(); ++i) if (all[i].k == k && all[i].fate != 6) ids.push_back(i);
		std::vector<mpreal> key(all.size());
		for (size_t i : ids) key[i] = keyOf(all[i]);
		std::sort(ids.begin(), ids.end(), [&](size_t a, size_t b) { return key[a] < key[b]; });
		for (size_t s = 0; s < ids.size();)
		{
			size_t e = s + 1;
			while (e < ids.size() && key[ids[e]] - key[ids[e - 1]] <= tol * key[ids[e]]) ++e;
			if (e - s > 1)
			{
				bool hasRef = false;
				size_t best = ids[s];
				bool bestSet = false;
				for (size_t j = s; j < e; ++j)
				{
					const Entry& en = all[ids[j]];
					if (en.ref || en.level < minLevel) hasRef = true;
					else if (!bestSet || rank[ids[j]] < rank[best]) { best = ids[j]; bestSet = true; }
				}
				for (size_t j = s; j < e; ++j)
				{
					Entry& en = all[ids[j]];
					if (en.ref || en.level < minLevel) continue;
					if (!hasRef && ids[j] == best) continue;
					// direct (+-y) match or via 1/y? compare |v| with some reference/best entry
					bool direct = false;
					for (size_t q = s; q < e; ++q)
					{
						if (q == j) continue;
						const Entry& o = all[ids[q]];
						const bool oRef = o.ref || o.level < minLevel;
						if (hasRef ? !oRef : ids[q] != best) continue;
						if (mpfr::abs(mpfr::abs(o.v) - mpfr::abs(en.v)) <= tol * mpfr::abs(en.v)) { direct = true; break; }
					}
					if (hasRef) { en.fate = 3; ++(direct ? dupLower : dupLowerInv)[en.level]; }
					else { en.fate = 4; ++(direct ? dupRun : dupRunInv)[en.level]; }
				}
			}
			s = e;
		}
	}
	std::cout << "dedup done" << std::endl;
	// precision check on survivors: recompute at 2 precisions, need >= `digits` agreeing digits
	const mpfr_prec_t bitsD = mpfr::digits2bits(digits);
	std::map<std::pair<int, mpfr_prec_t>, mpreal> tCache;
	for (int k : ks)
		for (mpfr_prec_t bits : { bitsD + 50, bitsD + 150 })
		{
			mpreal::set_default_prec(bits);
			tCache[{ k, bits }] = tOf(k);
		}
	auto evalSeq = [&](const std::vector<int>& app, int k, mpfr_prec_t bits, mpreal& out) -> bool
	{
		mpreal::set_default_prec(bits);
		mpreal x = tCache.at({ k, bits });
		for (int i : app)
		{
			mpreal y;
			if (apply(i, x, y)) return false;
			x = y;
		}
		out = x;
		return true;
	};
	{
		std::vector<size_t> todo;
		for (size_t i = 0; i < all.size(); ++i)
			if (!all[i].ref && all[i].level >= minLevel && !all[i].fate) todo.push_back(i);
		std::atomic<size_t> nxt{ 0 };
		std::vector<std::thread> pool;
		for (int th = 0; th < 10; ++th)
			pool.emplace_back([&]()
				{
					for (size_t q; (q = nxt++) < todo.size();)
					{
						Entry& en = all[todo[q]];
						mpreal a, b;
						bool ok = evalSeq(en.app, en.k, bitsD + 50, a) && evalSeq(en.app, en.k, bitsD + 150, b);
						if (ok)
						{
							mpreal::set_default_prec(bitsD + 150);
							ok = mpfr::abs(a - b) < mpfr::abs(b) * mpfr::pow(mpreal(10), -digits);
						}
						if (!ok) en.fate = 5;
					}
				});
		for (auto& th : pool) th.join();
		for (size_t i : todo) if (all[i].fate == 5) ++precLoss[all[i].level];
	}
	auto seqName = [&](const Entry& e)
	{
		std::string s;
		for (size_t i = e.app.size(); i-- > 0;) s += fns[e.app[i]].name + "(";
		s += "zeta3*pi^" + std::to_string(e.k);
		s += std::string(e.app.size(), ')');
		return s;
	};
	std::cout << "NOTE: multilevel targets, k in {";
	for (int k : ks) std::cout << " " << k;
	std::cout << " }, levels " << minLevel << ".." << maxLevel << ", " << digits << " digits, degree " << degree << std::endl;
	std::vector<PolySearchTarget> main;
	for (size_t r = 0; r < ord.size(); ++r)
	{
		const Entry& en = all[ord[r]];
		if (en.ref || en.level < minLevel || en.fate) continue;
		auto cache = std::make_shared<std::map<mpfr_prec_t, mpreal>>(); // t at given bits (xAtPrec is called from one thread only)
		const std::vector<int> app = en.app;
		const int k = en.k;
		main.push_back({ seqName(en), [=](int d)
			{
				const mpfr_prec_t bits = mpfr::digits2bits(d) + 50;
				mpreal::set_default_prec(bits);
				auto it = cache->find(bits);
				if (it == cache->end()) it = cache->emplace(bits, tOf(k)).first;
				mpreal x = it->second;
				for (int i : app)
				{
					x = fns[i].f(x);
				}
				mpreal::set_default_prec(mpfr::digits2bits(d));
				mpreal r(0, mpfr::digits2bits(d));
				r = x;
				return r;
			} });
	}
	// funcs captured by reference above would dangle; fns is local to this function, which outlives reportPolynomialSearchParallel. OK.
	for (int L = minLevel; L <= maxLevel; ++L)
	{
		double kept = 0;
		for (const auto& en : all) if (!en.ref && en.level == L && !en.fate) ++kept;
		std::cout << "LEVEL " << L << ": all sequences " << total[L] << ", rule-pruned (equal/inverse neighbours) " << total[L] - valid[L]
			<< ", algebraic outermost pair " << algTop[L] << ", domain errors " << dropDomain[L] << ", range (|y|>1e30 or <1e-30) " << dropRange[L]
			<< ", dup of lower-level composite/t (+-y) " << dupLower[L] << " (via 1/y) " << dupLowerInv[L]
			<< ", dup within run (+-y) " << dupRun[L] << " (via 1/y) " << dupRunInv[L]
			<< ", precision loss " << precLoss[L] << ", FINAL " << kept << std::endl;
	}
	std::cout << "TOTAL TARGETS BEFORE RANGE/SKIP: " << main.size() << std::endl;
	if (const char* r = std::getenv("ZML_RANGE")) // a:b -> main targets with index in [a, b)
	{
		size_t lo = 0, hi = main.size();
		if (std::sscanf(r, "%zu:%zu", &lo, &hi) == 2)
		{
			hi = std::min(hi, main.size());
			lo = std::min(lo, hi);
			main = std::vector<PolySearchTarget>(main.begin() + lo, main.begin() + hi);
			std::cout << "RESUME: ZML_RANGE=" << lo << ":" << hi << std::endl;
		}
	}
	if (const char* f = std::getenv("ZML_SKIP_DONE")) // previous output: skip targets whose job line shows a non-Found status at this degree
	{
		std::ifstream in(f);
		std::set<std::string> done;
		std::string line;
		const std::string tag = ", degree " + std::to_string(degree) + ": ";
		while (std::getline(in, line))
		{
			const size_t b = line.find("] ");
			const size_t e = line.rfind(tag);
			if (line.find("  [") != 0 || b == std::string::npos || e == std::string::npos || e < b) continue;
			const std::string st = line.substr(e + tag.size());
			if (st.rfind("Found", 0) == 0 || st.rfind("Precision", 0) == 0) continue; // re-run (need verification / unreliable)
			done.insert(line.substr(b + 2, e - b - 2));
		}
		std::vector<PolySearchTarget> rest;
		for (auto& t : main) if (!done.count(t.name)) rest.push_back(t);
		std::cout << "RESUME: ZML_SKIP_DONE=" << f << ": skipped " << main.size() - rest.size() << " done targets, " << rest.size() << " remain" << std::endl;
		main = rest;
	}
	applyStrideEnv(main);
	std::vector<PolySearchTarget> targets = main;
	const size_t nMain = targets.size();
	std::vector<int> checkMinDeg;
	auto addChk = [&](const std::string& nm, std::function<mpreal()> v, int minDeg)
	{
		targets.push_back({ "[CHECK] " + nm, [v](int d) { mpreal::set_default_prec(mpfr::digits2bits(d)); return v(); } });
		checkMinDeg.push_back(minDeg);
	};
	addChk("cos(arctan(sinh(arcsinh(2)))) = 1/sqrt5 (deg 2)", []() { return mpfr::cos(mpfr::atan(mpfr::sinh(mpfr::asinh(mpreal(2))))); }, 2);
	addChk("tanh(log(cosh(arccosh(3)))) = 4/5 (deg 1)", []() { return mpfr::tanh(mpfr::log(mpfr::cosh(mpfr::acosh(mpreal(3))))); }, 1);
	addChk("cosh(arcsinh(sin(arcsin(2^(-1/3))))) = sqrt(1+2^(-2/3)) (deg 6)", []() { return mpfr::cosh(mpfr::asinh(mpfr::sin(mpfr::asin(mpfr::cbrt(mpreal(1) / 2))))); }, 6);
	addChk("sin(arctan(sinh(arcsinh(3^(1/8))))) = c/sqrt(1+c^2), c=3^(1/8) (deg 8)", []() { return mpfr::sin(mpfr::atan(mpfr::sinh(mpfr::asinh(mpfr::root(mpreal(3), 8))))); }, 8);
	{
		const int D = degree;
		addChk("tanh(log(cosh(arccosh(2^(1/D)+1)))) = (a^2-1)/(a^2+1), D=" + std::to_string(D) + " (deg D)",
			[D]() { return mpfr::tanh(mpfr::log(mpfr::cosh(mpfr::acosh(mpfr::root(mpreal(2), D) + 1)))); }, D);
		addChk("tan(arctan(cos(arccos(2^(-1/D))))) = 2^(-1/D), D=" + std::to_string(D) + " (deg D)",
			[D]() { return mpfr::tan(mpfr::atan(mpfr::cos(mpfr::acos(mpreal(1) / mpfr::root(mpreal(2), D))))); }, D);
	}
	std::cout << "TARGET COUNT = " << nMain << " main + " << checkMinDeg.size() << " validation" << std::endl;
	const auto sm = reportPolynomialSearchParallel(targets, digits, degree, true, 10);
	reportHighDegSummary(targets, nMain, sm, checkMinDeg, digits, degree);
}

// y = f(zeta3 * (a x)^(b w)) algebraic (= alpha), a algebraic, b rational, x and w transcendental.  u = g(alpha) > 0,  g = f^-1.
//   log zeta3 = log u - b w (log a + log x).
// Forms (what is absorbed into the PSLQ basis, so that no alpha grid is needed):
//   GEN (sin, cos, tan, sinh, ..., alpha inside a transcendental g; alpha from a grid):
//        (L = log u - log zeta3, w*{log2,log3,log5,log7,log(1+sqrt2),log phi,1,log pi,log gamma,log log2,log log3} [, 1])         n = 12 [13]
//   identity (u = alpha):  (log zeta3, log2,log3,log5,log7,log(1+sqrt2),log phi, w*{...11})                                       n = 18 [19]
//   log (u = exp(alpha)):  (log zeta3, 1,sqrt2,sqrt3,sqrt5,sqrt6, w*{...11})                                                      n = 17
//   exp (u = log(alpha)):  (log zeta3, log2,log3,log5,log7, loglog2,loglog3,loglog5,loglog7, w*{...11} [, 1])                    n = 20 [21]
// When w = 1, the w*X columns that duplicate plain columns are dropped.
// Env vars (smoke tests): ZETA3_TE_THREADS (default 10), ZETA3_TE_JOBS (evenly spaced sample of that many production jobs),
//   ZETA3_TE_PMAX / ZETA3_TE_QMAX (rational grid, default 60 / 30), ZETA3_TE_QUAD (0/1, default 1), ZETA3_TE_COUNT_ONLY (1: print job counts only),
//   ZETA3_TE_SKIP_DONE=<old output path> (resume: skip production jobs whose old line has status none; others and CHECK jobs rerun).
// variant: 0 = without constant column, 1 = with constant column only, 2 = both.
static void zeta3TranscendentalExponentImpl(int digits, int variant, bool absorbedMode, int maxPrime, int maxPrimeLogLog)
{
	std::cout.precision(30);
	const auto tStart = std::chrono::steady_clock::now();
	auto envInt = [](const char* n, int d) { const char* e = std::getenv(n); return e ? std::atoi(e) : d; };
	const int pMax = envInt("ZETA3_TE_PMAX", 60), qMax = envInt("ZETA3_TE_QMAX", 30), useQuad = envInt("ZETA3_TE_QUAD", 1);
	const int nThreads = std::max(1, envInt("ZETA3_TE_THREADS", 10)), sampleJobs = envInt("ZETA3_TE_JOBS", 0), countOnly = envInt("ZETA3_TE_COUNT_ONLY", 0);
	auto nanv = []() { mpreal r; r.setNan(); return r; };
	const auto savedPrec = mpreal::get_default_prec();
	mpreal::set_default_prec(mpfr::digits2bits(50));

	// ---- functions f with the branches of g = f^-1 (u values; NaN when outside the domain) ----
	typedef std::function<std::vector<mpreal>(const mpreal&)> InvFn;
	struct TFn { std::string name; int form; std::vector<std::string> bl; InvFn inv; };
	enum { GEN = 0, ID = 1, LOG = 2, EXP = 3 };
	auto newton = [](const std::function<mpreal(const mpreal&)>& g, const std::function<mpreal(const mpreal&)>& dg, mpreal x)
	{
		const mpreal eps = mpfr::pow(mpreal(2), -(long)mpreal::get_default_prec() + 20);
		for (int i = 0; i < 300; ++i)
		{
			const mpreal dx = g(x) / dg(x);
			x -= dx;
			if (mpfr::abs(dx) < eps) break;
		}
		return x;
	};
	auto single = [](const std::function<mpreal(const mpreal&)>& g) { return InvFn([g](const mpreal& a) { return std::vector<mpreal>{ g(a) }; }); };
	const std::vector<TFn> fns = {
		{ "identity", ID, { "u=alpha" }, single([](const mpreal& a) { return a; }) },
		{ "log", LOG, { "exp(alpha)" }, single([](const mpreal& a) { return mpfr::exp(a); }) },
		{ "exp", EXP, { "log(alpha)" }, single([](const mpreal& a) { return mpfr::log(a); }) },
		{ "sin", GEN, { "u0", "u0+2pi", "pi-u0", "3pi-u0" }, [](const mpreal& a) {
				const mpreal u0 = mpfr::asin(a), p = mpfr::const_pi();
				return std::vector<mpreal>{ u0, u0 + 2 * p, p - u0, 3 * p - u0 }; } },
		{ "cos", GEN, { "u0", "u0+2pi", "2pi-u0", "4pi-u0" }, [](const mpreal& a) {
				const mpreal u0 = mpfr::acos(a), p = mpfr::const_pi();
				return std::vector<mpreal>{ u0, u0 + 2 * p, 2 * p - u0, 4 * p - u0 }; } },
		{ "tan", GEN, { "u0", "u0+pi", "u0+2pi" }, [](const mpreal& a) {
				const mpreal u0 = mpfr::atan(a), p = mpfr::const_pi();
				return std::vector<mpreal>{ u0, u0 + p, u0 + 2 * p }; } },
		{ "sinh", GEN, { "asinh(alpha)" }, single([](const mpreal& a) { return mpfr::asinh(a); }) },
		{ "cosh", GEN, { "acosh(alpha)" }, single([](const mpreal& a) { return mpfr::acosh(a); }) },
		{ "tanh", GEN, { "atanh(alpha)" }, single([](const mpreal& a) { return mpfr::atanh(a); }) },
		{ "arctan", GEN, { "tan(alpha)" }, single([nanv](const mpreal& a) { return (a > 0 && a < mpfr::const_pi() / 2) ? mpfr::tan(a) : nanv(); }) },
		{ "arcsinh", GEN, { "sinh(alpha)" }, single([](const mpreal& a) { return mpfr::sinh(a); }) },
		{ "arccosh", GEN, { "cosh(alpha)" }, single([](const mpreal& a) { return mpfr::cosh(a); }) },
		{ "arcsin(1/t)", GEN, { "1/sin(alpha)" }, single([nanv](const mpreal& a) { return (a > 0 && a <= mpfr::const_pi() / 2) ? 1 / mpfr::sin(a) : nanv(); }) },
		{ "arccos(1/t)", GEN, { "1/cos(alpha)" }, single([nanv](const mpreal& a) { return (a > 0 && a < mpfr::const_pi() / 2) ? 1 / mpfr::cos(a) : nanv(); }) },
		{ "log(t-1)", GEN, { "1+exp(alpha)" }, single([](const mpreal& a) { return 1 + mpfr::exp(a); }) },
		{ "t^t", GEN, { "u^u=alpha, u>1/e" }, single([nanv, newton](const mpreal& a) {
				if (!(a >= mpfr::exp(-1 / mpfr::exp(mpreal(1))))) return nanv();
				if (a == 1) return mpreal(1);
				const mpreal la = mpfr::log(a);
				return newton([&](const mpreal& c) { return c * mpfr::log(c) - la; }, [](const mpreal& c) { return mpfr::log(c) + 1; }, a > 1 ? a + 1 : mpreal(1)); }) },
		{ "Gamma", GEN, { "Gamma(u)=alpha, u>1.4616" }, single([nanv, newton](const mpreal& a) {
				if (!(a >= mpfr::gamma(mpreal("1.46163214496836234126265954232572132846819620400644")))) return nanv();
				const mpreal la = mpfr::log(a);
				return newton([&](const mpreal& c) { return mpfr::lngamma(c) - la; }, [](const mpreal& c) { return mpfr::digamma(c); }, a + 2); }) },
		{ "exp(1/t)", GEN, { "1/log(alpha)" }, single([](const mpreal& a) { return 1 / mpfr::log(a); }) },
	};
	auto findF = [&](const std::string& n) { for (size_t i = 0; i < fns.size(); ++i) if (fns[i].name == n) return (int)i; return -1; };

	// ---- alpha grid (GEN functions only) ----
	struct Alpha { std::string name; std::function<mpreal()> val; };
	std::vector<Alpha> alphas;
	for (int p = 1; !absorbedMode && p <= pMax; ++p)
		for (int q = 1; q <= qMax; ++q)
			if (std::gcd(p, q) == 1)
			{
				alphas.push_back({ std::to_string(p) + "/" + std::to_string(q), [p, q]() { return mpreal(p) / q; } });
				alphas.push_back({ "-" + std::to_string(p) + "/" + std::to_string(q), [p, q]() { return -mpreal(p) / q; } });
			}
	if (useQuad && !absorbedMode)
	{
		for (int d : { 2, 3, 5 })
			for (int u = -3; u <= 3; ++u)
				for (int v = -3; v <= 3; ++v)
					for (int c = 1; c <= 4; ++c)
					{
						if (v == 0 || std::gcd(std::gcd(std::abs(u), std::abs(v)), c) != 1) continue;
						alphas.push_back({ "(" + std::to_string(u) + "+" + std::to_string(v) + "*sqrt" + std::to_string(d) + ")/" + std::to_string(c),
							[=]() { return (mpreal(u) + v * mpfr::sqrt(mpreal(d))) / c; } });
					}
		for (int b : { 2, 3 })
			for (int inv = 0; inv < 2; ++inv)
				for (int sg = 0; sg < 2; ++sg)
					alphas.push_back({ std::string(sg ? "-" : "") + (inv ? "1/" : "") + std::to_string(b) + "^(1/3)",
						[=]() { mpreal r = mpfr::cbrt(mpreal(b)); if (inv) r = 1 / r; return sg ? -r : r; } });
	}

	// ---- w ----
	struct WSpec { std::string name; std::function<mpreal()> val; };
	const std::vector<WSpec> ws = {
		{ "1", []() { return mpreal(1); } }, { "pi", []() { return mpfr::const_pi(); } }, { "e", []() { return mpfr::exp(mpreal(1)); } },
		{ "gamma", []() { return mpfr::const_euler(); } }, { "log2", []() { return mpfr::log(mpreal(2)); } },
		{ "log3", []() { return mpfr::log(mpreal(3)); } }, { "pi^2", []() { return mpfr::const_pi() * mpfr::const_pi(); } },
		{ "sqrt2*pi", []() { return mpfr::sqrt(mpreal(2)) * mpfr::const_pi(); } },
	};
	auto findW = [&](const std::string& n) { for (size_t i = 0; i < ws.size(); ++i) if (ws[i].name == n) return (int)i; return -1; };

	// ---- columns ----
	std::vector<int> primes;
	for (int q = 2; q <= std::max(maxPrime, maxPrimeLogLog); ++q)
	{
		bool pr = true;
		for (int d = 2; d * d <= q; ++d) if (q % d == 0) pr = false;
		if (pr) primes.push_back(q);
	}
	std::map<std::string, std::function<mpreal()>> base = {
		{ "1", []() { return mpreal(1); } },
		{ "log(1+sqrt2)", []() { return mpfr::log(1 + mpfr::sqrt(mpreal(2))); } },
		{ "log(phi)", []() { return mpfr::log((1 + mpfr::sqrt(mpreal(5))) / 2); } },
		{ "log(pi)", []() { return mpfr::log(mpfr::const_pi()); } }, { "log(gamma)", []() { return mpfr::log(mpfr::const_euler()); } },
		{ "sqrt2", []() { return mpfr::sqrt(mpreal(2)); } }, { "sqrt3", []() { return mpfr::sqrt(mpreal(3)); } },
		{ "sqrt5", []() { return mpfr::sqrt(mpreal(5)); } }, { "sqrt6", []() { return mpfr::sqrt(mpreal(6)); } },
	};
	for (int q : primes)
	{
		base["log" + std::to_string(q)] = [q]() { return mpfr::log(mpreal(q)); };
		base["log(log" + std::to_string(q) + ")"] = [q]() { return mpfr::log(mpfr::log(mpreal(q))); };
	}
	typedef std::vector<std::pair<std::string, bool>> ColSpec; // (base name, multiplied by w)
	auto colSpec = [&](int form, const std::string& wname, bool wc)
	{
		std::vector<std::string> plain, logP, llP, wf;
		for (int q : primes) { if (q <= maxPrime) logP.push_back("log" + std::to_string(q)); if (q <= maxPrimeLogLog) llP.push_back("log(log" + std::to_string(q) + ")"); }
		if (form == ID) { plain = logP; plain.push_back("log(1+sqrt2)"); plain.push_back("log(phi)"); }
		if (form == LOG) plain = { "1", "sqrt2", "sqrt3", "sqrt5", "sqrt6" };
		if (form == EXP) { plain = logP; plain.insert(plain.end(), llP.begin(), llP.end()); }
		wf = logP;
		for (const char* n : { "log(1+sqrt2)", "log(phi)", "1", "log(pi)", "log(gamma)", "log(log2)", "log(log3)" }) wf.push_back(n);
		auto has = [&](const std::string& n) { return std::find(plain.begin(), plain.end(), n) != plain.end(); };
		ColSpec cs;
		for (const auto& p : plain) cs.push_back({ p, false });
		for (const auto& n : wf)
		{
			// drop w*X when it numerically equals a plain column
			if (wname == "1" && has(n)) continue;
			if (n == "1" && wname == "log2" && has("log2")) continue;
			if (n == "1" && wname == "log3" && has("log3")) continue;
			cs.push_back({ n, true });
		}
		if (wc) cs.push_back({ "1", false });
		return cs;
	};

	// ---- jobs ----
	struct Job
	{
		int fi = 0, ai = -1, br = 0, wi = 0, wc = 0, ci = -1;
		LinearSearchOutcome oc;
		double secs = 0;
		std::string kind, name;
		bool checkPass = false;
	};
	auto probeBranches = [&](const TFn& fn, const mpreal& a)
	{
		std::vector<int> keep;
		const std::vector<mpreal> us = fn.inv(a);
		for (size_t b = 0; b < us.size(); ++b)
		{
			const mpreal& u = us[b];
			if (!mpfr::isfinite(u) || !(u > mpreal("1e-8"))) continue;
			bool dup = false;
			for (int k : keep) if (mpfr::abs(us[k] - u) < mpreal("1e-30")) dup = true;
			if (!dup) keep.push_back((int)b);
		}
		return keep;
	};
	// validation cases: c = g(alpha0) / (a0 x0)^(b0 w0); expected relation given by labels
	struct ChkSpec
	{
		std::string f; Alpha al; int br; std::string w, desc; std::function<mpreal()> a0, x0; long bp, bq;
		std::vector<std::pair<std::string, long>> exp;
		int needPrime = 0;
	};
	const std::vector<ChkSpec> specs = {
		{ "tan", { "1/2", []() { return mpreal(1) / 2; } }, 0, "e", "a0=3, x0=pi, b0=2/5",
			[]() { return mpreal(3); }, []() { return mpfr::const_pi(); }, 2, 5, { { "L", 5 }, { "w*log3", -2 }, { "w*log(pi)", -2 } } },
		{ "sin", { "1/3", []() { return mpreal(1) / 3; } }, 2, "gamma", "reflected branch pi-asin(1/3), a0=5/7, x0=log3, b0=1",
			[]() { return mpreal(5) / 7; }, []() { return mpfr::log(mpreal(3)); }, 1, 1, { { "L", 1 }, { "w*log5", -1 }, { "w*log7", 1 }, { "w*log(log3)", -1 } } },
		{ "exp", { "3", []() { return mpreal(3); } }, 0, "log2", "alpha0=3 (on grid), a0=sqrt2, x0=gamma, b0=-3/4",
			[]() { return mpfr::sqrt(mpreal(2)); }, []() { return mpfr::const_euler(); }, -3, 4, { { "logz", 8 }, { "log(log3)", -8 }, { "w*log2", -3 }, { "w*log(gamma)", -6 } } },
		{ "identity", { "1/sqrt2", []() { return 1 / mpfr::sqrt(mpreal(2)); } }, 0, "pi", "alpha0=1/sqrt2 (on grid), a0=phi, x0=e, b0=1/3",
			[]() { return (1 + mpfr::sqrt(mpreal(5))) / 2; }, []() { return mpfr::exp(mpreal(1)); }, 1, 3, { { "logz", 6 }, { "log2", 3 }, { "w*log(phi)", 2 }, { "w*1", 2 } } },
		{ "identity", { "2^(3/7)*3^(-1/5)", []() { return mpfr::pow(mpreal(2), mpreal(3) / 7) * mpfr::pow(mpreal(3), -mpreal(1) / 5); } }, 0, "e",
			"alpha0=2^(3/7)3^(-1/5) (off grid), a0=1+sqrt2, x0=gamma, b0=3/4",
			[]() { return 1 + mpfr::sqrt(mpreal(2)); }, []() { return mpfr::const_euler(); }, 3, 4,
			{ { "logz", 140 }, { "log2", -60 }, { "log3", 28 }, { "w*log(1+sqrt2)", 105 }, { "w*log(gamma)", 105 } } },
		{ "log", { "1/3+sqrt5/7", []() { return mpreal(1) / 3 + mpfr::sqrt(mpreal(5)) / 7; } }, 0, "log3",
			"alpha0=1/3+sqrt5/7 (off grid), a0=2, x0=pi, b0=-2/3",
			[]() { return mpreal(2); }, []() { return mpfr::const_pi(); }, -2, 3, { { "logz", 21 }, { "1", -7 }, { "sqrt5", -3 }, { "w*log2", -14 }, { "w*log(pi)", -14 } } },
		{ "exp", { "5^(2/9)", []() { return mpfr::pow(mpreal(5), mpreal(2) / 9); } }, 0, "gamma",
			"alpha0=5^(2/9) (off grid), a0=3, x0=e, b0=1/2",
			[]() { return mpreal(3); }, []() { return mpfr::exp(mpreal(1)); }, 1, 2,
			{ { "logz", 2 }, { "log2", -2 }, { "log3", 4 }, { "log(log5)", -2 }, { "w*log3", 1 }, { "w*1", 1 } } },
		{ "identity", { "29^(1/3)/2", []() { return mpfr::cbrt(mpreal(29)) / 2; } }, 0, "pi",
			"alpha0=29^(1/3)/2 (prime 29, off grid), a0=phi, x0=e, b0=2/3",
			[]() { return (1 + mpfr::sqrt(mpreal(5))) / 2; }, []() { return mpfr::exp(mpreal(1)); }, 2, 3,
			{ { "logz", 3 }, { "log29", -1 }, { "log2", 3 }, { "w*log(phi)", 2 }, { "w*1", 2 } }, 29 },
		{ "identity", { "3^(1/4)*13^(-1/2)", []() { return mpfr::pow(mpreal(3), mpreal(1) / 4) * mpfr::pow(mpreal(13), -mpreal(1) / 2); } }, 0, "e",
			"alpha0=3^(1/4)13^(-1/2), a0=23 (prime 23 in w*log), x0=gamma, b0=-1/2",
			[]() { return mpreal(23); }, []() { return mpfr::const_euler(); }, -1, 2,
			{ { "logz", 4 }, { "log3", -1 }, { "log13", 2 }, { "w*log23", -2 }, { "w*log(gamma)", -2 } }, 23 },
	};

	std::vector<Job> checks, prod;
	for (size_t ci = 0; ci < specs.size(); ++ci)
	{
		const int fi = findF(specs[ci].f);
		if ((fns[fi].form == GEN) == absorbedMode || specs[ci].needPrime > maxPrime) continue;
		for (int wc = 0; wc < 2; ++wc)
		{
			if (wc && fns[fi].form == LOG) continue; // plain 1 is already in the basis
			Job j; j.fi = fi; j.br = specs[ci].br; j.wi = findW(specs[ci].w); j.wc = wc; j.ci = (int)ci;
			checks.push_back(j);
		}
	}
	std::vector<long> jobsPerF(fns.size(), 0);
	for (size_t fi = 0; fi < fns.size(); ++fi)
	{
		if ((fns[fi].form == GEN) == absorbedMode) continue;
		std::vector<std::pair<int, int>> abList; // (alpha index, branch)
		if (fns[fi].form == GEN)
		{
			for (size_t ai = 0; ai < alphas.size(); ++ai)
				for (int b : probeBranches(fns[fi], alphas[ai].val()))
					abList.push_back({ (int)ai, b });
		}
		else abList.push_back({ -1, 0 });
		for (const auto& ab : abList)
			for (size_t wi = 0; wi < ws.size(); ++wi)
				for (int wc = 0; wc < 2; ++wc)
				{
					if (wc == 0 && variant == 1) continue;
					if (wc == 1 && variant == 0) continue;
					if (wc == 1 && (ws[wi].name == "1" || fns[fi].form == LOG)) continue;
					Job j; j.fi = (int)fi; j.ai = ab.first; j.br = ab.second; j.wi = (int)wi; j.wc = wc;
					prod.push_back(j);
					++jobsPerF[fi];
				}
	}
	const size_t totalProd = prod.size();
	std::cout << "NOTE: y = f(zeta3 (a x)^(b w)) algebraic. log zeta3 = log u - b w (log a + log x), u = g(alpha) > 0. digits = " << digits
		<< ", variant = " << variant << ", threads = " << nThreads << ", rational grid p<=" << pMax << " q<=" << qMax << ", quad/cuberoot grid = " << useQuad
		<< ", alphas = " << alphas.size() << std::endl;
	std::cout << "TOTAL production jobs = " << totalProd << " + " << checks.size() << " validation jobs" << std::endl;
	for (size_t i = 0; i < fns.size(); ++i) std::cout << "  jobs for f=" << fns[i].name << ": " << jobsPerF[i] << std::endl;
	mpreal::set_default_prec(savedPrec);
	struct Built { std::string name; std::vector<std::string> labels; std::function<std::vector<mpreal>(int)> compute; };
	auto build = [&](const Job& j)
	{
		Built b;
		const TFn& fn = fns[j.fi];
		const int form = fn.form;
		const WSpec& w = ws[j.wi];
		const bool isChk = j.ci >= 0;
		const Alpha* al = isChk ? &specs[j.ci].al : (j.ai >= 0 ? &alphas[j.ai] : nullptr);
		const ColSpec cs = colSpec(form, w.name, j.wc != 0);
		b.labels.push_back(form == GEN ? "L" : "logz");
		for (const auto& c : cs) b.labels.push_back((c.second ? "w*" : "") + c.first);
		b.name = std::string(isChk ? "[CHECK] " : "") + "f=" + fn.name + (al ? " alpha=" + al->name + " branch=" + fn.bl[j.br] : std::string(" (alpha absorbed)"))
			+ " w=" + w.name + (j.wc ? " [const]" : "") + " n=" + std::to_string(b.labels.size()) + (isChk ? " " + specs[j.ci].desc : std::string());
		const InvFn inv = fn.inv;
		const auto av = al ? al->val : std::function<mpreal()>();
		const auto wv = w.val;
		const int br = j.br;
		std::function<mpreal()> zfn = []() { return mpfr::zeta(mpreal(3)); };
		if (isChk)
		{
			const ChkSpec& sp = specs[j.ci];
			const auto a0 = sp.a0, x0 = sp.x0;
			const long bp = sp.bp, bq = sp.bq;
			zfn = [=]() { return inv(av())[br] / mpfr::pow(a0() * x0(), mpreal(bp) / bq * wv()); };
		}
		const std::map<std::string, std::function<mpreal()>>* basep = &base;
		b.compute = [=](int)
		{
			mpreal first = mpfr::log(zfn());
			if (form == GEN) first = mpfr::log(inv(av())[br]) - first;
			std::vector<mpreal> v = { first };
			const mpreal wvv = wv();
			for (const auto& c : cs)
			{
				const mpreal x = basep->at(c.first)();
				v.push_back(c.second ? wvv * x : x);
			}
			return v;
		};
		return b;
	};

	// ---- resume: skip production jobs that a previous output shows as finished with status "none" ----
	// Key = full job name (f, alpha, branch, w, [const], n) = text between "JOB i/N " and " | kind | normBound ...". The counter is ignored.
	// Lines with [CHECK] and lines with any status other than none are never skipped.
	if (const char* skipPath = std::getenv("ZETA3_TE_SKIP_DONE"))
	{
		std::map<std::string, int> nameCount;
		for (const auto& j : prod) ++nameCount[build(j).name];
		size_t dupNames = 0;
		for (const auto& kv : nameCount) if (kv.second > 1) ++dupNames;
		std::map<std::string, double> done; // name -> normBound
		std::set<std::string> ambiguous;
		std::ifstream in(skipPath);
		std::string line;
		size_t oldLines = 0;
		while (std::getline(in, line))
		{
			if (line.compare(0, 4, "JOB ") != 0) continue;
			const size_t sp = line.find(' ', 4);
			if (sp == std::string::npos) continue;
			const size_t sep = line.find(" | ", sp + 1);
			if (sep == std::string::npos) continue;
			const std::string name = line.substr(sp + 1, sep - sp - 1);
			if (name.compare(0, 8, "[CHECK] ") == 0) continue;
			const size_t sep2 = line.find(" | ", sep + 3);
			if (sep2 == std::string::npos) continue;
			const std::string kind = line.substr(sep + 3, sep2 - sep - 3);
			const size_t nb = line.find("normBound ", sep2);
			if (nb == std::string::npos) continue;
			++oldLines;
			if (kind != "none") { ambiguous.insert(name); continue; }
			if (done.count(name)) ambiguous.insert(name);
			done[name] = std::strtod(line.c_str() + nb + 10, nullptr);
		}
		// conservative: a name that also appears with a non-none status (or twice) is rerun
		size_t skipped = 0;
		double minSkippedNb = 1e300;
		std::vector<Job> rest;
		for (const auto& j : prod)
		{
			const std::string nm = build(j).name;
			const auto it = done.find(nm);
			if (it != done.end() && !ambiguous.count(nm) && nameCount[nm] == 1) { ++skipped; minSkippedNb = std::min(minSkippedNb, it->second); }
			else rest.push_back(j);
		}
		prod.swap(rest);
		std::cout << "RESUME from " << skipPath << ": old job lines parsed = " << oldLines << ", duplicate production names in this run = " << dupNames
			<< ", skipped = " << skipped << ", remaining = " << prod.size() << std::endl;
		std::cout << "RESUME skipped jobs: min normBound = " << minSkippedNb << " (combine with the min of this run's summary)" << std::endl;
	}
	if (countOnly) return;
	if (sampleJobs > 0 && (size_t)sampleJobs < prod.size())
	{
		std::vector<Job> sub;
		for (size_t i = 0; i < (size_t)sampleJobs; ++i) sub.push_back(prod[i * prod.size() / sampleJobs]);
		prod.swap(sub);
	}

	// ---- run: validation first, then production ----
	std::vector<Job*> all;
	for (auto& j : checks) all.push_back(&j);
	for (auto& j : prod) all.push_back(&j);
	std::mutex mtx;
	std::atomic<size_t> next(0);
	auto worker = [&]()
	{
		for (;;)
		{
			const size_t i = next++;
			if (i >= all.size()) return;
			Job& j = *all[i];
			const Built bt = build(j);
			j.name = bt.name;
			const auto t0 = std::chrono::steady_clock::now();
			std::ostringstream os;
			os.precision(30);
			j.oc = reportLinearRelationSearch(os, bt.name, bt.compute, bt.labels, digits);
			j.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
			bool l0 = true;
			if (j.oc.found)
			{
				l0 = j.oc.coeffs[0].find_first_not_of("-0") == std::string::npos;
				j.kind = !j.oc.verified ? "unverified" : (l0 ? "degenerate" : "real");
			}
			else j.kind = "none";
			std::string chkMsg;
			if (j.ci >= 0)
			{
				bool ok = j.oc.found && j.oc.verified && !l0 && j.oc.coeffs.size() == bt.labels.size();
				if (ok)
				{
					std::map<std::string, long> ex(specs[j.ci].exp.begin(), specs[j.ci].exp.end());
					const long s = std::atol(j.oc.coeffs[0].c_str()) > 0 ? 1 : -1;
					for (size_t k = 0; k < bt.labels.size(); ++k)
					{
						const auto it = ex.find(bt.labels[k]);
						if (std::atol(j.oc.coeffs[k].c_str()) != s * (it == ex.end() ? 0 : it->second)) ok = false;
					}
				}
				j.checkPass = ok;
				chkMsg = std::string("  check: ") + (ok ? "DETECTED, VERIFIED, coefficients match expected" : "FAILED (not detected or coefficients differ)");
			}
			std::lock_guard<std::mutex> lk(mtx);
			if (j.ci >= 0 || j.oc.found) std::cout << os.str();
			std::cout << "JOB " << (i + 1) << "/" << all.size() << " " << bt.name << " | " << j.kind << " | normBound " << j.oc.normBound
				<< " | " << j.secs << " s" << std::endl;
			if (!chkMsg.empty()) std::cout << chkMsg << std::endl;
		}
	};
	std::vector<std::thread> th;
	for (int t = 0; t < nThreads; ++t) th.emplace_back(worker);
	for (auto& t : th) t.join();

	// ---- summary ----
	std::cout << "\n===== SUMMARY (" << digits << " digits, variant " << variant << ") =====" << std::endl;
	int nReal = 0, nDeg = 0, nUnv = 0;
	double minAll = 1e300;
	std::vector<double> minF(fns.size(), 1e300), secF(fns.size(), 0);
	std::vector<int> cntF(fns.size(), 0), realF(fns.size(), 0), degF(fns.size(), 0);
	for (const auto& j : prod)
	{
		++cntF[j.fi]; secF[j.fi] += j.secs;
		if (j.kind == "real") { ++nReal; ++realF[j.fi]; std::cout << "REAL RELATION: " << j.name << std::endl; }
		else if (j.kind == "degenerate") { ++nDeg; ++degF[j.fi]; std::cout << "DEGENERATE (coefficient of first column is 0): " << j.name << std::endl; }
		else if (j.kind == "unverified") { ++nUnv; std::cout << "UNVERIFIED: " << j.name << std::endl; }
		else
		{
			const double nb = std::strtod(j.oc.normBound.c_str(), nullptr);
			minAll = std::min(minAll, nb); minF[j.fi] = std::min(minF[j.fi], nb);
		}
	}
	std::cout << "jobs run = " << prod.size() << ", real relations = " << nReal << ", degenerate = " << nDeg << ", unverified = " << nUnv << std::endl;
	std::cout << "min normBound (jobs without relation) overall = " << minAll << std::endl;
	for (size_t i = 0; i < fns.size(); ++i)
		std::cout << "  f=" << fns[i].name << ": jobs " << cntF[i] << ", real " << realF[i] << ", degenerate " << degF[i]
		<< ", min normBound " << minF[i] << ", mean time " << (cntF[i] ? secF[i] / cntF[i] : 0.0) << " s" << std::endl;
	std::cout << "--- validation ---" << std::endl;
	int chkFail = 0;
	for (const auto& j : checks)
	{
		if (!j.checkPass) ++chkFail;
		std::cout << (j.checkPass ? "PASS " : "FAIL ") << j.name << std::endl;
	}
	std::cout << "validation: " << (chkFail == 0 ? "all PASS" : "some FAIL") << std::endl;
	std::cout << "total elapsed " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count() << " s" << std::endl;
}

void test_pslq_zeta3_transcendental_exponent_absorbed(int digits, int maxPrime, int maxPrimeLogLog)
{
	zeta3TranscendentalExponentImpl(digits, 2, true, maxPrime, maxPrimeLogLog);
}

void test_pslq_zeta3_transcendental_exponent_grid(int digits, int maxPrime, int variant)
{
	zeta3TranscendentalExponentImpl(digits, variant, false, maxPrime, 3);
}
