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

		PolySearchSummary reportPolynomialSearchParallel(const std::vector<PolySearchTarget>& targets, int pslqDigits, int maxDegree)
	{
		mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
		const size_t nt = targets.size();
		std::vector<mpreal> xs(nt);
		for (size_t t = 0; t < nt; ++t)
			xs[t] = targets[t].xAtPrec(pslqDigits);

		// jobs: (target, degree), biggest degree first for load balancing
		std::vector<std::pair<size_t, int>> jobs;
		for (int d = maxDegree; d >= 1; --d)
			for (size_t t = 0; t < nt; ++t)
				jobs.push_back({ t, d });
		std::vector<std::vector<PslqResult>> results(nt, std::vector<PslqResult>(maxDegree + 1));
		std::vector<std::vector<double>> seconds(nt, std::vector<double>(maxDegree + 1));

		std::atomic<size_t> next{ 0 }, done{ 0 };
		const unsigned nThreads = std::max(1u, std::thread::hardware_concurrency());
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
							<< toString(results[t][d].status) << " (" << static_cast<long>(seconds[t][d] * 1000) << " ms)" << std::endl;
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
			for (int d = 1; d <= maxDegree; ++d)
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
				std::cout << (genuine ? "  VERIFIED" : "  SPURIOUS") << "\n    " << polynomialToString(r.relation)
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
struct LinearSearchOutcome { bool found = false; bool verified = false; };
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
