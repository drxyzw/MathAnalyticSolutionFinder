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
	void reportAlgebraicExpPiSearch(const char* name, XAtPrec xAtPrec, int pslqDigits, int degree, int maxDen, int maxAbs)
	{
		std::cout << "=== " << name << " ===" << std::endl;
		mpreal::set_default_prec(mpfr::digits2bits(pslqDigits));
		const mpreal x = xAtPrec(pslqDigits);
		const mpreal pi = mpfr::const_pi();
		const std::vector<Rational> rats = enumerateRationals(maxDen, maxAbs);
		const size_t nr = rats.size();
		std::vector<mpreal> expMinus(nr), piMinus(nr);
		for (size_t i = 0; i < nr; ++i)
		{
			const mpreal q = mpreal(rats[i].num) / rats[i].den;
			expMinus[i] = mpfr::exp(-q);
			piMinus[i] = mpfr::pow(pi, -q);
		}
		const size_t total = nr * nr;
		std::cout << "degree <= " << degree << ", " << pslqDigits << " digits, denominators <= " << maxDen
			<< ", |r|,|s| <= " << maxAbs << ": " << total << " (r, s) pairs" << std::endl;

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
						const mpreal y = x * expMinus[k / nr] * piMinus[k % nr];
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
			const Rational& er = rats[k / nr];
			const Rational& ps = rats[k % nr];
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
				const mpreal y = x * expMinus[k / nr] * piMinus[k % nr];
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
				<< " with |coeffs|_2 < " << minBound.toString("%.3Re") << " (weakest at r = " << rats[minBoundIdx / nr].str()
				<< ", s = " << rats[minBoundIdx % nr].str() << ")" << std::endl;
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
