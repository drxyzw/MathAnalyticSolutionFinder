#include <filesystem>
#include <string>
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
