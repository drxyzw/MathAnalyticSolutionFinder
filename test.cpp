#include <filesystem>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <chrono>
#include "check_recurring_decimal.h"
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
