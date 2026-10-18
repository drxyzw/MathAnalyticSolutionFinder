// HelloGMP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <filesystem>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#if __has_include(<gmp_util.h>)
#include <gmp_util.h>
#else
#include <mpreal.h>
#endif

using mutex_t = std::mutex;
mutex_t mutex;
const int digits = 500;
using mpfr::mpreal;

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

bool checkRecurringDecimal(const mpreal& x, unsigned& period) {
	// suppose x = M + 0.a(1)a(2)a(3).....a(n-1) a(n) a(n+1)...a(n+p-1) a(n+p)....a(N)
	// M: integer
	// n-1: non-recurring digits
	// repeated part: a(n) ... a(n+p-1)
	// p: period
	// N: digits of fraction part: precision digits minus # of digits of M
	//
	// 10^p x = M1 + 0.a(p+1) a(p+2) a(p+3)..a(n-1+p) a(n+p) a(n+p+1)..a(n+2p+1)a(n+2p)...aN     0    ....   0
	//      x = M  + 0.a1     a2     a3..... a(n-1)   an     a(n+1)... a(n+p-1) a(n+p)....a(N-p) a(N-p+1)... aN
	//                                                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ cancel each other
	// where M1 is an integer
	// 
	// (10^p-1)x = M2+ M3 * 10^(-n+1) - a(N-p+1)*10^(-N+p-1)  - a(N-p+2)*10^(-N+p-2)  ... - aN*10^(-N)
	// (10^p-1)x*10^(n-1) = M4        - a(N-p+1)*10^(n-N+p-2) - a(N-p+2)*10^(n-N+p-3) ... - aN*10^(n-N-1)
	// where M2, M3, and M4 are integers
	// {(10^p-1)x*10(n-1)} = - a(N-p+1)*10^(n-N+p-2) - a(N-p+2)*10^(n-N+p-3) ... - aN*10^(n-N-1)
	// (as long as n-N+p-2 <= -1 ie n <= N-p+1 (1))
	// |{(10^p-1)x*10^(n-1)}| < 10^(n-N+p-1) (2)
	// This is a condition for a recurring decimal
	// at least, we require 2 periods to be detected as recurring decimal:
	// (non-recurring digits) + (2 periods)
	// n-1 + 2p <= N ie n <= N-2p+1 (3)
	// When (1) and (3) are compared, because p >= 1, the condition (2) is strictor: n_max = N-2p
	// When n = n_max = N-2p+1,
	mpreal zero = (mpreal)(0.);
	const auto eff_digits = x >= 1. ? (digits - 1 - (int)(mpfr::max(zero, mpfr::floor(mpfr::log10(x))))) : digits;
	const auto N = (unsigned)(eff_digits);
	mpreal one = (mpreal)(1.);
	mpreal ten_N = mpfr::exp10(mpreal(N));
	mpreal ten_p;
	mpreal ten_nm1;
	mpreal ten_nmNpm1;
	for (unsigned p_i = 1; p_i <= digits/2; ++p_i) {
		ten_p = mpfr::exp10(mpreal(p_i));
		unsigned n_max = N - 2*p_i + 1;
		for(unsigned n = 0; n <= n_max; ++n) {
			ten_nm1 = mpfr::exp10(n-1);
			ten_nmNpm1 = mpfr::exp10(n-N+p_i-1);
			const auto z = (ten_p - one) * x * ten_nm1;
			//std::cout << z << std::endl;
			const auto z_fraction = z - mpfr::trunc(z); // remove integer part
			if (-ten_nmNpm1 < z_fraction && z_fraction < ten_nmNpm1) {
				period = p_i;
				return true;
			}
		}
	}
	period = 0;
	return false;
}
int main()
{
	mpreal::set_default_prec(mpfr::digits2bits(digits));
	std::cout.precision(digits);

	//const mpreal x = (mpreal)(7777768117.11) + (mpreal)(2.) / (mpreal)(71.) / (mpreal)(100.);
	const mpreal x = mpfr::const_pi();
	unsigned p;
	bool r = checkRecurringDecimal(x, p);
	std::cout << "x: " << x << std::endl;
	std::cout << "recurring decimal?: " << (r ? "true" : "false") << ", period: " << p << std::endl;
	//test_mpreal();

	//test_multithreading_trigonometric_functions();
}

