// HelloGMP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <filesystem>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <gmp_util.h>

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

bool checkRecurringDecimal(const mpreal x, unsigned &period) {
	// suppose x = M + 0.a1a2a3.....a(n-1) an a(n+1)...a(n+p-1) a(n+p)....aN
	// M: integer
	// repeated part: an ... a(n+p-1)
	// p: period
	// N: precision digits minus # of digits of M
	//
	// 10^p x = M1 + ap.a(p+1) a(p+2) a(p+3)..a(n-1+p) a(n+p) a(n+p+1)..a(n+2p+1)a(n+2p)...aN     0    ....   0
	//      x = M  +  0.a1     a2     a3..... a(n-1)   an     a(n+1)... a(n+p-1) a(n+p)....a(N-p) a(N-p+1)... aN
	//                                                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ cancel each other
	// where M1 is an integer
	// 
	// (10^p-1)x=M2+ M3 * 10^(-n+1) - a(N-p+1)*10^(-N+p-1)  - a(N-p+2)*10^(-N+p-2)  ... - aN*10^(-N)
	// (10^p-1)x*10(n-1) = M4       - a(N-p+1)*10^(n-N+p-2) - a(N-p+2)*10^(n-N+p-3) ... - aN*10^(n-N-1)
	// where M2, M3, and M4 are integers
	// {(10^p-1)x*10(n-1)} = - a(N-p+1)*10^(n-N+p-2) - a(N-p+2)*10^(n-N+p-3) ... - aN*10^(n-N-1)
	// (as long as n-N+p-2 <= -1 ie n <= N-p+1 (1))
	// |{(10^p-1)x*10^(n-1)}| < 10^(n-N+p-1)
	// This is a condition for a recurring decimal
	// at least, we require 2 periods: n + 2p <= N ie n <= N-2p (2)
	// Because p >= 1, the condition (2) is strictor: n_max = N-2p
	// When n = n_max,
	// |{(10^p-1)x*10^(N-2p-1)}| < 10^(-p-1)
	const auto eff_digits = digits - (int)(mpfr::floor(mpfr::log10(x)));
	const auto N = (mpreal)(eff_digits);
	mpreal one = (mpreal)(1.);
	mpreal ten = (mpreal)(10.);
	mpreal ten_p = (mpreal)(1.);
	mpreal ten_Nm1 = mpfr::pow((mpreal)(10.), (mpreal)(eff_digits - 1));
	//mpreal ten_Nm2pm1 = ((mpreal)(10.))^((mpreal)(digits - 1));
	mpreal zero_one = (mpreal)(0.1);
	mpreal one_hundred = (mpreal)(100.);
	for (unsigned p_i = 1; p_i <= digits/2; ++p_i) {
		ten_p *= ten;
		mpreal ten_Nm2pm1 = ten_Nm1 / (ten_p * ten_p);
		mpreal threshold = zero_one / ten_p;
		//const auto ten_p_m_one = (ten_p - one);
		//std::cout << ten_p_m_one << std::endl;
		//const auto x_ten_p_m_one = (ten_p - one) * x;
		//std::cout << x_ten_p_m_one << std::endl;
		const auto z = (ten_p - one) * x * ten_Nm2pm1;
		//std::cout << z << std::endl;
		const auto z_fraction = z - mpfr::trunc(z); // remove integer part
		if (-threshold < z_fraction && z_fraction < threshold) {
			period = p_i;
			return true;
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

