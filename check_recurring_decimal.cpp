#include <string>
#include <algorithm>
#include <gmpxx.h>
#include "check_recurring_decimal.h"

bool checkRecurringDecimal(const mpreal& x, int& period, int& n) {
	// suppose x = M + 0.a(1)a(2)a(3).....a(n-1) a(n) a(n+1)...a(n+p-1) a(n+p)....a(N)
	// M: integer
	// n-1: non-recurring digits
	// repeated part: a(n) ... a(n+p-1)
	// p: period
	// N: digits of fraction part: precision digits minus # of digits of M minus ERROR_BUFFER
	//
	// Transformation:
	// x --> x / 10^[# of digits of M]
	// N --> [precision digits]
	// n-1 --> non-recurring digits of updated x
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
	// where {X} = X - [X]* and [X]* is the closest integer
	// This is a condition for a recurring decimal
	//
	// At least, we require 2 periods to be detected as recurring decimal:
	// (non-recurring digits) + (2 periods)
	// n-1 + 2p <= N ie n <= N-2p+1 (3)
	// When (1) and (3) are compared, because p >= 1, the condition (2) is strictor: n_max = N-2p
	// Hence n_max = N-2p+1
	// Another condition is MIN_RECURRING_DIGITS
	// n + MIN_RECURRING_DIGITS <= N
	// n_max = N - MIN_RECURRING_DIGITS
	// n_max = min(N-2p+1, N - MIN_RECURRING_DIGITS)
	//
	// p <= MAX_PERIOD
	//
	// Because the final digit has truncation error and + multiplication / division
	// --> erro second last digit --> keep 3 digits for a buffer
	int ERROR_BUFFER = 3;
	const mpreal zero = (mpreal)(0.);
	const mpreal one = (mpreal)(1.);
	const mpreal ten = (mpreal)(10.);
	const mpreal zero_one = (mpreal)(0.1);
	const mpreal digits_x_int = mpfr::max(zero, one + mpfr::floor(mpfr::log10(mpfr::abs(x))));
	const mpreal x_tf = x * mpfr::exp10(-digits_x_int);
	const int N = digits - ERROR_BUFFER;
	mpreal ten_p = one; // initial value
	mpreal ten_nm1;
	mpreal ten_nmNpm1;
	const int MIN_RECURRING_DIGITS = (int)(N/2);
	const int MAX_PERIOD = (int)(N/2);
	for (int p_i = 1; p_i <= MAX_PERIOD; ++p_i) {
		ten_p *= ten;
		const int n_max = std::min(N - MIN_RECURRING_DIGITS, N - 2*p_i + 1);
		ten_nm1 = zero_one;
		ten_nmNpm1 = mpfr::exp10(-N+p_i-1);
		for(int n_j = 1; n_j <= n_max; ++n_j) {
			ten_nm1 *= ten;
			ten_nmNpm1 *= ten;
			const mpreal z = (ten_p - one) * x_tf * ten_nm1;
			//std::cout << z << std::endl;
			const mpreal z_fraction = z - mpfr::round(z); // remove integer part
			if (-ten_nmNpm1 < z_fraction && z_fraction < ten_nmNpm1) {
				period = p_i;
				n = n_j;
				return true;
			}
		}
	}
	period = 0;
	n = 0;
	return false;
}

// Same test as checkRecurringDecimal, done in exact integer arithmetic.
// In checkRecurringDecimal, x_tf = x / 10^D can be negative and z = (10^p-1) * x_tf * 10^(n-1).
// Here y = |x| / 10^D >= 0, so z' = (10^p-1) * y * 10^(n-1) = |z|.
// |{z}| = |{z'}| because the distance to the nearest integer is the same for z and -z,
// so the test can use z' >= 0.
// A has target digits as integer
// A = round(y * 10^digits) >= 0, B = (10^p-1) * A >= 0.
// z' * 10^(digits-n+1) = (10^p-1) * y * 10^digits = B (up to the rounding of A), so
// |{z'}| < 10^(n-N+p-1)  <=>  B is within 10^j of a multiple of 10^k
// where j = p + ERROR_BUFFER and k = digits - n + 1,
// i.e. the decimal digits of B at positions j..k-1 (from the right) are all 0 or all 9.
// The smallest n comes from the longest such run, so no loop over n is needed.
bool checkRecurringDecimalFast(const mpreal& x, int& period, int& n) {
	const int ERROR_BUFFER = 3;
	const mpreal zero = (mpreal)(0.);
	const mpreal one = (mpreal)(1.);
	const mpreal digits_x_int = mpfr::max(zero, one + mpfr::floor(mpfr::log10(mpfr::abs(x))));
	const mpreal y = mpfr::abs(x) * mpfr::exp10(-digits_x_int);
	const mpreal scaled = y * mpfr::exp10(mpreal(digits));
	mpz_class A;
	mpfr_get_z(A.get_mpz_t(), scaled.mpfr_srcptr(), MPFR_RNDN);
	const int N = digits - ERROR_BUFFER;
	const int MIN_RECURRING_DIGITS = (int)(N/2);
	const int MAX_PERIOD = (int)(N/2);
	mpz_class ten_p = 1;
	for (int p_i = 1; p_i <= MAX_PERIOD; ++p_i) {
		ten_p *= 10;
		const int n_max = std::min(N - MIN_RECURRING_DIGITS, N - 2*p_i + 1);
		const mpz_class B = A * ten_p - A;
		const std::string s = B.get_str();
		// j-th digit from the right; positions above the top digit are 0
		auto digit_at = [&s](int j) { return j < (int)s.size() ? s[s.size() - 1 - j] : '0'; };
		const int j = p_i + ERROR_BUFFER;
		const char c = digit_at(j);
		if (c != '0' && c != '9')
			continue;
		int k = j;
		while (k < digits && digit_at(k) == c)
			++k;
		const int n_j = digits + 1 - k;
		if (n_j <= n_max) {
			period = p_i;
			n = n_j;
			return true;
		}
	}
	period = 0;
	n = 0;
	return false;
}
