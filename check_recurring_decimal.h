#pragma once

#if __has_include(<gmp_util.h>)
#include <gmp_util.h>
#else
#include <mpreal.h>
#endif

using mpfr::mpreal;

inline const int digits = 50;

bool checkRecurringDecimal(const mpreal& x, int& period, int& n);
bool checkRecurringDecimalFast(const mpreal& x, int& period, int& n);
