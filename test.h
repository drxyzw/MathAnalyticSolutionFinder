#pragma once

#include <vector>

void test_mpreal();
void test_multithreading_trigonometric_functions();
void test_check_recurring_decimal_performance();
void test_pslq_apery_over_pi_cubed();
void test_pslq_apery_multiplicative();
void test_pslq_apery_algebraic_exp_pi_grid();
void test_pslq_zeta3_zeta5_ratio(int p, int q); // x = zeta(3)^(1/p) / zeta(5)^(1/q)
void test_pslq_zeta3_zeta5_ratio_parallel(int maxDegree); // (p, q) = (3, 5) and (2, 4), degrees 1..maxDegree in parallel
void test_pslq_zeta3_zeta7_over_zeta5_squared(int pslqDigits, int maxDegree);
void test_pslq_zeta3_zeta7_over_zeta5_squared_exp_grid(int pslqDigits, int degree, int maxDen, int maxAbs); // x = alpha * e^r
void test_pslq_zeta3_zeta7_over_zeta5_squared_multiplicative(int pslqDigits, int maxPrime); // log x = r + sum c_i log(base_i)
void test_pslq_zeta357_multiplicative(int pslqDigits, int maxPrime); // zeta(3)^a zeta(5)^b zeta(7)^c = e^r * prod p^(a_p) * units
void test_pslq_zeta_multiplicative(int pslqDigits, int maxPrime, const std::vector<int>& zetaArgs, bool includePi, bool includeEuler = false);
