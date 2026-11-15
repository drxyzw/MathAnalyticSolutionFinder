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
void test_pslq_zeta_weight_homogeneous(int digits, int minWeight, int maxWeight, int maxZetaArg); // homogeneous weight-w monomials in pi^2, zeta(odd)
void test_pslq_zeta3_additive_elementary(int digits, bool includeLog5); // zeta(3) vs weight-3 monomials in pi, log2, log3, (log5), Catalan G
void test_pslq_zeta_central_binomial(int digits, int targetZeta); // zeta(3) or zeta(5) vs sum s^n/(n^w C(2n,n)^k) and weight-w constants
void test_pslq_zeta3_elementary_functions_relations(int digits); // multiplicative/additive relations of f(zeta3*pi^k)
void test_pslq_zeta3_elementary_functions(int digits, int maxDegree); // f(zeta(3))*pi^k algebraic? f in exp, log, trig, hyperbolic, ...
void test_pslq_zeta3_elementary_functions_highdeg(int digits, int degree, std::vector<int> ks); // f(zeta3*pi^k), only degree==`degree` is run
void test_pslq_zeta3_composite_functions_highdeg(int digits, int degree, std::vector<int> ks); // f(g(zeta3*pi^k)), algebraicity only, degree==`degree` only
void test_pslq_zeta3_multilevel_composite(int digits, int degree, std::vector<int> ks, int minLevel, int maxLevel); // y = f1(f2(..fL(zeta3*pi^k))), L=minLevel..maxLevel, algebraicity, degree==`degree` only
void test_pslq_zeta3_composite_functions(int digits, int maxDegree, bool doAlgebraic, bool doRelations); // y = f(g(zeta3*pi^k)): algebraic? and multiplicative/additive relations
