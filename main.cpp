// HelloGMP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "test.h"

int main()
{
	//test_check_recurring_decimal_performance();

	//test_pslq_apery_over_pi_cubed();

	//test_pslq_apery_multiplicative();

	//test_pslq_apery_algebraic_exp_pi_grid();

	//test_pslq_zeta3_zeta5_ratio(3, 5);

	//test_pslq_zeta3_zeta5_ratio(2, 4);

	//test_pslq_zeta3_zeta5_ratio_parallel(100);

	//test_pslq_zeta3_zeta7_over_zeta5_squared(700, 50);

	//test_pslq_zeta3_zeta7_over_zeta5_squared_exp_grid(700, 12, 12, 20);

	//test_pslq_zeta3_zeta7_over_zeta5_squared_multiplicative(700, 23);

	//test_pslq_zeta3_zeta7_over_zeta5_squared_multiplicative(2000, 100);

	//test_pslq_zeta357_multiplicative(2000, 100);

	//test_pslq_zeta_multiplicative(2000, 100, { 3, 5, 7, 9, 11 }, true);

	//test_pslq_zeta_multiplicative(2000, 100, { 3, 5, 7, 9, 11 }, true, true);

	//test_pslq_zeta3_additive_elementary(3000, false);

	//test_pslq_zeta3_additive_elementary(3000, true);

	//test_pslq_zeta_weight_homogeneous(2000, 6, 24, 23);

	//test_pslq_zeta_central_binomial(3000, 3);

	//test_pslq_zeta_central_binomial(3000, 5);

	//test_pslq_zeta3_elementary_functions(1000, 12);

	//test_pslq_zeta3_elementary_functions_relations(2000);

	//test_pslq_zeta3_composite_functions(1000, 8, true, true);

	//test_pslq_zeta3_elementary_functions_highdeg(2000, 24, {0, -2, -3});

	//test_pslq_zeta3_composite_functions_highdeg(1500, 16, {0, -2, -3});

	//test_pslq_zeta3_multilevel_composite(1800, 12, {0, -1, -2, -3}, 3, 4);

	//test_pslq_zeta3_transcendental_exponent_absorbed(3000, 31, 13);

	//test_pslq_zeta3_transcendental_exponent_grid(800, 13);

	// test_pslq_zeta3_transcendental_exponent_polynomial(600, 12);
	test_pslq_zeta3_transcendental_exponent_polynomial_p37(400, 12);
	//test_pslq_zeta3_transcendental_exponent_absorbed_ext(3000, 31, 31);

	//test_mpreal();

	//test_multithreading_trigonometric_functions();
}
