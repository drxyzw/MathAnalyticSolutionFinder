# 実験ログ

目的：1 = ζ(3) の表式、2 = 奇数ゼータ値の間の関係。

| 日付 | 目的 | 調べた形 | 手法 | パラメータ | 結果 | どこまで除外できたか | 結果ファイル |
|---|---|---|---|---|---|---|---|
| 〜2026-11 | 1 | ζ(3)/π³ が代数的数 | 多項式 PSLQ | 1000桁、次数≤20 | なし（前セッション、結果ファイルなし） | – | test.cpp: test_pslq_apery_over_pi_cubed |
| 〜2026-11 | 1 | ζ(3)/π³ = α·e^r·π^s | グリッド + 多項式 PSLQ | 500桁、次数≤8、分母≤6、\|r\|,\|s\|≤6 | なし（前セッション、結果ファイルなし） | – | test.cpp: test_pslq_apery_algebraic_exp_pi_grid |
| 〜2026-11 | 1 | ζ(3)/π³ = α·e^r·π^s（大規模） | grid_search.cpp | 100桁、次数≤10、分母≤10、\|r\|,\|s\|≤100（\|log10 α\|≤7.3 の 5,979,257 組） | ヒット0件 | 約10⁷ | results/log.txt, results/results.tsv, results/hits.txt |
| 〜2026-11 | 1 | ζ(3)/π³ = e^r·π^s·∏p^a·単数 | log PSLQ | 1000桁、素数≤23 | なし（前セッション、結果ファイルなし） | – | test.cpp: test_pslq_apery_multiplicative |
| 2026-11-14 | 2 | ζ(3)^(1/3)/ζ(5)^(1/5) が代数的数 | 多項式 PSLQ（並列） | 500桁、次数≤100 | なし | 次数4で10⁸⁰、次数24で4.9×10¹⁵、次数100で7 | results/zeta3_zeta5_ratio_deg100.txt |
| 2026-11-14 | 2 | ζ(3)^(1/2)/ζ(5)^(1/4) が代数的数 | 多項式 PSLQ（並列） | 500桁、次数≤100 | なし | 次数4で10⁸⁰、次数24で6.3×10¹⁵、次数100で6.6 | results/zeta3_zeta5_ratio_deg100.txt |
| 2026-11-14 | 2 | ζ(3)ζ(7)/ζ(5)² が代数的数 | 多項式 PSLQ（並列） | 700桁、次数≤50 | なし | 次数1で2×10²⁸⁰、次数50で2.9×10⁹ | results/zeta3_zeta7_over_zeta5sq_700digits_deg50.txt |
| 2026-11-14 | 2 | ζ(3)ζ(7)/ζ(5)² = α·e^r | グリッド + 多項式 PSLQ | 700桁、次数≤12、分母≤12、\|r\|≤6 | なし | 1.2×10⁴³ | results/zeta3_zeta7_over_zeta5sq_exp_grid_700digits.txt |
| 2026-11-14 | 2 | 同上 | 同上 | \|r\|≤20 | なし | 1.2×10⁴³ | results/zeta3_zeta7_over_zeta5sq_exp_grid_700digits_r20.txt |
| 2026-11-14 | 2 | ζ(3)ζ(7)/ζ(5)² = e^r·∏p^a·単数 | log PSLQ | 700桁、素数≤23 | なし | 1.0×10⁴⁰ | results/zeta3_zeta7_over_zeta5sq_multiplicative_700digits.txt |
| 2026-11-14 | 2 | 同上 | log PSLQ | 2000桁、素数≤100 | なし | 2.4×10⁵³ | results/zeta3_zeta7_over_zeta5sq_multiplicative_2000digits_p100.txt |
| 2026-11-14 | 2 | ζ(3)^a ζ(5)^b ζ(7)^c = e^r·∏p^a·単数 | log PSLQ | 2000桁、素数≤100 | なし | 1.1×10⁵⁰ | results/zeta357_multiplicative_2000digits_p100.txt |
| 2026-11-14 | 2 | ζ(3,5,7,11)、π の積の関係 | log PSLQ | 2000桁、素数≤100 | なし | 1.2×10⁴⁷ | results/zeta3_5_7_11_pi_multiplicative_2000digits_p100.txt |
| 2026-11-14 | 2 | ζ(3,5,7,9,11)、π の積の関係 | log PSLQ | 2000桁、素数≤100 | なし | 5.2×10⁴⁵ | results/zeta3_5_7_9_11_pi_multiplicative_2000digits_p100.txt |
| 2026-11-14 | 2 | ζ(3,5,7,9,11)、π、γ の積の関係 | log PSLQ | 2000桁、素数≤100 | なし | 3.0×10⁴⁴ | results/zeta3_5_7_9_11_pi_gamma_multiplicative_2000digits_p100.txt |
| 2026-11-14 | 1 | ζ(3) = π, log2, log3, G の重さ3の積（13個）の有理1次結合 | 線形 PSLQ | 3000桁、n=14 | なし（検証：Li₃(1/2)、Li₃(±1/3) の既知恒等式を再発見・6000桁で検算済み） | 2.8×10¹⁷¹ | results/zeta3_additive_elementary_3000digits_log23.txt |
| 2026-11-14 | 1 | 同上 + log5（24個） | 線形 PSLQ | 3000桁、n=25 | なし（検証同上） | 1.0×10⁹⁶ | results/zeta3_additive_elementary_3000digits_log235.txt |
| 2026-11-14 | 2 | π², ζ(3), …, ζ(23) の重さ w の同次多項式関係（w=6..24） | 線形 PSLQ（並列） | 2000桁、n=2〜66 | なし（全 w で検証ケース PASS：偶数 w は ζ(w)∝π^w、奇数 w は構成値） | w=10で1.2×10³²⁰、w=16で1.4×10⁹⁴、w=20で5.4×10⁴⁵、w=24で3.3×10²² | results/zeta_weight_homogeneous_2000digits_w6to24.txt |
| 2026-11-14 | 1 | ζ(3) の多項式連分数 F1：a_n=−n⁶、b_n=3次（\|c\|≤100）、値が (A+Bζ(3))/(C+Dζ(3))（\|A..D\|≤20） | cf_search（Ramanujan Machine 型） | 1.6×10⁹ 通り、210桁で検算 | 既知の Apéry 6/ζ(3) と Ramanujan Machine の 8/(7ζ(3)) のみ（自己テスト PASS） | 幾何収束（約0.07桁/項以上）で照合範囲内のもの | results/cf_search_zeta3_F1_C100_M20.txt |
| 2026-11-14 | 2 | ζ(5) の多項式連分数 F3：a_n=−n¹⁰、b_n=5次（\|c\|≤20） | cf_search | 4.75×10⁹ 通り | なし（このファミリーに既知例なし、照合機構は F1/F2 で検証済み） | 同上 | results/cf_search_zeta5_F3_C20_M20.txt |
| 2026-11-14 | 1 | ζ(3) の多項式連分数 F4：a_n=−(αn²+βn+γ)³（α≤3, \|β\|,\|γ\|≤3）、b_n=3次（\|c\|≤40） | cf_search | 5.5×10⁹ 通り | 既知の 8/(7ζ(3)) のみ（自己テスト：Apéry とそのシフト版 PASS） | 幾何収束で照合範囲内のもの | results/cf_search_zeta3_F4_C40_A3_3_M20.txt |
| 2026-11-14 | 1 | ζ(3) = Σ s^n/(n³C(2n,n)^k)（s∈{1,±2,±3}, k=1; s∈{±1,±4}, k=2）と π³, π²log2, log³2 の有理1次結合（Koecher 級数 S(−1,3,1) は除外） | 線形 PSLQ | 3000桁、n=13 | なし（検証：Koecher 2ζ(3)+5S(−1,3,1)=0、Σ1/(n²C)=π²/18 を再発見） | 4.9×10¹⁸⁴ | results/zeta3_central_binomial_3000digits.txt |
| 2026-11-14 | 2 | ζ(5) = 同型の重さ5の級数 + π²ζ(3), ζ(3)log²2, π⁴log2, π²log³2, log⁵2（Koecher 型を除外） | 線形 PSLQ | 3000桁、n=15 | なし（検証：Koecher の ζ(5) 式を再発見） | 1.0×10¹⁶⁰ | results/zeta5_central_binomial_3000digits.txt |
| 2026-11-14 | 1 | f(ζ(3)/π^k) が代数的数（f = exp, log, sin, cos, tan, sinh, cosh, tanh, arctan, arcsinh, arccosh, arcsin(1/t), arccos(1/t), log(t−1), t^t, Γ, exp(1/t)、k=−2..3、90ターゲット） | 多項式 PSLQ（並列） | 1000桁、次数≤12 | なし（検証17件すべて正しい次数で検出・検算） | 次数12で約3.5×10⁶¹、次数8で約8×10⁸⁸、次数4で約10¹⁶⁰ | results/zeta3_elementary_functions_1000digits_deg12.txt |
