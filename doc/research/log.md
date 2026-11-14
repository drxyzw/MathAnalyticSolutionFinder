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
