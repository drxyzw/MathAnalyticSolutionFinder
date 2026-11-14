# MathAnalyticSolutionFinder

## 研究の目的
1. **第一目的：ζ(3) の解析的表式を求める。**
   例えば π³ × (何らかの代数的数) のような形。ただし、この例は望みが薄いことが判明している。
   解析的計算ではなく、高精度の小数表示を使って数値的に探索する。
2. **第二目的（ζ(3) 単独が難しい場合）：複数の奇数ゼータ値の間の関係式を探す。**
   例えば ζ(3)ζ(7) = (代数的数)·ζ(5)²。ただし、この例も望みが薄いことが判明している。

手法は PSLQ に限らない。

## エージェントの役割
- **メインセッション（リード, Opus/Fable）：** ユーザーと対話して、研究方針と実験を決める。詳細は下記の `doc/research/lead.md` を参照。
- **zeta-coder（Sonnet）：** 実装と短い試し実行を担当する。定義は `.claude/agents/zeta-coder.md`。

## ビルド
```
g++ -O2 -std=c++17 main.cpp test.cpp check_recurring_decimal.cpp pslq.cpp -o main -lmpfr -lgmp -pthread
```
- ビルド設定は `.vscode/tasks.json` と同じ（ただし tasks.json は `-g`）。
- 新しいソースファイルを追加したら、このコマンドと tasks.json の両方を更新すること。
- ライブラリは MPFR と GMP。mpreal ラッパーは `#include <mpreal.h>`、または gmp_util.h を使う。

## データ
- `data/zeta/zeta3.txt` 〜 `zeta99.txt`：奇数ゼータ値、10000桁。元データは `data/zeta/original/zeta3to99.txt`。
- `mpfr::zeta(mpreal(s))`、`mpfr::const_pi()`、`mpfr::const_euler()`、`mpfr::li2` などで、任意の精度の値を直接計算することもできる。

## 既存の道具（作り直さずに再利用する）
- [pslq.h](pslq.h)：
  - `pslq(x, maxNormBound, maxIter)`
  - `findIntegerPolynomial(x, degree)`：(1, x, …, x^d) に対して PSLQ を行う。
  - 戻り値の `status` は Found / PrecisionExhausted / NormBoundReached / MaxIterReached のいずれか。
  - `normBound` は「これより小さいノルムの整数関係は存在しない」ことを意味する。
  - 信頼する桁は精度の80%まで。偽の関係はガードで除外される。
- [test.cpp](test.cpp) のヘルパー関数：
  - `reportPolynomialSearch(name, x, maxDegree)`：次数 1..maxDegree を順に調べる。
  - `reportPolynomialSearchParallel(targets, digits, maxDegree)`：(ターゲット, 次数) の組を並列に調べる。見つかった関係は2倍の精度で検算する。
  - `reportAlgebraicExpPiSearch(name, xAtPrec, digits, degree, maxDen, maxAbs, piGrid)`：x = α·e^r·π^s の形をグリッドで探索する。piGrid=false にすると s=0 に固定する。
  - `reportMultiplicativeSearch(name, x, basis)`：(log x, log 基底…) に対して PSLQ を行う。
  - `reportLogRelationSearch(name, basis)`：∏ b_i^(c_i) = 1 の形の関係を探す。
  - `test_pslq_zeta_multiplicative(digits, maxPrime, zetaArgs, includePi, includeEuler)`：ζ値・π・γ・e・素数・単数の間の積の関係を探す。
- [check_recurring_decimal.h](check_recurring_decimal.h)：循環小数かどうかを判定する。
- [grid_search.cpp](grid_search.cpp)：グリッド探索（独立したプログラム）。
- [cf_search.cpp](cf_search.cpp)：多項式連分数（Ramanujan Machine 型）の探索。独立プログラム。`./cf_search <F1|F2|F3> <C> <M> <threads> [stride]`、`./cf_search F4 <C> <M> <threads> <Aalpha> <Ab> [stride]`（a_n=-(αn²+βn+γ)³）。ビルド：`g++ -O2 -std=c++17 cf_search.cpp -o cf_search -lmpfr -lgmp -pthread`

## 決まりごと
- 新しい実験は関数として追加する。
  - PSLQ 系の実験は `test_pslq_*` として [test.cpp](test.cpp) に書き、[test.h](test.h) で宣言する。
  - PSLQ 以外の手法は、必要に応じて新しいファイルにしてよい。
- [main.cpp](main.cpp) では、前回の呼び出しをコメントアウトして残し、今回の呼び出しを1行だけ有効にする。
- **必ず、正解が分かっている検証用のケースを同じ形で含める。** 正しく検出できることを確認してから本題を調べる。
- 結果は `results/<内容>_<digits>digits_<条件>.txt` に保存する。

## 解釈のルール
- PSLQ の判定力の目安：精度 D 桁、ベクトル長 n のとき、係数ノルムがおよそ 10^(0.8·D/n) までの関係を判定できる。
- 「見つからなかった」は証明ではない。ノルムの下限など「どこまで除外できたか」と一緒に報告する。
- 関係が見つかった場合は、必ずより高い精度で検算する。
- PSLQ は「ζ の係数が少なくとも1つ非ゼロ」という条件を課さない。ただし、ノルムの下限はすべての関係に対して成り立つ。したがって「関係なし」という結論はこの条件で弱まらない。

## 計算資源
- 10コア（`nproc`）。
- `std::thread::hardware_concurrency()` は16を返すため、スレッドの詰め込みすぎに注意する。
- 長い計算は、出力をファイルにリダイレクトしてバックグラウンドで実行する。

## リード用の指示
以下はメインセッション（リード）だけに向けた指示。サブエージェントは無視してよい。

@doc/research/lead.md
