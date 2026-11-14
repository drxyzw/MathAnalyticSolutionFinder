// Ramanujan-Machine style search for polynomial continued fractions
//   x = b_0 + a_1/(b_1 + a_2/(b_2 + ...)),   a_n = sa * n^ea,  b_n = sum_k c_k n^k,  b_0 = c_0
// matched against Mobius transforms (A + B t)/(C + D t) of a target constant t.
//
// Usage: cf_search <F1|F2|F3> <C> <M> <threads> [stride]
//   F1: t = zeta(3), a_n = -n^6,        b cubic,   |c_i| <= C
//   F2: t = zeta(2), a_n = +-n^4,       b quadratic
//   F3: t = zeta(5), a_n = -n^10,       b quintic
//   stride (optional, default 1): run only every stride-th task (for timing estimates; not a full search)
// Build: g++ -O2 -std=c++17 cf_search.cpp -o cf_search -lmpfr -lgmp -pthread

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <tuple>
#include <vector>
#include <mpreal.h>

using mpfr::mpreal;
using ld = long double;

struct Family {
    std::string name;
    int K;          // degree of b
    int ea;         // exponent in a_n
    int sa;         // sign of a_n
    int tz;         // target zeta(tz)
};

struct Entry {
    double val;
    int8_t A, B, C, D;
};

static const int N1 = 200;   // terms, stage 1
static const int N2 = 600;   // terms, stage 2

static std::vector<ld> g_pw[8];   // g_pw[k][n] = n^k
static std::vector<ld> g_a;       // a_n
static std::vector<Entry> g_table;
static ld g_t;

static void init_powers(const Family& f) {
    for (int k = 0; k <= f.K; k++) {
        g_pw[k].assign(N2 + 1, 0);
        for (int n = 0; n <= N2; n++) g_pw[k][n] = powl((ld)n, (ld)k);
    }
    g_a.assign(N2 + 1, 0);
    for (int n = 0; n <= N2; n++) g_a[n] = f.sa * powl((ld)n, (ld)f.ea);
}

// Evaluate the CF with N terms. Returns false if rejected. Early-prunes slow convergence when prune=true.
static bool eval_cf(const ld* c, int K, int N, bool prune, ld& out, ld tol) {
    ld pm1 = 1, qm1 = 0;             // p_{-1}, q_{-1}
    ld p0 = c[0], q0 = 1;            // p_0, q_0
    ld prevRatio = 0;
    for (int n = 1; n <= N; n++) {
        ld b = 0;
        for (int k = 0; k <= K; k++) b += c[k] * g_pw[k][n];
        ld a = g_a[n];
        ld p = b * p0 + a * pm1;
        ld q = b * q0 + a * qm1;
        pm1 = p0; qm1 = q0; p0 = p; q0 = q;
        if ((n & 7) == 0) {
            ld m = fabsl(q0) > fabsl(p0) ? fabsl(q0) : fabsl(p0);
            if (m > 1e1000L) {
                ld s = 1 / m;
                p0 *= s; q0 *= s; pm1 *= s; qm1 *= s;
            }
        }
        if (prune && (n == 50 || n == 100)) {
            if (q0 == 0 || qm1 == 0) return false;
            ld r = p0 / q0, rp = pm1 / qm1;
            ld lim = (n == 50) ? 1e-2L : 1e-6L;
            if (!(fabsl(r - rp) <= lim * (1 + fabsl(r)))) return false;
        }
    }
    if (!(fabsl(q0) > 1e-300L) || !(fabsl(qm1) > 1e-300L)) return false;
    ld r = p0 / q0, rp = pm1 / qm1;
    if (!std::isfinite((double)r) || !std::isfinite((double)rp)) return false;
    if (!(fabsl(r - rp) <= tol * (fabsl(r) > 1 ? fabsl(r) : 1))) return false;
    out = r;
    (void)prevRatio;
    return true;
}

static void build_table(int M, ld t) {
    g_table.clear();
    for (int A = -M; A <= M; A++)
        for (int B = -M; B <= M; B++)
            for (int C = -M; C <= M; C++)
                for (int D = -M; D <= M; D++) {
                    if (C == 0 && D == 0) continue;
                    if ((long)A * D - (long)B * C == 0) continue;
                    int first = A ? A : (B ? B : (C ? C : D));
                    if (first < 0) continue;
                    int g = std::gcd(std::gcd(abs(A), abs(B)), std::gcd(abs(C), abs(D)));
                    if (g != 1) continue;
                    ld v = ((ld)A + B * t) / ((ld)C + D * t);
                    g_table.push_back({(double)v, (int8_t)A, (int8_t)B, (int8_t)C, (int8_t)D});
                }
    std::sort(g_table.begin(), g_table.end(), [](const Entry& x, const Entry& y) { return x.val < y.val; });
}

// Look up x in table with relative tolerance; returns matching entries.
static void lookup(ld x, double reltol, std::vector<Entry>& res) {
    double xd = (double)x;
    double tol = reltol * std::max(1.0, std::fabs(xd));
    auto it = std::lower_bound(g_table.begin(), g_table.end(), xd - tol,
                               [](const Entry& e, double v) { return e.val < v; });
    for (; it != g_table.end() && it->val <= xd + tol; ++it) res.push_back(*it);
}

struct Hit {
    std::vector<int> c;
    Entry e;
};

// High-precision verification.
struct VerifyResult {
    mpreal diff;
    double digitsPerTerm;
    bool ok;
};

static VerifyResult verify(const Family& f, const std::vector<int>& c, const Entry& e, const mpreal& t) {
    const int NN = 2000;
    mpreal pm1 = 1, qm1 = 0, p0 = c[0], q0 = 1;
    std::vector<std::pair<int, mpreal>> chk;
    for (int n = 1; n <= NN; n++) {
        long long bn = 0;
        long long pw = 1;
        for (int k = 0; k <= f.K; k++) { bn += (long long)c[k] * pw; pw *= n; }
        mpreal a = mpreal(f.sa) * pow(mpreal(n), f.ea);
        mpreal p = bn * p0 + a * pm1;
        mpreal q = bn * q0 + a * qm1;
        pm1 = p0; qm1 = q0; p0 = p; q0 = q;
        if (n % 10 == 0 && n <= 200 && q0 != 0) chk.push_back({n, p0 / q0});
    }
    VerifyResult r;
    r.digitsPerTerm = NAN;
    if (q0 == 0) { r.diff = 1; r.ok = false; return r; }
    mpreal X = p0 / q0;
    mpreal target = (mpreal(e.A) + e.B * t) / (mpreal(e.C) + e.D * t);
    r.diff = abs(X - target);
    r.ok = r.diff < mpreal("1e-150");
    int n1 = -1, n2 = -1;
    double l1 = 0, l2 = 0;
    for (auto& pr : chk) {
        mpreal err = abs(X - pr.second);
        if (err > 0 && err < mpreal("1e-1") && err > mpreal("1e-180")) {
            double l = (double)log10(err);
            if (n1 < 0) { n1 = pr.first; l1 = l; }
            n2 = pr.first; l2 = l;
        }
    }
    if (n1 >= 0 && n2 > n1) r.digitsPerTerm = (l1 - l2) / (n2 - n1);
    return r;
}

static std::string cstr(const std::vector<int>& c) {
    std::string s = "b=[";
    for (size_t i = 0; i < c.size(); i++) s += (i ? "," : "") + std::to_string(c[i]);
    return s + "] (c0..cK)";
}

static std::string fmtd(const mpreal& d) {
    if (d == 0) return "0";
    char buf[64];
    snprintf(buf, sizeof buf, "%.2f", (double)log10(d));
    return std::string("1e") + buf;
}

static double reltol1 = 1e-11;

// Stage-1 + stage-2 on a coefficient vector. Appends hits.
static void test_candidate(const Family& f, const std::vector<int>& cv, std::vector<Hit>& hits) {
    ld c[8];
    for (int k = 0; k <= f.K; k++) c[k] = cv[k];
    ld x;
    if (!eval_cf(c, f.K, N1, true, x, 1e-14L)) return;
    std::vector<Entry> m;
    lookup(x, reltol1, m);
    if (m.empty()) return;
    ld x2;
    if (!eval_cf(c, f.K, N2, false, x2, 1e-12L)) return;
    for (auto& e : m) {
        ld v = ((ld)e.A + e.B * g_t) / ((ld)e.C + e.D * g_t);
        if (fabsl(x2 - v) <= 1e-14L * (fabsl(v) > 1 ? fabsl(v) : 1)) hits.push_back({cv, e});
    }
}

static bool selftest(const Family& f, const std::vector<std::vector<int>>& cases,
                     const std::vector<std::string>& desc, const mpreal& t, bool required) {
    bool all = true;
    for (size_t i = 0; i < cases.size(); i++) {
        std::vector<Hit> hits;
        test_candidate(f, cases[i], hits);
        std::cout << "  selftest case " << desc[i] << " " << cstr(cases[i]) << ": ";
        if (hits.empty()) { std::cout << "NO HIT\n"; if (required) all = false; continue; }
        for (auto& h : hits) {
            VerifyResult v = verify(f, h.c, h.e, t);
            std::cout << "hit (A,B,C,D)=(" << (int)h.e.A << "," << (int)h.e.B << "," << (int)h.e.C << "," << (int)h.e.D
                      << ") diff=" << fmtd(v.diff) << " digits/term=" << v.digitsPerTerm
                      << (v.ok ? " VERIFIED" : " NOT-VERIFIED") << "; ";
            if (!v.ok && required) all = false;
        }
        std::cout << "\n";
    }
    return all;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: cf_search <F1|F2|F3> <C> <M> <threads> [stride]\n";
        return 1;
    }
    std::string fam = argv[1];
    int C = atoi(argv[2]);
    int M = atoi(argv[3]);
    int nthreads = atoi(argv[4]);
    int stride = argc > 5 ? atoi(argv[5]) : 1;
    if (M > 100) { std::cerr << "M must be <= 100 (int8)\n"; return 1; }

    mpreal::set_default_prec(700);
    std::vector<Family> fams;
    if (fam == "F1") fams.push_back({"F1", 3, 6, -1, 3});
    else if (fam == "F2") { fams.push_back({"F2+", 2, 4, +1, 2}); fams.push_back({"F2-", 2, 4, -1, 2}); }
    else if (fam == "F3") fams.push_back({"F3", 5, 10, -1, 5});
    else { std::cerr << "unknown family\n"; return 1; }

    mpreal tmp = mpfr::zeta(mpreal(fams[0].tz));
    mpreal tmpt = (fams[0].tz == 2) ? mpfr::const_pi() * mpfr::const_pi() / 6 : tmp;
    g_t = (ld)tmpt.toLDouble();
    std::cout << "cf_search family=" << fam << " C=" << C << " M=" << M << " threads=" << nthreads
              << " stride=" << stride << "\n";
    auto T0 = std::chrono::steady_clock::now();
    build_table(M, g_t);
    std::cout << "mobius table size = " << g_table.size() << "\n" << std::flush;

    // self-tests
    bool selfOK = true;
    if (fam == "F1") {
        init_powers(fams[0]);
        std::cout << "SELFTEST F1 (required: Apery 34n^3+51n^2+27n+5 -> 6/zeta(3))\n";
        bool r1 = selftest(fams[0], {{5, 27, 51, 34}}, {"Apery"}, tmpt, true);
        std::cout << "SELFTEST F1 (informational: 2n^3+3n^2+3n+1)\n";
        selftest(fams[0], {{1, 3, 3, 2}}, {"n^3+(n+1)^3"}, tmpt, false);
        selfOK = r1;
    } else if (fam == "F2") {
        std::cout << "SELFTEST F2 (required: 11n^2+11n+3 with a_n = +-n^4 -> 5/zeta(2) type; one sign must hit)\n";
        bool any = false;
        for (auto& ff : fams) {
            init_powers(ff);
            std::cout << " family " << ff.name << "\n";
            any |= selftest(ff, {{3, 11, 11}}, {"Apery2"}, tmpt, true);
        }
        selfOK = any;
    } else {
        init_powers(fams[0]);
        std::cout << "SELFTEST F3: none known (skipped)\n";
    }
    std::cout << (selfOK ? "SELFTEST PASS" : "SELFTEST FAIL") << "\n" << std::flush;

    // main search
    for (auto& f : fams) {
        init_powers(f);
        int W = 2 * C + 1;
        long long ntasks = (long long)W * W;   // top two coefficients
        std::atomic<long long> next(0);
        std::vector<Hit> allHits;
        std::mutex mu;
        std::atomic<long long> ncand(0), nconv(0);
        auto worker = [&]() {
            std::vector<Hit> hits;
            long long local = 0;
            while (true) {
                long long id = next.fetch_add(stride);
                if (id >= ntasks) break;
                std::vector<int> cv(f.K + 1, -C);
                cv[f.K] = (int)(id / W) - C;
                cv[f.K - 1] = (int)(id % W) - C;
                // odometer over c_0..c_{K-2}
                while (true) {
                    local++;
                    test_candidate(f, cv, hits);
                    int k = 0;
                    while (k <= f.K - 2) {
                        if (cv[k] < C) { cv[k]++; break; }
                        cv[k] = -C; k++;
                    }
                    if (k > f.K - 2) break;
                }
            }
            std::lock_guard<std::mutex> lk(mu);
            allHits.insert(allHits.end(), hits.begin(), hits.end());
            ncand += local;
        };
        std::cout << "\n=== search " << f.name << ": a_n=" << (f.sa > 0 ? "+" : "-") << "n^" << f.ea
                  << ", deg b=" << f.K << ", tasks=" << ntasks << " ===\n" << std::flush;
        auto T1 = std::chrono::steady_clock::now();
        std::vector<std::thread> th;
        for (int i = 0; i < nthreads; i++) th.emplace_back(worker);
        for (auto& t : th) t.join();
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - T1).count();
        std::cout << "candidates tested = " << ncand << "  search time = " << sec << " s ("
                  << (sec * 1e6 * nthreads / std::max<long long>(1, ncand)) << " us/candidate/thread)\n";
        std::cout << "stage-2 hits = " << allHits.size() << "\n";
        std::sort(allHits.begin(), allHits.end(), [](const Hit& a, const Hit& b) { return a.c < b.c; });
        mpreal t = (f.tz == 2) ? mpfr::const_pi() * mpfr::const_pi() / 6 : mpfr::zeta(mpreal(f.tz));
        std::map<std::tuple<int, int, int, int>, int> cnt;
        int nver = 0;
        for (auto& h : allHits) {
            VerifyResult v = verify(f, h.c, h.e, t);
            std::cout << f.name << " " << cstr(h.c) << " (A,B,C,D)=(" << (int)h.e.A << "," << (int)h.e.B << ","
                      << (int)h.e.C << "," << (int)h.e.D << ") diff=" << fmtd(v.diff)
                      << " digits/term=" << v.digitsPerTerm << (v.ok ? " VERIFIED" : " rejected") << "\n";
            if (v.ok) {
                nver++;
                cnt[{h.e.A, h.e.B, h.e.C, h.e.D}]++;
            }
        }
        std::cout << "SUMMARY " << f.name << ": stage-2 hits=" << allHits.size() << " verified=" << nver << "\n";
        for (auto& kv : cnt)
            std::cout << "  (A,B,C,D)=(" << std::get<0>(kv.first) << "," << std::get<1>(kv.first) << ","
                      << std::get<2>(kv.first) << "," << std::get<3>(kv.first) << ") x" << kv.second << "\n";
    }
    double tot = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    std::cout << "elapsed " << tot << " s\n";
    return 0;
}
