// Ramanujan-Machine style search for polynomial continued fractions
//   x = b_0 + a_1/(b_1 + a_2/(b_2 + ...)),   b_n = sum_k c_k n^k,  b_0 = c_0
// matched against Mobius transforms (A + B t)/(C + D t) of one or many target constants t.
//
// Usage:
//   cf_search <F1|F2|F3>[multi] <C> <M> <threads> [stride]
//   cf_search F4[multi] <C> <M> <threads> <Aalpha> <Ab> [stride]
//   cf_search G1 <Ca> <Cb> <M> <threads> [stride]
//
//   F1: t = zeta(3), a_n = -n^6,                  b cubic,     |c_i| <= C
//   F2: t = zeta(2), a_n = +-n^4,                 b quadratic
//   F3: t = zeta(5), a_n = -n^10,                 b quintic
//   F4: t = zeta(3), a_n = -(al n^2 + be n + ga)^3, 1 <= al <= Aalpha, |be|,|ga| <= Ab, b cubic
//   G1: general a_n = a2 n^2 + a1 n + a0 (|.| <= Ca), b_n quadratic (|.| <= Cb); always multi-target.
//       Polynomially convergent CFs are handled by Richardson extrapolation (error series in 1/n).
//   "multi": instead of the single target, use f(zeta(3) pi^k) for 17 elementary functions f,
//            k in {-3,-2,-1,0}, plus the controls zeta3, e, pi.
//   stride (optional): run only every stride-th task (timing estimates only; not a full search)
// Build: g++ -O2 -std=c++17 cf_search.cpp -o cf_search -lmpfr -lgmp -pthread

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
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
    int ea;         // exponent in a_n (F1-F3)
    int sa;         // sign of a_n (F1-F3)
    int tz;         // target zeta(tz) in single-target mode
    int poly = 0;   // 1: F4 a_n = -(al n^2+be n+ga)^3 ; 2: G1 a_n = al n^2 + be n + ga
    int al = 0, be = 0, ga = 0;
    bool accel = false;   // G1: allow polynomial convergence via Richardson extrapolation
};

struct Entry {
    double val;
    int8_t A, B, C, D;
    uint8_t tg;
};

static const int N1 = 200;   // terms, stage 1 (non-accel)
static const int N2 = 600;   // terms, stage 2 (non-accel)
static const int NMAXACC = 8192;   // accel stage 2: n0=1024, 4 levels

static int g_nmax = N2;
static std::vector<ld> g_pw[8];   // g_pw[k][n] = n^k
static std::vector<Entry> g_table;
static std::vector<std::string> g_tname;
static std::vector<ld> g_tl;
static std::vector<mpreal> g_tm;

static void init_powers(const Family& f) {
    g_nmax = f.accel ? NMAXACC : N2;
    for (int k = 0; k <= f.K; k++) {
        g_pw[k].assign(g_nmax + 1, 0);
        for (int n = 0; n <= g_nmax; n++) g_pw[k][n] = powl((ld)n, (ld)k);
    }
}

// a_n for n = 0..g_nmax; returns false if some a_n (n>=1) is zero
static bool make_a(const Family& f, std::vector<ld>& a) {
    a.assign(g_nmax + 1, 0);
    for (int n = 0; n <= g_nmax; n++) {
        ld base = (ld)f.al * n * n + (ld)f.be * n + f.ga;
        if (f.poly == 1) {
            if (n >= 1 && base == 0) return false;
            a[n] = -base * base * base;
        } else if (f.poly == 2) {
            if (n >= 1 && base == 0) return false;
            a[n] = base;
        } else a[n] = f.sa * powl((ld)n, (ld)f.ea);
    }
    return true;
}

// Evaluate the CF with N terms. Returns false if rejected. Early-prunes slow convergence when prune=true.
static bool eval_cf(const ld* g_a, const ld* c, int K, int N, bool prune, ld& out, ld tol) {
    ld pm1 = 1, qm1 = 0;             // p_{-1}, q_{-1}
    ld p0 = c[0], q0 = 1;            // p_0, q_0
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
    return true;
}

// Convergents at n0, 2 n0, 4 n0, 8 n0.
static bool eval_levels(const ld* av, const ld* c, int K, int n0, ld* r) {
    ld pm1 = 1, qm1 = 0, p0 = c[0], q0 = 1;
    int idx = 0, target = n0;
    const int Nmax = n0 * 8;
    for (int n = 1; n <= Nmax; n++) {
        ld b = 0;
        for (int k = 0; k <= K; k++) b += c[k] * g_pw[k][n];
        ld a = av[n];
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
        if (n == target) {
            if (!(fabsl(q0) > 1e-300L)) return false;
            r[idx] = p0 / q0;
            if (!std::isfinite((double)r[idx])) return false;
            idx++;
            target *= 2;
        }
    }
    return true;
}

// From 4 convergents at n, 2n, 4n, 8n. Geometric: x = r3, tol = geomTol.
// Polynomial (error ~ sum c_k/n^k): Richardson with exponents 1,2,3; tol from the last two estimates.
static bool accel_value(const ld* r, ld geomTol, ld maxPolyTol, ld& x, ld& tol, bool& geom) {
    ld scale = fabsl(r[3]) > 1 ? fabsl(r[3]) : 1;
    ld d1 = fabsl(r[1] - r[0]), d2 = fabsl(r[2] - r[1]), d3 = fabsl(r[3] - r[2]);
    if (d3 <= 1e-14L * scale) { x = r[3]; tol = geomTol; geom = true; return true; }
    geom = false;
    if (!(d3 < 0.75L * d2 && d2 < 0.9L * d1 + 1e-300L)) return false;
    ld t1[3], t2[2];
    for (int i = 0; i < 3; i++) t1[i] = 2 * r[i + 1] - r[i];
    for (int i = 0; i < 2; i++) t2[i] = (4 * t1[i + 1] - t1[i]) / 3;
    ld t3 = (8 * t2[1] - t2[0]) / 7;
    x = t3;
    tol = 4 * fabsl(t3 - t2[1]) / scale;
    if (tol < 1e-12L) tol = 1e-12L;
    if (!(tol <= maxPolyTol)) return false;
    return true;
}

static long long nDropped = 0;

static void build_table(int M) {
    g_table.clear();
    for (size_t tg = 0; tg < g_tl.size(); tg++) {
        ld t = g_tl[tg];
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
                        if (!std::isfinite((double)v) || fabsl(v) > 1e150L) continue;
                        // Drop images indistinguishable from the rationals B/D (|t| huge) or A/C (|t| tiny):
                        // such entries would match any rational-valued CF.
                        ld sc = fabsl(v) > 1 ? fabsl(v) : 1;
                        if (D != 0 && fabsl(v - (ld)B / D) < 1e-6L * sc) { nDropped++; continue; }
                        if (C != 0 && fabsl(v - (ld)A / C) < 1e-6L * sc) { nDropped++; continue; }
                        g_table.push_back({(double)v, (int8_t)A, (int8_t)B, (int8_t)C, (int8_t)D, (uint8_t)tg});
                    }
    }
    std::sort(g_table.begin(), g_table.end(), [](const Entry& x, const Entry& y) { return x.val < y.val; });
}

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
    int al = 0, be = 0, ga = 0;
};

struct VerifyResult {
    mpreal diff;
    double digitsPerTerm;
    bool ok;
    bool accelUsed;
};

static std::vector<mpreal> mp_ratios(const Family& f, const std::vector<int>& c, const std::vector<int>& ns) {
    mpreal pm1 = 1, qm1 = 0, p0 = c[0], q0 = 1;
    std::vector<mpreal> out;
    size_t idx = 0;
    int NN = ns.back();
    for (int n = 1; n <= NN; n++) {
        long long bn = 0, pw = 1;
        for (int k = 0; k <= f.K; k++) { bn += (long long)c[k] * pw; pw *= n; }
        mpreal a;
        if (f.poly == 1) { mpreal base = mpreal(f.al) * n * n + mpreal(f.be) * n + f.ga; a = -base * base * base; }
        else if (f.poly == 2) { a = mpreal(f.al) * n * n + mpreal(f.be) * n + f.ga; }
        else a = mpreal(f.sa) * pow(mpreal(n), f.ea);
        mpreal p = bn * p0 + a * pm1;
        mpreal q = bn * q0 + a * qm1;
        pm1 = p0; qm1 = q0; p0 = p; q0 = q;
        if (n == ns[idx]) {
            out.push_back(q0 == 0 ? mpreal(0) : mpreal(p0 / q0));
            idx++;
        }
    }
    return out;
}

static VerifyResult verify(const Family& f, const std::vector<int>& c, const Entry& e) {
    const mpreal& t = g_tm[e.tg];
    std::vector<int> ns;
    for (int n = 10; n <= 200; n += 10) ns.push_back(n);
    ns.push_back(2000);
    std::vector<mpreal> rs = mp_ratios(f, c, ns);
    VerifyResult r;
    r.digitsPerTerm = NAN;
    r.accelUsed = false;
    mpreal X = rs.back();
    mpreal target = (mpreal(e.A) + e.B * t) / (mpreal(e.C) + e.D * t);
    r.diff = abs(X - target);
    r.ok = r.diff < mpreal("1e-150");
    int n1 = -1, n2 = -1;
    double l1 = 0, l2 = 0;
    for (size_t i = 0; i + 1 < rs.size(); i++) {
        mpreal err = abs(X - rs[i]);
        if (err > 0 && err < mpreal("1e-1") && err > mpreal("1e-180")) {
            double l = (double)log10(err);
            if (n1 < 0) { n1 = ns[i]; l1 = l; }
            n2 = ns[i];
            l2 = l;
        }
    }
    if (n1 >= 0 && n2 > n1) r.digitsPerTerm = (l1 - l2) / (n2 - n1);
    if (!r.ok && f.accel) {
        // polynomial convergence: Neville/Richardson with exponents 1..7 on n = 250 * 2^j, j = 0..7
        std::vector<int> ns2;
        for (int j = 0; j < 8; j++) ns2.push_back(250 << j);
        std::vector<mpreal> T = mp_ratios(f, c, ns2);
        for (int k = 1; k < 8; k++) {
            mpreal w = pow(mpreal(2), k);
            for (int i = 0; i + k < 8; i++) T[i] = (w * T[i + 1] - T[i]) / (w - 1);
        }
        mpreal d2 = abs(T[0] - target);
        r.accelUsed = true;
        r.diff = d2;
        r.ok = d2 < mpreal("1e-9");
        r.digitsPerTerm = NAN;
    }
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
static void test_candidate(const ld* av, const Family& f, const std::vector<int>& cv, std::vector<Hit>& hits) {
    ld c[8];
    for (int k = 0; k <= f.K; k++) c[k] = cv[k];
    if (!f.accel) {
        ld x;
        if (!eval_cf(av, c, f.K, N1, true, x, 1e-14L)) return;
        std::vector<Entry> m;
        lookup(x, reltol1, m);
        if (m.empty()) return;
        ld x2;
        if (!eval_cf(av, c, f.K, N2, false, x2, 1e-12L)) return;
        for (auto& e : m) {
            ld t = g_tl[e.tg];
            ld v = ((ld)e.A + e.B * t) / ((ld)e.C + e.D * t);
            if (fabsl(x2 - v) <= 1e-14L * (fabsl(v) > 1 ? fabsl(v) : 1)) hits.push_back({cv, e});
        }
        return;
    }
    ld r[4], x, tol;
    bool geom;
    if (!eval_levels(av, c, f.K, 128, r)) return;
    if (!accel_value(r, 1e-11L, 1e-5L, x, tol, geom)) return;
    std::vector<Entry> m;
    lookup(x, (double)tol, m);
    if (m.empty()) return;
    ld r2[4], x2, tol2;
    bool geom2;
    if (!eval_levels(av, c, f.K, 1024, r2)) return;
    if (!accel_value(r2, 1e-13L, 1e-9L, x2, tol2, geom2)) return;
    for (auto& e : m) {
        ld t = g_tl[e.tg];
        ld v = ((ld)e.A + e.B * t) / ((ld)e.C + e.D * t);
        ld sc = fabsl(v) > 1 ? fabsl(v) : 1;
        if (fabsl(x2 - v) <= tol2 * sc) hits.push_back({cv, e});
    }
}

// Passes if some VERIFIED hit has target name expectTarget (or, if expectTarget is empty, nothing is required).
static bool selftest(const Family& f, const std::vector<int>& c, const std::string& desc,
                     const std::string& expectTarget) {
    std::vector<Hit> hits;
    std::vector<ld> av;
    make_a(f, av);
    test_candidate(av.data(), f, c, hits);
    std::cout << "  selftest case " << desc << " " << cstr(c) << ": ";
    bool found = false;
    if (hits.empty()) std::cout << "NO HIT";
    for (auto& h : hits) {
        VerifyResult v = verify(f, c, h.e);
        std::cout << "[" << g_tname[h.e.tg] << "] (A,B,C,D)=(" << (int)h.e.A << "," << (int)h.e.B << ","
                  << (int)h.e.C << "," << (int)h.e.D << ") diff=" << fmtd(v.diff)
                  << " digits/term=" << v.digitsPerTerm
                  << (v.ok ? (v.accelUsed ? " VERIFIED-ACCEL" : " VERIFIED") : " NOT-VERIFIED") << "; ";
        if (v.ok && g_tname[h.e.tg] == expectTarget) found = true;
    }
    std::cout << "\n";
    return expectTarget.empty() ? true : found;
}

static void add_target(const std::string& name, const mpreal& v) {
    if (mpfr::isnan(v) || mpfr::isinf(v)) {
        std::cout << "skip target " << name << " (outside the domain)\n";
        return;
    }
    ld l = (ld)v.toLDouble();
    if (!std::isfinite((double)l) || fabsl(l) > 1e100L) {
        std::cout << "skip target " << name << " (too large)\n";
        return;
    }
    g_tname.push_back(name);
    g_tl.push_back(l);
    g_tm.push_back(v);
}

static void build_multi_targets() {
    struct Fn { std::string name; std::function<mpreal(const mpreal&)> f; };
    const std::vector<Fn> fns = {
        {"exp", [](const mpreal& t) { return mpfr::exp(t); }},
        {"log", [](const mpreal& t) { return mpfr::log(t); }},
        {"sin", [](const mpreal& t) { return mpfr::sin(t); }},
        {"cos", [](const mpreal& t) { return mpfr::cos(t); }},
        {"tan", [](const mpreal& t) { return mpfr::tan(t); }},
        {"sinh", [](const mpreal& t) { return mpfr::sinh(t); }},
        {"cosh", [](const mpreal& t) { return mpfr::cosh(t); }},
        {"tanh", [](const mpreal& t) { return mpfr::tanh(t); }},
        {"arctan", [](const mpreal& t) { return mpfr::atan(t); }},
        {"arcsinh", [](const mpreal& t) { return mpfr::asinh(t); }},
        {"arccosh", [](const mpreal& t) { return mpfr::acosh(t); }},
        {"arcsin(1/t)", [](const mpreal& t) { return mpfr::asin(1 / t); }},
        {"arccos(1/t)", [](const mpreal& t) { return mpfr::acos(1 / t); }},
        {"log(t-1)", [](const mpreal& t) { return mpfr::log(t - 1); }},
        {"t^t", [](const mpreal& t) { return mpfr::pow(t, t); }},
        {"Gamma", [](const mpreal& t) { return mpfr::gamma(t); }},
        {"exp(1/t)", [](const mpreal& t) { return mpfr::exp(1 / t); }},
    };
    mpreal z3 = mpfr::zeta(mpreal(3)), pi = mpfr::const_pi();
    for (auto& fn : fns)
        for (int k : {-3, -2, -1, 0}) {
            mpreal t = z3 * pow(pi, mpreal(k));
            add_target(fn.name + "(zeta3*pi^" + std::to_string(k) + ")", fn.f(t));
        }
    add_target("zeta3", z3);
    add_target("e", mpfr::exp(mpreal(1)));
    add_target("pi", pi);
}

int main(int argc, char** argv) {
    auto usage = []() {
        std::cerr << "Usage: cf_search <F1|F2|F3>[multi] <C> <M> <threads> [stride]\n"
                  << "       cf_search F4[multi] <C> <M> <threads> <Aalpha> <Ab> [stride]\n"
                  << "       cf_search G1 <Ca> <Cb> <M> <threads> [stride]\n";
    };
    if (argc < 5) { usage(); return 1; }
    std::string fam = argv[1];
    bool multi = false;
    if (fam.size() > 5 && fam.substr(fam.size() - 5) == "multi") { multi = true; fam = fam.substr(0, fam.size() - 5); }
    int C = 0, Ca = 0, M = 0, nthreads = 0, Aal = 0, Ab = 0, stride = 1;
    if (fam == "G1") {
        if (argc < 6) { usage(); return 1; }
        multi = true;
        Ca = atoi(argv[2]); C = atoi(argv[3]); M = atoi(argv[4]); nthreads = atoi(argv[5]);
        stride = argc > 6 ? atoi(argv[6]) : 1;
    } else if (fam == "F4") {
        if (argc < 7) { usage(); return 1; }
        C = atoi(argv[2]); M = atoi(argv[3]); nthreads = atoi(argv[4]);
        Aal = atoi(argv[5]); Ab = atoi(argv[6]);
        stride = argc > 7 ? atoi(argv[7]) : 1;
    } else {
        C = atoi(argv[2]); M = atoi(argv[3]); nthreads = atoi(argv[4]);
        stride = argc > 5 ? atoi(argv[5]) : 1;
    }
    if (M > 100) { std::cerr << "M must be <= 100 (int8)\n"; return 1; }

    mpreal::set_default_prec(700);
    std::vector<Family> fams;
    if (fam == "F1") fams.push_back({"F1", 3, 6, -1, 3});
    else if (fam == "F2") { fams.push_back({"F2+", 2, 4, +1, 2}); fams.push_back({"F2-", 2, 4, -1, 2}); }
    else if (fam == "F4") { Family f4{"F4", 3, 0, -1, 3}; f4.poly = 1; fams.push_back(f4); }
    else if (fam == "F3") fams.push_back({"F3", 5, 10, -1, 5});
    else if (fam == "G1") { Family g{"G1", 2, 0, 1, 3}; g.poly = 2; g.accel = true; fams.push_back(g); }
    else { usage(); return 1; }

    std::cout << "cf_search family=" << fam << (multi ? " (multi-target)" : "") << " C=" << C
              << (fam == "G1" ? " Ca=" + std::to_string(Ca) : std::string()) << " M=" << M << " threads=" << nthreads
              << " stride=" << stride << "\n";
    auto T0 = std::chrono::steady_clock::now();
    if (multi) build_multi_targets();
    else {
        int tz = fams[0].tz;
        mpreal t = (tz == 2) ? mpfr::const_pi() * mpfr::const_pi() / 6 : mpfr::zeta(mpreal(tz));
        add_target("zeta" + std::to_string(tz), t);
    }
    std::cout << "targets = " << g_tname.size() << "\n";
    build_table(M);
    std::cout << "mobius table size = " << g_table.size() << " entries ("
              << (g_table.size() * sizeof(Entry) / 1048576) << " MB); dropped near-rational images: " << nDropped << "\n" << std::flush;

    // self-tests
    bool selfOK = true;
    if (fam == "F1") {
        init_powers(fams[0]);
        std::cout << "SELFTEST F1 (required: Apery 34n^3+51n^2+27n+5 -> 6/zeta(3))\n";
        selfOK = selftest(fams[0], {5, 27, 51, 34}, "Apery", "zeta3");
        std::cout << "SELFTEST F1 (informational: 2n^3+3n^2+3n+1)\n";
        selftest(fams[0], {1, 3, 3, 2}, "n^3+(n+1)^3", "");
    } else if (fam == "F2") {
        std::cout << "SELFTEST F2 (required: 11n^2+11n+3 with a_n = +-n^4 -> 5/zeta(2) type; one sign must hit)\n";
        bool any = false;
        for (auto& ff : fams) {
            init_powers(ff);
            std::cout << " family " << ff.name << "\n";
            any |= selftest(ff, {3, 11, 11}, "Apery2", "zeta2");
        }
        selfOK = any;
    } else if (fam == "F4") {
        init_powers(fams[0]);
        Family fa = fams[0], fb = fams[0];
        fa.al = 1; fa.be = 0; fa.ga = 0;
        fb.al = 1; fb.be = 2; fb.ga = 1;
        std::cout << "SELFTEST F4 (required: (al,be,ga)=(1,0,0) Apery -> 6/zeta(3);\n"
                  << "  shifted Apery (1,2,1), b'_n=b_{n+1}=[117,231,153,34] -> zeta3/(5 zeta3-6) = (0,1,-6,5))\n";
        bool r1 = selftest(fa, {5, 27, 51, 34}, "Apery(1,0,0)", "zeta3");
        bool r2 = selftest(fb, {117, 231, 153, 34}, "shifted(1,2,1)", "zeta3");
        selfOK = r1 && r2;
    } else if (fam == "G1") {
        init_powers(fams[0]);
        Family fp = fams[0], fe = fams[0];
        fp.al = 1; fp.be = 0; fp.ga = 0;     // a_n = n^2
        fe.al = 0; fe.be = 1; fe.ga = 0;     // a_n = n
        std::cout << "SELFTEST G1 (required: a_n=n^2, b_n=2n+1 -> 4/pi (target pi, polynomial convergence);\n"
                  << "  a_n=n, b_n=n+1 -> 1/(e-2) (target e); a_n=(2n+1)^2, b_n=2 -> pi/(4-pi) via Richardson)\n";
        Family fq = fams[0];
        fq.al = 4; fq.be = 4; fq.ga = 1;     // a_n = (2n+1)^2, b_n = 2: pi/(4-pi) (polynomial convergence, accel path)
        bool r1 = selftest(fp, {1, 2, 0}, "4/pi (geometric, 0.77 digits/term)", "pi");
        bool r2 = selftest(fe, {1, 1, 0}, "1/(e-2)", "e");
        bool r3 = selftest(fq, {2, 0, 0}, "Brouncker tail pi/(4-pi), a_n=(2n+1)^2, b_n=2 (polynomial conv.)", "pi");
        selfOK = r1 && r2 && r3;
    } else {
        init_powers(fams[0]);
        std::cout << "SELFTEST F3: none known (skipped)\n";
    }
    if (fam == "F1" && multi) {
        // already covered: expected target zeta3
    }
    std::cout << (selfOK ? "SELFTEST PASS" : "SELFTEST FAIL") << "\n" << std::flush;

    // main search
    for (auto& f : fams) {
        init_powers(f);
        int W = 2 * C + 1;
        const bool isPoly = f.poly != 0;
        const int na = 2 * Ca + 1, nb = 2 * Ab + 1;
        long long nabg = 1;
        if (f.poly == 1) nabg = (long long)Aal * nb * nb;
        if (f.poly == 2) nabg = (long long)na * na * na;
        const int nfix = isPoly ? 1 : 2;   // number of fixed top coefficients per task
        long long ntasks = isPoly ? nabg * W : (long long)W * W;
        std::atomic<long long> next(0);
        std::vector<Hit> allHits;
        std::mutex mu;
        std::atomic<long long> ncand(0);
        auto worker = [&]() {
            std::vector<Hit> hits, tmp;
            std::vector<ld> av;
            long long local = 0;
            while (true) {
                long long id = next.fetch_add(stride);
                if (id >= ntasks) break;
                std::vector<int> cv(f.K + 1, -C);
                Family ff = f;
                if (isPoly) {
                    long long g = id / W;
                    cv[f.K] = (int)(id % W) - C;
                    if (f.poly == 1) {
                        ff.ga = (int)(g % nb) - Ab; g /= nb;
                        ff.be = (int)(g % nb) - Ab; g /= nb;
                        ff.al = (int)g + 1;
                    } else {
                        ff.ga = (int)(g % na) - Ca; g /= na;
                        ff.be = (int)(g % na) - Ca; g /= na;
                        ff.al = (int)g - Ca;
                    }
                    if (!make_a(ff, av)) continue;
                } else {
                    cv[f.K] = (int)(id / W) - C;
                    cv[f.K - 1] = (int)(id % W) - C;
                    make_a(ff, av);
                }
                while (true) {
                    local++;
                    tmp.clear();
                    test_candidate(av.data(), ff, cv, tmp);
                    for (auto& h : tmp) { h.al = ff.al; h.be = ff.be; h.ga = ff.ga; hits.push_back(h); }
                    int k = 0;
                    while (k <= f.K - nfix) {
                        if (cv[k] < C) { cv[k]++; break; }
                        cv[k] = -C; k++;
                    }
                    if (k > f.K - nfix) break;
                }
            }
            std::lock_guard<std::mutex> lk(mu);
            allHits.insert(allHits.end(), hits.begin(), hits.end());
            ncand += local;
        };
        std::cout << "\n=== search " << f.name << ": a_n=";
        if (f.poly == 1) std::cout << "-(al n^2+be n+ga)^3, al<=" << Aal << ", |be|,|ga|<=" << Ab;
        else if (f.poly == 2) std::cout << "a2 n^2+a1 n+a0, |a_i|<=" << Ca;
        else std::cout << (f.sa > 0 ? "+" : "-") << "n^" << f.ea;
        std::cout << ", deg b=" << f.K << ", |c_i|<=" << C << ", tasks=" << ntasks << " ===\n" << std::flush;
        auto T1 = std::chrono::steady_clock::now();
        std::vector<std::thread> th;
        for (int i = 0; i < nthreads; i++) th.emplace_back(worker);
        for (auto& t : th) t.join();
        double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - T1).count();
        std::cout << "candidates tested = " << ncand << "  search time = " << sec << " s ("
                  << (sec * 1e6 * nthreads / std::max<long long>(1, ncand)) << " us/candidate/thread)\n";
        std::cout << "stage-2 hits = " << allHits.size() << "\n";
        std::sort(allHits.begin(), allHits.end(), [](const Hit& a, const Hit& b) {
            return std::tie(a.al, a.be, a.ga, a.c) < std::tie(b.al, b.be, b.ga, b.c);
        });
        std::map<std::string, std::pair<int, int>> perTarget;   // stage2, verified
        std::map<std::tuple<std::string, int, int, int, int>, int> cnt;
        int nver = 0;
        for (auto& h : allHits) {
            Family ff = f; ff.al = h.al; ff.be = h.be; ff.ga = h.ga;
            VerifyResult v = verify(ff, h.c, h.e);
            const std::string& tn = g_tname[h.e.tg];
            std::cout << f.name;
            if (f.poly) std::cout << " (a:" << h.al << "," << h.be << "," << h.ga << ")";
            std::cout << " [" << tn << "] " << cstr(h.c) << " (A,B,C,D)=(" << (int)h.e.A << "," << (int)h.e.B << ","
                      << (int)h.e.C << "," << (int)h.e.D << ") diff=" << fmtd(v.diff)
                      << " digits/term=" << v.digitsPerTerm
                      << (v.ok ? (v.accelUsed ? " VERIFIED-ACCEL" : " VERIFIED") : " rejected") << "\n";
            perTarget[tn].first++;
            if (v.ok) {
                nver++;
                perTarget[tn].second++;
                cnt[{tn, h.e.A, h.e.B, h.e.C, h.e.D}]++;
            }
        }
        std::cout << "SUMMARY " << f.name << ": stage-2 hits=" << allHits.size() << " verified=" << nver << "\n";
        for (auto& kv : cnt)
            std::cout << "  [" << std::get<0>(kv.first) << "] (A,B,C,D)=(" << std::get<1>(kv.first) << ","
                      << std::get<2>(kv.first) << "," << std::get<3>(kv.first) << "," << std::get<4>(kv.first)
                      << ") x" << kv.second << "\n";
        std::cout << "PER-TARGET (stage-2 hits / verified):\n";
        for (auto& tn : g_tname) {
            auto it = perTarget.find(tn);
            std::cout << "  " << tn << ": " << (it == perTarget.end() ? 0 : it->second.first) << " / "
                      << (it == perTarget.end() ? 0 : it->second.second) << "\n";
        }
    }
    double tot = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    std::cout << "elapsed " << tot << " s\n";
    return 0;
}
