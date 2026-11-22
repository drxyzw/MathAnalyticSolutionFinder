// Standalone validation / benchmark driver for pslqStandard vs pslqMultipair2L.
// Build: g++ -O2 -std=c++17 pslq_m2l_check.cpp pslq.cpp -o pslq_m2l_check -lmpfr -lgmp -pthread
// Usage: ./pslq_m2l_check known | nore <digits> <n> | planted <digits> <n> <log10coef> <trials> | bench <digits> <n>
#include "pslq.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>

static double now()
{
	return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static std::string relStr(const std::vector<mpreal>& r)
{
	std::string s = "[";
	for (size_t i = 0; i < r.size(); ++i)
		s += (i ? "," : "") + r[i].toString("%.0Rf");
	return s + "]";
}

// normalize sign so that the first nonzero coefficient is positive
static std::vector<mpreal> canon(std::vector<mpreal> r)
{
	for (auto& v : r) { if (v != 0) { if (v < 0) for (auto& w : r) w = -w; break; } }
	return r;
}

struct Run { PslqResult r; double t; };
static Run run(bool m2l, const std::vector<mpreal>& v)
{
	const double t0 = now();
	Run o;
	o.r = m2l ? pslqMultipair2L(v) : pslqStandard(v);
	o.t = now() - t0;
	return o;
}

static void line(const char* name, const std::vector<mpreal>& v)
{
	Run a = run(false, v), b = run(true, v);
	auto fmt = [](const Run& o) {
		std::string s = toString(o.r.status);
		s += " nb=" + o.r.normBound.toString(3) + " it=" + std::to_string(o.r.iterations) + " t=" + std::to_string(o.t) + "s";
		if (o.r.status == PslqResult::Status::Found) s += " rel=" + relStr(canon(o.r.relation));
		return s; };
	for (const Run* o : { &a, &b })
		if (o->r.status == PslqResult::Status::Found)
		{
			mpreal n2 = 0; for (auto& v : o->r.relation) n2 += v * v;
			if (o->r.normBound > mpfr::sqrt(n2) * (1 + mpreal(1e-30))) printf("   !!! normBound exceeds the found relation's norm\n");
		}
	bool same = a.r.status == b.r.status;
	if (same && a.r.status == PslqResult::Status::Found) {
		auto ra = canon(a.r.relation), rb = canon(b.r.relation);
		same = ra == rb;
	}
	printf("%-28s\n   std: %s\n   m2l: %s\n   %s\n", name, fmt(a).c_str(), fmt(b).c_str(), same ? "SAME" : "DIFFERENT");
	fflush(stdout);
}

static mpreal li3half(int digits)
{
	mpreal s = 0, p = 1;
	for (int k = 1; k < 3.33 * digits + 20; ++k) { p /= 2; s += p / mpreal(k) / mpreal(k) / mpreal(k); }
	return s;
}

static mpreal polyRoot(const std::vector<int>& c) // real root by scan + Newton
{
	auto f = [&](const mpreal& x) { mpreal v = 0; for (size_t i = c.size(); i-- > 0;) v = v * x + c[i]; return v; };
	auto df = [&](const mpreal& x) { mpreal v = 0; for (size_t i = c.size(); i-- > 1;) v = v * x + mpreal(c[i]) * mpreal((int)i); return v; };
	mpreal x0 = -6, prev = f(x0);
	for (double xx = -5.9; xx <= 6; xx += 0.05)
	{
		mpreal v = f(mpreal(xx));
		if ((v > 0) != (prev > 0)) { x0 = xx; break; }
		prev = v;
	}
	mpreal x = x0;
	for (int i = 0; i < 400; ++i) x -= f(x) / df(x);
	return x;
}

static std::vector<mpreal> powers(const mpreal& x, int deg)
{
	std::vector<mpreal> p(deg + 1); p[0] = 1;
	for (int k = 1; k <= deg; ++k) p[k] = p[k - 1] * x;
	return p;
}

static mpreal randomReal(std::mt19937_64& g, int digits)
{
	mpreal v = 0, sc = 1;
	for (int k = 0; k < digits / 9.0 + 3; ++k) { sc /= 1e9; v += mpreal((double)(g() % 1000000000ULL)) * sc; }
	return v + mpreal((double)(g() % 1000)) / 1000 + 0.5;
}

int main(int argc, char** argv)
{
	const std::string mode = argc > 1 ? argv[1] : "known";
	if (mode == "known")
	{
		const int D = 800;
		mpreal::set_default_prec(mpfr::digits2bits(D));
		const mpreal pi = mpfr::const_pi(), l2 = mpfr::log(mpreal(2)), l3 = mpfr::log(mpreal(3)), z3 = mpfr::zeta(mpreal(3));
		line("Li3(1/2) identity", { li3half(D), z3, pi * pi * l2, l2 * l2 * l2 });
		line("2^(1/4)+3^(1/4) deg16", powers(mpfr::root(mpreal(2), 4) + mpfr::root(mpreal(3), 4), 16));
		line("cos(2pi/51) deg16", powers(mpfr::cos(2 * pi / 51), 16));
		line("2^(1/12) Eisenstein deg12", powers(mpfr::root(mpreal(2), 12), 12));
		line("Eisenstein(2) random deg12", powers(polyRoot({ 2, 4, 6, -8, 10, 12, -14, 6, -4, 8, 2, -6, 1 }), 12));
		line("log2,log3,log6", { l2, l3, mpfr::log(mpreal(6)) });
		line("log 2,3,5,6,10,15", { l2, l3, mpfr::log(mpreal(5)), mpfr::log(mpreal(6)), mpfr::log(mpreal(10)), mpfr::log(mpreal(15)) });
	}
	else if (mode == "nore")
	{
		const int D = atoi(argv[2]); const int n = atoi(argv[3]);
		mpreal::set_default_prec(mpfr::digits2bits(D));
		const mpreal pi = mpfr::const_pi(), z3 = mpfr::zeta(mpreal(3));
		const int maxT = argc > 4 ? atoi(argv[4]) : 99;
		int cnt = 0;
		for (int k = -3; k <= 0; ++k)
		{
			if (cnt++ >= maxT) break;
			mpreal t = z3 * mpfr::pow(pi, mpreal(k));
			line(("zeta3*pi^" + std::to_string(k)).c_str(), powers(t, n - 1));
		}
		if (cnt++ < maxT) line("zeta5", powers(mpfr::zeta(mpreal(5)), n - 1));
		if (cnt++ < maxT) line("exp(zeta3)", powers(mpfr::exp(z3), n - 1));
	}
	else if (mode == "planted")
	{
		const int D = atoi(argv[2]); const int n = atoi(argv[3]); const double lc = atof(argv[4]); const int trials = atoi(argv[5]);
		mpreal::set_default_prec(mpfr::digits2bits(D));
		std::mt19937_64 g(12345 + n * 1000 + D);
		int foundS = 0, foundM = 0, same = 0; double tS = 0, tM = 0;
		for (int tr = 0; tr < trials; ++tr)
		{
			std::vector<mpreal> c(n), x(n);
			for (int i = 0; i < n; ++i)
			{
				double sc = std::pow(10.0, lc) * (0.5 + (double)(g() % 1000) / 1000.0);
				mpreal v = mpfr::floor(randomReal(g, 40) * sc) * ((g() & 1) ? 1 : -1);
				c[i] = v == 0 ? mpreal(1) : v;
			}
			mpreal acc = 0;
			for (int i = 0; i + 1 < n; ++i) { x[i] = randomReal(g, D); acc += c[i] * x[i]; }
			x[n - 1] = -acc / c[n - 1];
			Run a = run(false, x), b = run(true, x);
			const bool fa = a.r.status == PslqResult::Status::Found, fb = b.r.status == PslqResult::Status::Found;
			foundS += fa; foundM += fb; tS += a.t; tM += b.t;
			for (const Run* o : { &a, &b })
				if (o->r.status == PslqResult::Status::Found)
				{
					mpreal n2 = 0; for (auto& v : o->r.relation) n2 += v * v;
					if (o->r.normBound > mpfr::sqrt(n2) * (1 + mpreal(1e-30))) printf("  trial %d: !!! %s normBound %s exceeds relation norm %s\n", tr, o == &a ? "std" : "m2l", o->r.normBound.toString(12).c_str(), mpfr::sqrt(n2).toString(12).c_str());
				}
			bool eq = fa && fb && canon(a.r.relation) == canon(b.r.relation);
			same += eq;
			if (!(fa && fb && eq))
				printf("  trial %d: std %s m2l %s (nb std %s m2l %s) planted norm %s\n", tr, toString(a.r.status), toString(b.r.status),
					a.r.normBound.toString(3).c_str(), b.r.normBound.toString(3).c_str(),
					mpfr::sqrt([&] { mpreal s = 0; for (auto& v : c) s += v * v; return s; }()).toString(3).c_str());
		}
		printf("planted D=%d n=%d coef~1e%.1f: found std %d/%d, m2l %d/%d, identical relation %d, time std %.2fs m2l %.2fs\n", D, n, lc, foundS, trials, foundM, trials, same, tS, tM);
	}
	else if (mode == "bench")
	{
		const int D = atoi(argv[2]); const int n = atoi(argv[3]);
		mpreal::set_default_prec(mpfr::digits2bits(D));
		mpreal t = mpfr::zeta(mpreal(3)) / mpfr::pow(mpfr::const_pi(), mpreal(3));
		line(("bench zeta3/pi^3 D=" + std::to_string(D) + " n=" + std::to_string(n)).c_str(), powers(t, n - 1));
	}
	printf("m2l invariant violations: %ld\n", pslqMultipair2LInvariantViolations());
}
