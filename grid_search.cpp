// Long-running, resumable search for x = alpha * e^r * pi^s
// (alpha algebraic of degree <= degree, r and s rational) with x = zeta(3) / pi^3.
//
// For every (r, s) with denominators <= maxden and |r|, |s| <= maxabs, PSLQ is run on
// (1, y, ..., y^degree) with y = x * e^(-r) * pi^(-s). Each finished pair is appended to
// <out>/results.tsv immediately, so the program can be killed and restarted at any time;
// pairs already in results.tsv are skipped. Pairs are processed in order of increasing
// max(denominator), then increasing max(|r|, |s|), so simple exponents are done first.
//
// Pruning: with D digits and degree d, PSLQ can only detect polynomials with
// |coeffs|_2 < B = 10^(0.8 D / (d + 1)). For the minimal polynomial P of alpha,
// |P|_2 >= M(P) >= max(|alpha|, 1/|alpha|) (Mahler measure), so pairs with
// |log10 alpha| > log10 B can never yield a detectable relation and are skipped.
//
// Usage:
//   grid_search [--digits 100] [--degree 10] [--maxden 10] [--maxabs 100] [--threads N] [--out results]
//   grid_search --summary [--out results]       summarize results.tsv
//   grid_search --sanity                         small self-test with a known alpha * e^r * pi^s

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <unordered_set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "pslq.h"

namespace
{
	struct Rational
	{
		int num, den;
		std::string str() const { return den == 1 ? std::to_string(num) : std::to_string(num) + "/" + std::to_string(den); }
		mpreal value() const { return mpreal(num) / den; }
	};

	struct Pair
	{
		Rational r, s;
		std::string key() const { return r.str() + "\t" + s.str(); }
		uint64_t code() const { return pack(r.num, r.den, s.num, s.den); }
		static uint64_t pack(int rn, int rd, int sn, int sd)
		{
			auto u = [](int v) { return static_cast<uint64_t>(static_cast<uint16_t>(static_cast<int16_t>(v))); };
			return (u(rn) << 48) | (u(rd) << 32) | (u(sn) << 16) | u(sd);
		}
	};

	bool parseRational(const std::string& t, int& num, int& den)
	{
		const auto slash = t.find('/');
		try
		{
			num = std::stoi(t.substr(0, slash));
			den = slash == std::string::npos ? 1 : std::stoi(t.substr(slash + 1));
		}
		catch (...)
		{
			return false;
		}
		return true;
	}

	struct Options
	{
		int digits = 100;
		int degree = 10;
		int maxDen = 10;
		int maxAbs = 100;
		unsigned threads = std::max(1u, std::thread::hardware_concurrency());
		std::string out = "results";
		bool summary = false;
		bool sanity = false;
	};

	std::vector<Rational> enumerateRationals(int maxDen, int maxAbs)
	{
		std::vector<Rational> out;
		for (int q = 1; q <= maxDen; ++q)
			for (int p = -maxAbs * q; p <= maxAbs * q; ++p)
				if (std::gcd(p < 0 ? -p : p, q) == 1)
					out.push_back({ p, q });
		return out;
	}

	// Pairs with |ln alpha| = |ln x - r - s ln pi| <= maxLnAlpha.
	std::vector<Pair> enumeratePairs(int maxDen, int maxAbs, double lnX, double maxLnAlpha, size_t& totalUnpruned)
	{
		const auto rats = enumerateRationals(maxDen, maxAbs);
		const double lnPi = std::log(3.14159265358979323846);
		totalUnpruned = rats.size() * rats.size();
		std::vector<Pair> pairs;
		for (const auto& r : rats)
			for (const auto& s : rats)
				if (std::fabs(lnX - static_cast<double>(r.num) / r.den - lnPi * s.num / s.den) <= maxLnAlpha)
					pairs.push_back({ r, s });
		auto absCeil = [](const Rational& q) { return ((q.num < 0 ? -q.num : q.num) + q.den - 1) / q.den; };
		std::stable_sort(pairs.begin(), pairs.end(), [&](const Pair& a, const Pair& b)
			{
				const auto ka = std::make_pair(std::max(a.r.den, a.s.den), std::max(absCeil(a.r), absCeil(a.s)));
				const auto kb = std::make_pair(std::max(b.r.den, b.s.den), std::max(absCeil(b.r), absCeil(b.s)));
				return ka < kb;
			});
		return pairs;
	}

	std::string polynomialToString(const std::vector<mpreal>& c)
	{
		std::string out;
		for (size_t k = c.size(); k-- > 0;)
		{
			if (c[k] == 0)
				continue;
			const bool neg = c[k] < 0;
			const std::string a = mpfr::abs(c[k]).toString("%.0Rf");
			if (!out.empty())
				out += neg ? " - " : " + ";
			else if (neg)
				out += "-";
			if (k == 0 || a != "1")
				out += a;
			if (k >= 1)
				out += "x";
			if (k >= 2)
				out += "^" + std::to_string(k);
		}
		return out + " = 0";
	}

	// x at the current default precision.
	using Target = mpreal (*)();

	mpreal zeta3OverPiCubed()
	{
		const mpreal pi = mpfr::const_pi();
		return mpfr::zeta(mpreal(3)) / (pi * pi * pi);
	}

	mpreal sanityTarget() // (2^(1/3) + 1) * e^(1/2) * pi^(-2/3)
	{
		return (mpfr::cbrt(mpreal(2)) + 1) * mpfr::exp(mpreal(1) / 2) * mpfr::pow(mpfr::const_pi(), mpreal(-2) / 3);
	}

	mpreal yAt(Target target, const Pair& p)
	{
		return target() * mpfr::exp(-p.r.value()) * mpfr::pow(mpfr::const_pi(), -p.s.value());
	}

	// Returns a description line if the hit survives verification at doubled precision.
	bool verifyHit(Target target, const Pair& p, const Options& o, std::vector<mpreal>& relation, std::string& msg)
	{
		// reduce to the minimal-degree relation
		{
			const mpreal y = yAt(target, p);
			for (int d = 1; d < o.degree; ++d)
			{
				const PslqResult rd = findIntegerPolynomial(y, d);
				if (rd.status == PslqResult::Status::Found)
				{
					relation = rd.relation;
					break;
				}
			}
		}
		const mp_prec_t saved = mpreal::get_default_prec();
		mpreal::set_default_prec(mpfr::digits2bits(2 * o.digits));
		const mpreal y2 = yAt(target, p);
		mpreal v = 0;
		for (size_t i = relation.size(); i-- > 0;)
			v = v * y2 + mpreal(relation[i]);
		v = mpfr::abs(v);
		const bool genuine = v < mpfr::pow(mpreal(10), -1.5 * o.digits);
		msg = std::string(genuine ? "HIT" : "SPURIOUS") + ": r = " + p.r.str() + ", s = " + p.s.str()
			+ ", alpha = x*e^(-r)*pi^(-s) is a root of " + polynomialToString(relation)
			+ "   (residual at " + std::to_string(2 * o.digits) + " digits: " + v.toString("%.3Re") + ")";
		mpreal::set_default_prec(saved);
		return genuine;
	}

	int runSearch(const Options& o, Target target, const std::string& name)
	{
		namespace fs = std::filesystem;
		fs::create_directories(o.out);
		const std::string resultsPath = o.out + "/results.tsv";
		const std::string hitsPath = o.out + "/hits.txt";

		// resume: collect finished pairs (only those computed with the same digits/degree)
		std::unordered_set<uint64_t> done;
		{
			std::ifstream in(resultsPath);
			std::string line;
			while (std::getline(in, line))
			{
				std::istringstream ss(line);
				std::string r, s, digits, degree;
				int rn, rd, sn, sd;
				if (std::getline(ss, r, '\t') && std::getline(ss, s, '\t') && std::getline(ss, digits, '\t') && std::getline(ss, degree, '\t')
					&& digits == std::to_string(o.digits) && degree == std::to_string(o.degree)
					&& parseRational(r, rn, rd) && parseRational(s, sn, sd))
					done.insert(Pair::pack(rn, rd, sn, sd));
			}
		}

		const double log10Bound = 0.8 * o.digits / (o.degree + 1);
		double lnX;
		{
			mpreal::set_default_prec(mpfr::digits2bits(o.digits));
			lnX = mpfr::log(target()).toDouble();
		}
		size_t totalUnpruned = 0;
		std::vector<Pair> todo = enumeratePairs(o.maxDen, o.maxAbs, lnX, log10Bound * std::log(10.0), totalUnpruned);
		const size_t nPairs = todo.size();
		todo.erase(std::remove_if(todo.begin(), todo.end(), [&](const Pair& p) { return done.count(p.code()) > 0; }), todo.end());

		std::cout << "target: " << name << "\n"
			<< "digits " << o.digits << ", degree <= " << o.degree << ", denominators <= " << o.maxDen
			<< ", |r|,|s| <= " << o.maxAbs << ", threads " << o.threads << "\n"
			<< "detectable |coeffs|_2 < 10^" << log10Bound << ", so only pairs with |log10 alpha| <= " << log10Bound
			<< " are searched: " << nPairs << " of " << totalUnpruned << " pairs\n"
			<< (nPairs - todo.size()) << " already done, " << todo.size() << " to go" << std::endl;

		std::ofstream results(resultsPath, std::ios::app);
		std::ofstream hits(hitsPath, std::ios::app);
		std::mutex ioMutex;
		std::atomic<size_t> next{ 0 }, finished{ 0 };
		const auto t0 = std::chrono::steady_clock::now();
		const size_t reportEvery = std::max<size_t>(1, std::min<size_t>(1000, todo.size() / 200));

		std::vector<std::thread> threads;
		for (unsigned t = 0; t < o.threads; ++t)
		{
			threads.emplace_back([&]()
				{
					mpreal::set_default_prec(mpfr::digits2bits(o.digits)); // thread-local in MPFR
					for (size_t k; (k = next++) < todo.size();)
					{
						const Pair& p = todo[k];
						const PslqResult r = findIntegerPolynomial(yAt(target, p), o.degree);
						std::string hitMsg;
						bool hit = false;
						std::vector<mpreal> relation = r.relation;
						if (r.status == PslqResult::Status::Found)
							hit = verifyHit(target, p, o, relation, hitMsg);

						std::unique_lock lock(ioMutex);
						results << p.key() << '\t' << o.digits << '\t' << o.degree << '\t' << toString(r.status) << '\t'
							<< (r.status == PslqResult::Status::Found ? std::string("-") : r.normBound.toString("%.4Re"))
							<< '\t' << (r.status == PslqResult::Status::Found ? polynomialToString(relation) : std::string("-"))
							<< '\n';
						results.flush();
						if (!hitMsg.empty())
						{
							hits << hitMsg << std::endl;
							std::cout << hitMsg << std::endl;
						}
						(void)hit;
						const size_t f = ++finished;
						if (f % reportEvery == 0 || f == todo.size())
						{
							const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
							const double eta = el / f * (todo.size() - f);
							std::cout << "  " << f << " / " << todo.size() << "  (max denominator so far "
								<< std::max(p.r.den, p.s.den) << ")  elapsed " << static_cast<long>(el / 60) << " min, ETA "
								<< static_cast<long>(eta / 3600) << " h " << static_cast<long>(eta / 60) % 60 << " min" << std::endl;
						}
					}
				});
		}
		for (auto& th : threads)
			th.join();
		std::cout << "finished." << std::endl;
		return 0;
	}

	int runSummary(const Options& o)
	{
		std::ifstream in(o.out + "/results.tsv");
		if (!in)
		{
			std::cerr << "no " << o.out << "/results.tsv" << std::endl;
			return 1;
		}
		struct Agg { size_t n = 0, found = 0; double minLogBound = 1e300; std::string where; };
		std::map<std::pair<std::string, std::string>, Agg> byConfig; // (digits, degree)
		std::vector<std::string> foundLines;
		std::string line;
		while (std::getline(in, line))
		{
			std::vector<std::string> f;
			std::istringstream ss(line);
			for (std::string c; std::getline(ss, c, '\t');)
				f.push_back(c);
			if (f.size() < 7)
				continue;
			Agg& a = byConfig[{ f[2], f[3] }];
			++a.n;
			if (f[4] == "Found")
			{
				++a.found;
				foundLines.push_back(line);
				continue;
			}
			const double lb = std::log10(std::stod(f[5]));
			if (lb < a.minLogBound)
			{
				a.minLogBound = lb;
				a.where = "r = " + f[0] + ", s = " + f[1];
			}
		}
		for (const auto& [cfg, a] : byConfig)
		{
			std::cout << "digits " << cfg.first << ", degree <= " << cfg.second << ": " << a.n << " pairs done, "
				<< a.found << " with a relation\n";
			if (a.n > a.found)
				std::cout << "  for all others: no integer polynomial with |coeffs|_2 < 10^" << a.minLogBound
					<< " (weakest at " << a.where << ")\n";
		}
		for (const auto& l : foundLines)
			std::cout << "  relation: " << l << "\n";
		std::cout << "(verified hits are listed in " << o.out << "/hits.txt)" << std::endl;
		return 0;
	}
}

int main(int argc, char** argv)
{
	Options o;
	for (int i = 1; i < argc; ++i)
	{
		auto next = [&]() -> const char* {
			if (i + 1 >= argc)
			{
				std::cerr << "missing value for " << argv[i] << std::endl;
				std::exit(2);
			}
			return argv[++i];
		};
		if (!std::strcmp(argv[i], "--digits")) o.digits = std::atoi(next());
		else if (!std::strcmp(argv[i], "--degree")) o.degree = std::atoi(next());
		else if (!std::strcmp(argv[i], "--maxden")) o.maxDen = std::atoi(next());
		else if (!std::strcmp(argv[i], "--maxabs")) o.maxAbs = std::atoi(next());
		else if (!std::strcmp(argv[i], "--threads")) o.threads = std::max(1, std::atoi(next()));
		else if (!std::strcmp(argv[i], "--out")) o.out = next();
		else if (!std::strcmp(argv[i], "--summary")) o.summary = true;
		else if (!std::strcmp(argv[i], "--sanity")) o.sanity = true;
		else
		{
			std::cerr << "unknown option " << argv[i] << std::endl;
			return 2;
		}
	}

	if (o.summary)
		return runSummary(o);

	if (o.sanity)
	{
		Options s = o;
		s.maxDen = 3;
		s.maxAbs = 1;
		s.degree = 8;
		s.out = o.out + "_sanity";
		std::filesystem::remove_all(s.out);
		return runSearch(s, sanityTarget, "(2^(1/3) + 1) * e^(1/2) * pi^(-2/3)");
	}

	return runSearch(o, zeta3OverPiCubed, "zeta(3) / pi^3");
}
