// Experiment018: one plaintext character per cipher symbol, including SPACE.
// Homophones are allowed. There is no boundary slot and no distinguished label.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

struct Gram {
    int order = 0;
    double weight = 0;
    array<int, 4> codes{};
};

static void require(bool condition, const string& message) {
    if (!condition) throw runtime_error(message);
}

static long long integer(const char* text, const string& name) {
    size_t end = 0;
    const string value(text);
    const long long result = stoll(value, &end);
    require(end == value.size(), "invalid " + name);
    return result;
}

static double xlogx(int64_t count) {
    require(count >= 0, "negative plaintext count");
    return count ? double(count) * log(double(count)) : 0.0;
}

static bool blank(const string& line) {
    return line.find_first_not_of(" \t\r\n") == string::npos;
}

int main(int argc, char** argv) {
    try {
        if (argc != 12) {
            cerr << "MODEL GRAMS COUNTS INITIAL SEED RESTARTS STEPS CYCLES POLISH OUTPUT A\n"
                 << "MODEL: little-endian float64 P1/P2/P3/P4, sizes A/A^2/A^3/A^4\n"
                 << "GRAMS: n weight c0 ... c(n-1); COUNTS: G then G positive counts\n"
                 << "INITIAL: one row of G labels [0,A) per restart; no boundary slot\n";
            return 2;
        }
        const long long alphabet_arg = integer(argv[11], "alphabet size");
        const long long restarts_arg = integer(argv[6], "restarts");
        const long long steps_arg = integer(argv[7], "steps");
        const long long cycles_arg = integer(argv[8], "cycles");
        const long long polish_arg = integer(argv[9], "polish");
        require(alphabet_arg >= 2 && alphabet_arg <= 64, "alphabet size must be 2..64");
        require(restarts_arg > 0 && restarts_arg <= numeric_limits<int>::max(), "invalid restart count");
        require(steps_arg >= 0 && steps_arg <= numeric_limits<int>::max(), "invalid step count");
        require(cycles_arg >= 0 && cycles_arg <= numeric_limits<int>::max(), "invalid cycle count");
        require(polish_arg >= 0 && polish_arg <= 100, "polish must be 0..100");
        require(argv[5][0] != '-', "seed must be unsigned");
        size_t seed_end = 0;
        const uint64_t seed = stoull(argv[5], &seed_end);
        require(seed_end == strlen(argv[5]), "invalid seed");
        const int A = int(alphabet_arg), restarts = int(restarts_arg);
        const int steps = int(steps_arg), cycles = int(cycles_arg), polish = int(polish_arg);

        ifstream counts_file(argv[3]);
        require(bool(counts_file), "cannot open counts");
        long long glyph_count = 0;
        require(bool(counts_file >> glyph_count) && glyph_count > 0 && glyph_count <= numeric_limits<int>::max(), "invalid G");
        const int G = int(glyph_count);
        vector<int64_t> frequencies(G);
        int64_t N = 0;
        for (int c = 0; c < G; ++c) {
            require(bool(counts_file >> frequencies[c]) && frequencies[c] > 0, "counts must contain G positive integers");
            require(frequencies[c] <= numeric_limits<int64_t>::max() - N, "count total overflow");
            N += frequencies[c];
        }
        string extra;
        require(!(counts_file >> extra), "extra count values");

        array<size_t, 5> offsets{};
        size_t entries = 0, power = 1;
        for (int n = 1; n <= 4; ++n) {
            power *= size_t(A);
            offsets[n] = entries;
            entries += power;
        }
        ifstream model_file(argv[1], ios::binary);
        require(bool(model_file), "cannot open model");
        model_file.seekg(0, ios::end);
        require(model_file.tellg() == streamoff(entries * 8), "model byte count does not match P1/P2/P3/P4 and A");
        model_file.seekg(0);
        vector<double> model(entries);
        array<unsigned char, 8> bytes{};
        static_assert(sizeof(double) == 8 && numeric_limits<double>::is_iec559, "requires IEEE754 float64");
        for (size_t i = 0; i < entries; ++i) {
            require(bool(model_file.read(reinterpret_cast<char*>(bytes.data()), 8)), "truncated model");
            uint64_t bits = 0;
            for (int j = 0; j < 8; ++j) bits |= uint64_t(bytes[j]) << (8 * j);
            memcpy(&model[i], &bits, 8);
            require(isfinite(model[i]), "model log probabilities must be finite");
        }

        ifstream grams_file(argv[2]);
        require(bool(grams_file), "cannot open grams");
        vector<Gram> grams;
        vector<double> predicted_counts(G, 0);
        double mass = 0;
        string line;
        while (getline(grams_file, line)) {
            if (blank(line)) continue;
            istringstream row(line);
            Gram gram;
            require(bool(row >> gram.order >> gram.weight), "malformed gram row");
            require(gram.order >= 1 && gram.order <= 4 && isfinite(gram.weight) && gram.weight > 0, "invalid gram order or weight");
            for (int j = 0; j < gram.order; ++j) {
                require(bool(row >> gram.codes[j]) && gram.codes[j] >= 0 && gram.codes[j] < G, "gram cipher index out of range");
            }
            require(!(row >> extra), "extra gram fields");
            grams.push_back(gram);
            require(grams.size() <= size_t(numeric_limits<int>::max()), "too many grams");
            mass += gram.weight;
            predicted_counts[gram.codes[gram.order - 1]] += gram.weight;
        }
        require(!grams.empty() && isfinite(mass), "empty or invalid grams");
        require(abs(mass - double(N)) <= 1e-10 * max(1.0, double(N)), "gram prediction mass must equal N");
        for (int c = 0; c < G; ++c) {
            require(abs(predicted_counts[c] - double(frequencies[c])) <= 1e-10 * max(1.0, double(frequencies[c])),
                    "gram final-position weights must equal each cipher frequency");
        }

        ifstream initial_file(argv[4]);
        require(bool(initial_file), "cannot open initial keys");
        vector<vector<int>> initials;
        while (getline(initial_file, line)) {
            if (blank(line)) continue;
            istringstream row(line);
            vector<int> key(G);
            for (int& value : key) {
                require(bool(row >> value) && value >= 0 && value < A, "initial row must contain G labels in [0,A)");
            }
            require(!(row >> extra), "extra initial key fields: boundary slots are not allowed");
            initials.push_back(move(key));
        }
        require(initials.size() == size_t(restarts), "initial keys must contain exactly one row per restart");

        double constant = 0;
        for (int64_t n : frequencies) constant += xlogx(n);
        vector<vector<int>> affected(G);
        for (int i = 0; i < int(grams.size()); ++i) {
            auto glyphs = grams[i].codes;
            sort(glyphs.begin(), glyphs.begin() + grams[i].order);
            int previous = -1;
            for (int j = 0; j < grams[i].order; ++j) {
                if (glyphs[j] != previous) affected[glyphs[j]].push_back(i);
                previous = glyphs[j];
            }
        }
        auto value = [&](int index, const vector<int>& key) {
            const Gram& gram = grams[index];
            size_t flat = 0;
            for (int j = 0; j < gram.order; ++j) flat = flat * size_t(A) + size_t(key[gram.codes[j]]);
            return model[offsets[gram.order] + flat];
        };
        auto full = [&](const vector<int>& key) {
            double score = constant;
            vector<int64_t> totals(A, 0);
            for (int c = 0; c < G; ++c) totals[key[c]] += frequencies[c];
            for (int64_t n : totals) score -= xlogx(n);
            for (int i = 0; i < int(grams.size()); ++i) score += grams[i].weight * value(i, key);
            require(isfinite(score), "nonfinite objective");
            return score;
        };
        mt19937_64 rng(seed);
        uniform_real_distribution<double> uniform(0, 1);
        uniform_int_distribution<int> glyph(0, G - 1), label(0, A - 1);
        discrete_distribution<int> weighted(frequencies.begin(), frequencies.end());
        ofstream output(argv[10]);
        require(bool(output), "cannot open output");
        output << setprecision(17);
        cerr << setprecision(17);
        vector<uint64_t> marks(grams.size(), 0);
        uint64_t tick = 0;

        for (int restart = 0; restart < restarts; ++restart) {
            vector<int> key = initials[restart], best_key = key;
            vector<int64_t> totals(A, 0), increments(A, 0);
            vector<double> cache(grams.size());
            double current = 0, best = full(key), largest_error = 0;
            int cache_checks = 0, polish_passes = 0;
            array<uint64_t, 3> proposals{}, accepted{};
            vector<uint64_t> replacement_labels(A, 0), renamed_labels(A, 0);
            auto reset = [&]() {
                fill(totals.begin(), totals.end(), 0);
                for (int c = 0; c < G; ++c) totals[key[c]] += frequencies[c];
                for (int i = 0; i < int(grams.size()); ++i) cache[i] = value(i, key);
                current = full(key);
            };
            auto inspect_cache = [&]() {
                vector<int64_t> reconstructed(A, 0);
                for (int c = 0; c < G; ++c) reconstructed[key[c]] += frequencies[c];
                require(reconstructed == totals, "cache plaintext totals differ from full recount");
                for (int i = 0; i < int(grams.size()); ++i) {
                    require(cache[i] == value(i, key), "cached gram differs from full lookup");
                }
                const double exact = full(key), exact_best = full(best_key);
                const double error = max(abs(current - exact), abs(best - exact_best));
                largest_error = max(largest_error, error);
                require(error <= 1e-9 * max(1.0, double(N)), "incremental objective differs from full score");
                current = exact;
                best = exact_best;
                if (current > best) { best = current; best_key = key; }
                ++cache_checks;
            };
            reset();
            vector<int> old, impacted;
            auto proposal = [&](const vector<int>& changed, const vector<int>& replacement, double temperature, int type) {
                require(changed.size() == replacement.size(), "proposal size mismatch");
                ++proposals[type];
                old.clear(); impacted.clear();
                fill(increments.begin(), increments.end(), 0);
                if (++tick == 0) { fill(marks.begin(), marks.end(), 0); ++tick; }
                for (size_t j = 0; j < changed.size(); ++j) {
                    const int c = changed[j];
                    old.push_back(key[c]);
                    increments[key[c]] -= frequencies[c];
                    increments[replacement[j]] += frequencies[c];
                    key[c] = replacement[j];
                    for (int i : affected[c]) {
                        if (marks[i] != tick) { marks[i] = tick; impacted.push_back(i); }
                    }
                }
                double delta = 0;
                for (int t = 0; t < A; ++t) {
                    if (increments[t]) delta += xlogx(totals[t]) - xlogx(totals[t] + increments[t]);
                }
                for (int i : impacted) delta += grams[i].weight * (value(i, key) - cache[i]);
                require(isfinite(delta), "nonfinite proposal delta");
                const bool take = temperature > 0 ? delta >= 0 || uniform(rng) < exp(delta / temperature) : delta > 1e-8;
                if (take) {
                    ++accepted[type];
                    current += delta;
                    for (int t = 0; t < A; ++t) totals[t] += increments[t];
                    for (int i : impacted) cache[i] = value(i, key);
                    if (current > best) { best = current; best_key = key; }
                } else {
                    for (size_t j = 0; j < changed.size(); ++j) key[changed[j]] = old[j];
                }
                return take;
            };
            for (int cycle = 0; cycle < cycles; ++cycle) {
                key = best_key;
                reset();
                for (int step = 0; step < steps; ++step) {
                    const double progress = steps > 1 ? double(step) / double(steps - 1) : 0.0;
                    const double temperature = double(N) * .02 * pow(.00005 / .02, progress);
                    const double kind = uniform(rng);
                    if (kind < .8) {
                        const int c = uniform(rng) < .5 ? glyph(rng) : weighted(rng), next = label(rng);
                        ++replacement_labels[next];
                        if (key[c] != next) proposal({c}, {next}, temperature, 0);
                    } else if (kind < .95) {
                        const int a = glyph(rng), b = glyph(rng);
                        if (a != b && key[a] != key[b]) proposal({a, b}, {key[b], key[a]}, temperature, 1);
                    } else {
                        const int a = label(rng), b = label(rng);
                        if (a == b) continue;
                        ++renamed_labels[a]; ++renamed_labels[b];
                        vector<int> changed, replacements;
                        for (int c = 0; c < G; ++c) {
                            if (key[c] == a || key[c] == b) {
                                changed.push_back(c);
                                replacements.push_back(key[c] == a ? b : a);
                            }
                        }
                        if (!changed.empty()) proposal(changed, replacements, temperature, 2);
                    }
                }
                inspect_cache();
            }
            key = best_key;
            reset();
            for (int pass = 0; pass < polish; ++pass) {
                bool improved = false;
                ++polish_passes;
                for (int c = 0; c < G; ++c) {
                    for (int next = 0; next < A; ++next) {
                        if (next != key[c]) improved = proposal({c}, {next}, 0, 0) || improved;
                    }
                }
                if (!improved) break;
            }
            inspect_cache();
            const double score = full(best_key) / double(N);
            output << restart << ' ' << score;
            for (int next : best_key) output << ' ' << next;
            output << '\n'; output.flush();
            require(bool(output), "cannot write output");
            auto emit_array = [&](const auto& values) {
                cerr << '[';
                for (size_t i = 0; i < values.size(); ++i) { if (i) cerr << ','; cerr << values[i]; }
                cerr << ']';
            };
            cerr << "{\"restart\":" << restart << ",\"score\":" << score
                 << ",\"N\":" << N << ",\"prediction_weight\":" << mass
                 << ",\"cache_checks\":" << cache_checks << ",\"max_cache_error\":" << largest_error
                 << ",\"polish_passes\":" << polish_passes << ",\"proposals_replace_swap_rename\":";
            emit_array(proposals);
            cerr << ",\"accepted_replace_swap_rename\":"; emit_array(accepted);
            cerr << ",\"replacement_labels_considered\":"; emit_array(replacement_labels);
            cerr << ",\"rename_labels_considered\":"; emit_array(renamed_labels);
            cerr << "}\n";
        }
        return 0;
    } catch (const exception& error) {
        cerr << "gap_search: " << error.what() << '\n';
        return 1;
    }
}
