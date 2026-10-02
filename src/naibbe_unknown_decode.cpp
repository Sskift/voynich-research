#include "naibbe_unknown.hpp"
#include <map>
#include <functional>

namespace nu {
namespace {
constexpr int BOS_EOS = A;
constexpr int START_HISTORY = BOS_EOS * 24 + BOS_EOS;
using Outgoing = std::map<int, std::vector<int>>;

Outgoing outgoing(const Row& row, const Inventory* inventory = nullptr) {
    need(row.positions > 1 && row.tokens > 0,
         "DECODER_INVALID_INPUT: invalid row dimensions");
    need(row.edges.size() <= static_cast<size_t>(std::numeric_limits<int>::max()),
         "DECODER_INVALID_INPUT: too many edges");
    Outgoing result;
    for (size_t id = 0; id < row.edges.size(); ++id) {
        const Edge& e = row.edges[id];
        need(e.from >= 0 && e.to > e.from && e.to < row.positions &&
             e.first >= 0 && e.second >= -1 && e.token >= 0 &&
             e.token < row.tokens && std::isfinite(e.constant),
             "DECODER_INVALID_INPUT: invalid edge");
        if (inventory) {
            const Inventory& in = *inventory;
            need(e.first < in.G &&
                 ((e.second == -1 && in.role[e.first] == 0) ||
                  (e.second >= 0 && e.second < in.G &&
                   in.role[e.first] == 1 && in.role[e.second] == 2)),
                 "DECODER_INVALID_INPUT: invalid edge roles");
        }
        result[e.from].push_back(static_cast<int>(id));
    }
    return result;
}

void validate_scoring_inputs(const Inventory& in, const Model& model, const Key& key) {
    need(in.G > 0 && in.G <= S && in.role.size() == static_cast<size_t>(in.G) &&
         in.multiplicity.size() == static_cast<size_t>(in.G),
         "DECODER_INVALID_INPUT: invalid inventory");
    for (int g = 0; g < in.G; ++g)
        need(in.role[g] >= 0 && in.role[g] < R && in.multiplicity[g] > 0,
             "DECODER_INVALID_INPUT: invalid inventory member");
    std::set<std::pair<int, int>> seen;
    for (auto [p, s] : in.forbidden)
        need(p >= 0 && p < in.G && s >= 0 && s < in.G &&
             in.role[p] == 1 && in.role[s] == 2 && seen.emplace(p, s).second,
             "DECODER_INVALID_INPUT: invalid forbidden pair");
    Key rebuilt = key;
    need(rebuilt.rebuild(in) && rebuilt.letter == key.letter &&
         rebuilt.mass == key.mass && rebuilt.blocked == key.blocked,
         "DECODER_INVALID_INPUT: invalid or stale key cache");
    need(model.logp.size() == 24 * 24 * 24,
         "DECODER_INVALID_INPUT: invalid model size");
    for (double logp : model.logp)
        need(std::isfinite(logp) && logp <= 0,
             "DECODER_INVALID_INPUT: invalid model log probability");
}

struct State {
    double score = NEG;
    int previous_position = -1;
    int previous_history = -1;
    int previous_edge = -1;
};

double finite_sum(double first, double second) {
    double value = first + second;
    need(std::isfinite(value), "DECODER_NUMERICAL_FAILURE: nonfinite path score");
    return value;
}

DecodedRow decode_row(const Inventory& in, const Row& row, const Model& model,
                      const Key& key, long long max_states) {
    const Outgoing edges = outgoing(row, &in);
    std::vector<double> channels;
    channels.reserve(row.edges.size());
    for (const Edge& edge : row.edges) {
        double channel = edge_channel(edge, key);
        need(std::isfinite(channel),
             "DECODER_NUMERICAL_FAILURE: nonfinite channel score");
        channels.push_back(channel);
    }
    // Ordered sparse maps avoid allocating positions * 576 entries before the
    // reachable-state budget is checked. Insertion cannot invalidate iterators.
    std::map<int, std::map<int, State>> states;
    states[0][START_HISTORY].score = 0;
    DecodedRow result;
    result.states = 1;
    for (auto position = states.begin(); position != states.end(); ++position) {
        auto next_edges = edges.find(position->first);
        if (next_edges == edges.end()) continue;
        for (const auto& [history, state] : position->second) {
            for (int edge_id : next_edges->second) {
                need(result.transitions < std::numeric_limits<long long>::max(),
                     "DECODER_BUDGET_EXCEEDED: transition counter overflow");
                ++result.transitions;
                const Edge& edge = row.edges[edge_id];
                double score = finite_sum(state.score, channels[edge_id]);
                int new_history = history;
                for (int glyph : {edge.first, edge.second}) {
                    if (glyph < 0) continue;
                    int letter = key.letter[glyph];
                    score = finite_sum(score, model.at(new_history, letter));
                    new_history = (new_history % 24) * 24 + letter;
                }
                auto& destination = states[edge.to];
                auto found = destination.find(new_history);
                if (found == destination.end()) {
                    need(result.states < max_states,
                         "DECODER_BUDGET_EXCEEDED: reachable states exceed max_states");
                    ++result.states;
                    destination.emplace(new_history,
                                        State{score, position->first, history, edge_id});
                } else if (score > found->second.score) {
                    found->second = State{score, position->first, history, edge_id};
                }
            }
        }
    }
    auto terminal = states.find(row.positions - 1);
    need(terminal != states.end() && !terminal->second.empty(),
         "DECODER_NO_COMPLETE_PATH: no path to the row endpoint");
    int best_history = -1;
    for (const auto& [history, state] : terminal->second) {
        double score = finite_sum(state.score, model.at(history, BOS_EOS));
        if (score > result.score) {
            result.score = score;
            best_history = history;
        }
    }
    int position = row.positions - 1;
    int history = best_history;
    while (position != 0) {
        const State& state = states.at(position).at(history);
        need(state.previous_edge >= 0 && state.previous_position < position,
             "DECODER_INTERNAL_ERROR: invalid traceback");
        result.edges.push_back(state.previous_edge);
        position = state.previous_position;
        history = state.previous_history;
    }
    need(history == START_HISTORY, "DECODER_INTERNAL_ERROR: wrong traceback origin");
    std::reverse(result.edges.begin(), result.edges.end());
    for (int id : result.edges) {
        const Edge& edge = row.edges[id];
        result.atoms.push_back(edge.first);
        result.letters.push_back(key.letter[edge.first]);
        if (edge.second >= 0) {
            result.atoms.push_back(edge.second);
            result.letters.push_back(key.letter[edge.second]);
        }
    }
    return result;
}

// Synthetic helpers are intentionally independent of production edge scoring
// and Viterbi recurrences: enumerate all complete paths, then score each path.
struct Toy {
    Inventory inventory;
    Key key;
    Model model;
    int glyph(int role, int table, int letter) const {
        return key.slots[(role * T + table) * A + letter];
    }
};

Toy toy() {
    Toy out;
    for (int r = 0; r < R; ++r) {
        for (int t = 0; t < T; ++t) {
            for (int a = 0; a < A; ++a) {
                int g;
                if (r == 0 && t == 3 && a == 0) {
                    g = out.key.slots[(r * T + 2) * A];
                    ++out.inventory.multiplicity[g];
                    out.inventory.duplicate = g;
                } else {
                    g = out.inventory.G++;
                    out.inventory.role.push_back(r);
                    out.inventory.multiplicity.push_back(1);
                }
                out.key.slots[(r * T + t) * A + a] = g;
            }
        }
    }
    out.inventory.forbidden.emplace_back(out.glyph(1, 0, 0), out.glyph(2, 0, 0));
    need(out.key.rebuild(out.inventory), "Decoder self-check: invalid synthetic key");
    out.model.logp.resize(24 * 24 * 24);
    for (int h = 0; h < H; ++h) {
        double denominator = 0;
        for (int c = 0; c < 24; ++c) denominator += 1 + ((h * 13 + c * 7) % 29);
        for (int c = 0; c < 24; ++c)
            out.model.logp[h * 24 + c] =
                std::log((1 + ((h * 13 + c * 7) % 29)) / denominator);
    }
    return out;
}

struct Brute {
    double score = NEG;
    std::vector<int> shortest;
    std::vector<std::vector<int>> paths;
    std::set<std::pair<int, int>> states;
    long long transitions = 0;
};

double independent_score(const Row& row, const std::vector<int>& path,
                         const Key& key, const Model& model) {
    double score = 0;
    std::vector<int> text{BOS_EOS, BOS_EOS};
    for (int id : path) {
        const Edge& edge = row.edges[id];
        double channel = edge.constant + std::log(static_cast<double>(key.mass[edge.first]) / 52.0);
        if (edge.second >= 0) {
            int index = key.letter[edge.first] * A + key.letter[edge.second];
            channel += std::log(static_cast<double>(key.mass[edge.second]) / 52.0)
                     - std::log(static_cast<double>(2704 - key.blocked[index]) / 2704.0);
        }
        score += channel;
        for (int glyph : {edge.first, edge.second}) {
            if (glyph < 0) continue;
            int h = text[text.size() - 2] * 24 + text.back();
            text.push_back(key.letter[glyph]);
            score += model.logp[h * 24 + text.back()];
        }
    }
    score += model.logp[(text[text.size() - 2] * 24 + text.back()) * 24 + BOS_EOS];
    return score;
}

Brute enumerate(const Row& row, const Key& key, const Model& model) {
    Brute result;
    std::vector<int> path;
    std::function<void(int, int)> visit = [&](int position, int history) {
        result.states.emplace(position, history);
        if (position == row.positions - 1) {
            result.paths.push_back(path);
            double score = independent_score(row, path, key, model);
            if (score > result.score) result.score = score;
            if (result.shortest.empty() || path.size() < result.shortest.size() ||
                (path.size() == result.shortest.size() && path < result.shortest))
                result.shortest = path;
            return;
        }
        for (size_t id = 0; id < row.edges.size(); ++id) {
            const Edge& edge = row.edges[id];
            if (edge.from != position) continue;
            int new_history = history;
            for (int glyph : {edge.first, edge.second})
                if (glyph >= 0) new_history = (new_history % 24) * 24 + key.letter[glyph];
            path.push_back(static_cast<int>(id));
            visit(edge.to, new_history);
            path.pop_back();
        }
    };
    visit(0, START_HISTORY);
    for (const auto& state : result.states)
        for (const Edge& edge : row.edges)
            if (edge.from == state.first) ++result.transitions;
    return result;
}

void close_score(double observed, double expected) {
    need(std::isfinite(observed) && std::isfinite(expected) &&
         std::abs(observed - expected) <= 1e-11 * std::max(1.0, std::abs(expected)),
         "Decoder self-check: exhaustive score mismatch");
}

void compare_exhaustive(const Toy& data, const Row& row) {
    Brute reference = enumerate(row, data.key, data.model);
    need(!reference.paths.empty(), "Decoder self-check: toy unexpectedly has no path");
    Decoded decoded = decode(data.inventory, Lattice{{row}}, data.model, data.key);
    const DecodedRow& result = decoded.rows.at(0);
    close_score(result.score, reference.score);
    close_score(decoded.score, reference.score);
    close_score(independent_score(row, result.edges, data.key, data.model), reference.score);
    need(std::find(reference.paths.begin(), reference.paths.end(), result.edges) != reference.paths.end(),
         "Decoder self-check: traceback is not a complete path");
    need(result.states == static_cast<long long>(reference.states.size()) &&
         result.transitions == reference.transitions && decoded.max_states == result.states,
         "Decoder self-check: state or transition count mismatch");
    std::vector<int> atoms, letters;
    for (int id : result.edges)
        for (int g : {row.edges[id].first, row.edges[id].second})
            if (g >= 0) { atoms.push_back(g); letters.push_back(data.key.letter[g]); }
    need(result.atoms == atoms && result.letters == letters,
         "Decoder self-check: atom or letter traceback mismatch");
    need(initial_parse(Lattice{{row}}).at(0) == reference.shortest,
         "Decoder self-check: shortest-path canonical tie mismatch");
    Decoded exact_budget = decode(data.inventory, Lattice{{row}}, data.model, data.key, result.states);
    need(exact_budget.rows[0].edges == result.edges,
         "Decoder self-check: exact state budget changed decode");
}

template <class F> void expect_error(F function, const std::string& marker) {
    try { function(); }
    catch (const std::runtime_error& error) {
        need(std::string(error.what()).find(marker) != std::string::npos,
             "Decoder self-check: incorrect failure category");
        return;
    }
    throw std::runtime_error("Decoder self-check: required failure did not occur");
}
} // namespace

Decoded decode(const Inventory& in, const Lattice& lattice, const Model& model,
               const Key& key, long long max_states) {
    need(max_states > 0, "DECODER_INVALID_INPUT: max_states must be positive");
    need(!lattice.rows.empty(), "DECODER_INVALID_INPUT: empty lattice");
    validate_scoring_inputs(in, model, key);
    Decoded result;
    for (const Row& row : lattice.rows) {
        DecodedRow decoded = decode_row(in, row, model, key, max_states);
        result.score = finite_sum(result.score, decoded.score);
        result.max_states = std::max(result.max_states, decoded.states);
        result.rows.push_back(std::move(decoded));
    }
    return result;
}

std::vector<std::vector<int>> initial_parse(const Lattice& lattice) {
    need(!lattice.rows.empty(), "DECODER_INVALID_INPUT: empty lattice");
    std::vector<std::vector<int>> result;
    for (const Row& row : lattice.rows) {
        const Outgoing edges = outgoing(row);
        std::map<int, int> distance{{row.positions - 1, 0}};
        for (auto source = edges.rbegin(); source != edges.rend(); ++source) {
            int shortest = std::numeric_limits<int>::max();
            for (int id : source->second) {
                auto suffix = distance.find(row.edges[id].to);
                if (suffix != distance.end()) shortest = std::min(shortest, suffix->second + 1);
            }
            if (shortest != std::numeric_limits<int>::max()) distance[source->first] = shortest;
        }
        need(distance.count(0) != 0,
             "DECODER_NO_COMPLETE_PATH: initial_parse cannot reach the row endpoint");
        std::vector<int> path;
        int position = 0;
        while (position != row.positions - 1) {
            bool found = false;
            for (int id : edges.at(position)) {
                const Edge& edge = row.edges[id];
                auto suffix = distance.find(edge.to);
                if (suffix != distance.end() && suffix->second + 1 == distance.at(position)) {
                    path.push_back(id);
                    position = edge.to;
                    found = true;
                    break;
                }
            }
            need(found, "DECODER_INTERNAL_ERROR: shortest path traceback failed");
        }
        result.push_back(std::move(path));
    }
    return result;
}

void decoder_self_check() {
    Toy data = toy();
    auto u = [&](int a, int t = 0) { return data.glyph(0, t, a); };
    auto p = [&](int a, int t = 0) { return data.glyph(1, t, a); };
    auto s = [&](int a, int t = 0) { return data.glyph(2, t, a); };
    Row ambiguous{5, 2, {
        {2, 4, u(3), -1, 1, std::log(.97)},
        {0, 2, p(0), s(0), 0, std::log(19.0 / 36)},
        {0, 1, u(1), -1, 0, std::log(17.0 / 36)},
        {1, 2, u(2), -1, 0, std::log(.03)},
        {2, 3, u(1, 2), -1, 1, std::log(.03)},
        {3, 4, p(2, 1), s(3, 4), 1, std::log(19.0 / 36)},
        {0, 2, u(0, 2), -1, 0, std::log(17.0 / 36)},
        {2, 4, p(1, 1), s(2, 1), 1, std::log(19.0 / 36)}
    }};
    compare_exhaustive(data, ambiguous);
    // A varied deterministic family covers reconverging histories, duplicate
    // edges, multiple output lengths, dead ends, and shuffled edge positions.
    for (int seed = 0; seed < 48; ++seed) {
        Row row{6, 2, {}};
        for (int from = 0; from < 5; ++from) {
            int a = (seed + 3 * from) % 5;
            row.edges.push_back({from, from + 1, u(a, seed % T), -1,
                                 from < 3 ? 0 : 1, -.03 * (from + 1)});
            if (from + 2 <= 5)
                row.edges.push_back({from, from + 2, p((a + 1) % 5, (seed + 1) % T),
                                     s((a + 2) % 5, (seed + 2) % T),
                                     from < 3 ? 0 : 1, -.1 * (seed % 3)});
            if ((seed + from) % 3 == 0)
                row.edges.push_back({from, from + 1, u((a + 2) % 5, (seed + 3) % T),
                                     -1, from < 3 ? 0 : 1, -.2});
        }
        std::rotate(row.edges.begin(), row.edges.begin() + seed % row.edges.size(), row.edges.end());
        compare_exhaustive(data, row);
    }
    // Uniform scores give an exact same-state edge tie: first input edge wins.
    Toy tied = data;
    std::fill(tied.model.logp.begin(), tied.model.logp.end(), -std::log(24.0));
    Row duplicate_edges{2, 1, {{0, 1, u(0), -1, 0, 0}, {0, 1, u(0), -1, 0, 0}}};
    need(decode(tied.inventory, Lattice{{duplicate_edges}}, tied.model, tied.key).rows[0].edges ==
         std::vector<int>{0}, "Decoder self-check: input edge tie is unstable");
    // Final-state ties choose ascending numeric history, even when its edge is later.
    Row history_tie{2, 1, {{0, 1, u(1), -1, 0, 0}, {0, 1, u(0), -1, 0, 0}}};
    need(decode(tied.inventory, Lattice{{history_tie}}, tied.model, tied.key).rows[0].edges ==
         std::vector<int>{1}, "Decoder self-check: terminal history tie is unstable");
    // EOS alone reverses the pre-EOS ordering.
    Toy eos = tied;
    eos.model.logp[(BOS_EOS * 24 + 1) * 24 + BOS_EOS] = -.01;
    eos.model.logp[(BOS_EOS * 24 + 0) * 24 + BOS_EOS] = -20;
    compare_exhaustive(eos, history_tie);
    need(decode(eos.inventory, Lattice{{history_tie}}, eos.model, eos.key).rows[0].letters ==
         std::vector<int>{1}, "Decoder self-check: missing EOS");
    // The second token must use the two emitted letters, never BOS again.
    Row cross_token{3, 2, {{0, 1, p(0), s(1), 0, 0},
                          {1, 2, u(2), -1, 1, 0}, {1, 2, u(3), -1, 1, 0}}};
    Toy context = tied;
    context.model.logp[(0 * 24 + 1) * 24 + 2] = -.01;
    context.model.logp[(0 * 24 + 1) * 24 + 3] = -20;
    context.model.logp[START_HISTORY * 24 + 2] = -20;
    context.model.logp[START_HISTORY * 24 + 3] = -.01;
    compare_exhaustive(context, cross_token);
    need(decode(context.inventory, Lattice{{cross_token}}, context.model, context.key).rows[0].letters ==
         std::vector<int>({0, 1, 2}), "Decoder self-check: context reset at token boundary");
    // Greedy first-edge tie breaking must consider shortest suffix distance.
    Row shortest{5, 1, {{0, 1, u(0), -1, 0, 0}, {0, 3, u(1), -1, 0, 0},
                       {1, 2, u(2), -1, 0, 0}, {2, 4, u(3), -1, 0, 0},
                       {3, 4, u(4), -1, 0, 0}, {0, 2, u(5), -1, 0, 0}}};
    need(initial_parse(Lattice{{shortest}})[0] == std::vector<int>({1, 4}),
         "Decoder self-check: shortest parse canonical order is wrong");
    compare_exhaustive(data, shortest);
    Row disconnected{3, 1, {{0, 1, u(0), -1, 0, 0}}};
    expect_error([&] { decode(data.inventory, Lattice{{disconnected}}, data.model, data.key); },
                 "DECODER_NO_COMPLETE_PATH");
    expect_error([&] { initial_parse(Lattice{{disconnected}}); }, "DECODER_NO_COMPLETE_PATH");
    Row no_edges{2, 1, {}};
    expect_error([&] { decode(data.inventory, Lattice{{no_edges}}, data.model, data.key); },
                 "DECODER_NO_COMPLETE_PATH");
    expect_error([&] { initial_parse(Lattice{{no_edges}}); }, "DECODER_NO_COMPLETE_PATH");
    long long count = decode(data.inventory, Lattice{{ambiguous}}, data.model, data.key).rows[0].states;
    expect_error([&] { decode(data.inventory, Lattice{{ambiguous}}, data.model, data.key, count - 1); },
                 "DECODER_BUDGET_EXCEEDED");
    expect_error([&] { decode(data.inventory, Lattice{{duplicate_edges}}, data.model, data.key, 1); },
                 "DECODER_BUDGET_EXCEEDED");
    expect_error([&] { decode(data.inventory, Lattice{{duplicate_edges}}, data.model, data.key, 0); },
                 "DECODER_INVALID_INPUT");
    expect_error([&] { decode(data.inventory, Lattice{}, data.model, data.key); }, "DECODER_INVALID_INPUT");
    expect_error([&] { initial_parse(Lattice{}); }, "DECODER_INVALID_INPUT");
    Decoded two_rows = decode(data.inventory, Lattice{{ambiguous, ambiguous}}, data.model, data.key, count);
    need(two_rows.max_states == count && two_rows.rows[0].edges == two_rows.rows[1].edges,
         "Decoder self-check: budget or context leaked across rows");
    close_score(two_rows.score, 2 * two_rows.rows[0].score);
    Key stale = data.key;
    ++stale.mass[0];
    expect_error([&] { decode(data.inventory, Lattice{{duplicate_edges}}, data.model, stale); },
                 "DECODER_INVALID_INPUT");
}
} // namespace nu
