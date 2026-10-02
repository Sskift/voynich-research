"""Independent exhaustive-frontier OOV top-two audit for frozen020 v2.

The independent recovery audit supplies only its loader/hash checker. No
production module or production search routine is imported. Unlike production,
this graph has no lexical-trie state. Lexicon membership is tested at EOS only.
Exact Python-integer dyadic bounds replace floating-point A* heuristics: each
represented edge cost is enclosed on the 2**-60 lattice, reverse Dijkstra uses
lower integer costs, and distinct string prefixes are never merged.

This audits search under independently normalized represented log tables. It
does not certify infinite probability mass, the HMM, or full arithmetic error.
"""
from collections import Counter, deque
from functools import lru_cache
from pathlib import Path
import argparse
import datetime
import heapq
import itertools
import json
import math
import sys
import traceback

import numpy as np

from audit_omission_recovery_v2 import IndependentModel, Integrity, COUNTS as HASH_COUNTS


ROOT = Path(__file__).resolve().parents[1]
DATA = 'data/omission-oracle'
RESULTS = 'results/omission-oracle-v2'
INPUT = 'data/omission-oracle-v2/input-freeze.json'
OUTPUT = 'results/audit-016-omission-search-v2.json'
VOWELS = frozenset('aeiou')
SCALE = 1 << 60
RANK_TOLERANCE = 1e-12
LOG_AGREEMENT_TOLERANCE = 1e-10
EXPANSION_CAP = 500000
HEAP_CAP = 3000000


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def sha(path):
    import hashlib
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def enclosure(cost):
    """Exact enclosure of a represented nonnegative float, no FP multiplication."""
    need(math.isfinite(cost) and cost >= 0, 'Invalid represented edge cost')
    numerator, denominator = float(cost).as_integer_ratio()
    lower, remainder = divmod(numerator*SCALE, denominator)
    return lower, lower+bool(remainder)


TIE_UNITS = enclosure(RANK_TOLERANCE)[1]


def image(word):
    return ''.join(letter for letter in word if letter not in VOWELS) or word


class SearchModel:
    def __init__(self, model):
        self.model = model
        self.lexicon = frozenset(model.words)
        self.lexical_images = {}
        for word in sorted(self.lexicon):
            self.lexical_images.setdefault(image(word), []).append(word)
        self.vowels = tuple(c for c in model.alphabet if c in VOWELS)
        self.index = model.index
        self.log_nonempty = model.log_nonempty
        self.edge = lru_cache(maxsize=None)(self._edge)

    def _edge(self, history, letter):
        coordinates = tuple(self.index[c] for c in history+(letter,))
        cost = -float(self.model.tables[len(history)][coordinates])
        need(cost > 0, 'Search requires strictly positive represented edge costs')
        lower, upper = enclosure(cost)
        return cost, lower, upper

    def word(self, word):
        need(word and all(c in self.index and c != ' ' for c in word), 'Invalid word')
        history, costs, low, high = (' ',), [], 0, 0
        for letter in word+' ':
            cost, lower, upper = self.edge(history, letter)
            costs.append(cost)
            low += lower
            high += upper
            history = (history+(letter,))[-3:]
        raw = math.fsum(costs)
        return dict(word=word, raw_cost=raw, lower=low, upper=high,
                    log_char_probability=-raw-self.log_nonempty)


class SearchGraph:
    def __init__(self, model, skeleton):
        need(skeleton and not any(c in VOWELS for c in skeleton), 'Expected consonant skeleton')
        self.model, self.skeleton = model, skeleton
        start = (0, (' ',))
        self.states, indices, queue = [start], {start: 0}, deque([0])
        self.edges, self.eos = [], []
        while queue:
            index = queue.popleft()
            position, history = self.states[index]
            self.eos.append(model.edge(history, ' ') if position == len(skeleton) else None)
            edges = []
            letters = model.vowels+((skeleton[position],) if position < len(skeleton) else ())
            for letter in letters:
                target = (position+int(letter not in VOWELS), (history+(letter,))[-3:])
                if target not in indices:
                    indices[target] = len(self.states)
                    self.states.append(target)
                    queue.append(indices[target])
                edges.append((letter, indices[target], model.edge(history, letter)))
            self.edges.append(edges)
        reverse = [[] for _ in self.states]
        distances = [None]*len(self.states)
        heap = []
        for node, edges in enumerate(self.edges):
            if self.eos[node] is not None:
                distances[node] = self.eos[node][1]
                heapq.heappush(heap, (distances[node], node))
            for _, target, (_, lower, _) in edges:
                reverse[target].append((node, lower))
        while heap:
            distance, node = heapq.heappop(heap)
            if distance != distances[node]:
                continue
            for previous, cost in reverse[node]:
                candidate = distance+cost
                if distances[previous] is None or candidate < distances[previous]:
                    distances[previous] = candidate
                    heapq.heappush(heap, (candidate, previous))
        need(all(value is not None for value in distances), 'Reachable state cannot finish skeleton')
        # A posteriori exact integer inequalities additionally check Dijkstra.
        for node, edges in enumerate(self.edges):
            if self.eos[node] is not None:
                need(distances[node] <= self.eos[node][1], 'Invalid terminal lower heuristic')
            need(all(distances[node] <= edge[2][1]+distances[edge[1]] for edge in edges),
                 'Invalid integer Bellman lower heuristic')
        self.distance = distances

    def enumerate(self, expansion_cap=EXPANSION_CAP, heap_cap=HEAP_CAP):
        # Heap keys and all path/heuristic bounds are integers. EOS is popped,
        # never accepted merely because a prefix generated it.
        heap = [(self.distance[0], '', False, 0, 0, 0)]
        found, best_upper = [], []
        expanded, generated, peak, rejected = 0, 1, 1, 0
        stop, complete = 'FRONTIER_EXHAUSTED', False
        while heap:
            if len(best_upper) >= 2 and heap[0][0] > best_upper[1]+TIE_UNITS:
                stop, complete = 'TOP2_EXACT_INTEGER_FRONTIER', True
                break
            priority, prefix, terminal, state, lower, upper = heapq.heappop(heap)
            if terminal:
                if prefix in self.model.lexicon:
                    rejected += 1
                else:
                    word = self.model.word(prefix)
                    need((word['lower'], word['upper']) == (lower, upper), 'Path bound inconsistency')
                    found.append(word)
                    best_upper = sorted(best_upper+[upper])[:2]
                continue
            if expanded >= expansion_cap:
                heapq.heappush(heap, (priority, prefix, terminal, state, lower, upper))
                stop = 'EXPANDED_PREFIX_CAP'
                break
            expanded += 1
            additions = []
            if self.eos[state] is not None:
                _, a, b = self.eos[state]
                additions.append((lower+a, prefix, True, state, lower+a, upper+b))
            for letter, target, (_, a, b) in self.edges[state]:
                additions.append((lower+a+self.distance[target], prefix+letter, False,
                                  target, lower+a, upper+b))
            if len(heap)+len(additions) > heap_cap:
                heapq.heappush(heap, (priority, prefix, terminal, state, lower, upper))
                stop = 'HEAP_CAP'
                break
            for entry in additions:
                heapq.heappush(heap, entry)
            generated += len(additions)
            peak = max(peak, len(heap))
        else:
            complete = True
        found.sort(key=lambda entry: (entry['raw_cost'], entry['word']))
        frontier = heap[0][0] if heap else None
        return dict(search_complete=complete, found=found, stop_reason=stop,
                    expanded_prefixes=expanded, generated_entries=generated,
                    heap_peak=peak, lexicon_eos_rejected=rejected,
                    graph_states=len(self.states), graph_edges=sum(map(len, self.edges)),
                    expansion_cap=expansion_cap, heap_cap=heap_cap,
                    frontier_lower_integer=frontier,
                    second_terminal_upper_integer=best_upper[1] if len(best_upper) >= 2 else None,
                    dyadic_scale=SCALE, tie_margin_integer=TIE_UNITS)


def rank_interval_check(candidate, lower, upper):
    # Rank tolerance is the exact represented float, not its slightly wider
    # lattice ceiling. Rational integer cross-products decide the comparison.
    numerator, denominator = RANK_TOLERANCE.as_integer_ratio()
    worst = max(candidate['upper']-lower, upper-candidate['lower'], 0)
    best = max(candidate['lower']-upper, lower-candidate['upper'], 0)
    if worst*denominator <= numerator*SCALE:
        status = 'PASS'
    elif best*denominator > numerator*SCALE:
        status = 'FAIL'
    else:
        status = 'NUMERICALLY_AMBIGUOUS'
    return dict(status=status, candidate_lower_integer=candidate['lower'],
                candidate_upper_integer=candidate['upper'], rank_lower_integer=lower,
                rank_upper_integer=upper, maximum_cost_difference_bound=worst/SCALE,
                minimum_cost_difference_bound=best/SCALE)


def compare_archive(model, skeleton, rule, archive, search):
    need(archive['search_complete'] is True, 'Archived search is incomplete')
    expected_lexical = ([skeleton] if skeleton in model.lexicon else []) if rule == 'identity' else model.lexical_images.get(skeleton, [])
    need(archive['lexicon_candidates'] == expected_lexical, 'Incomplete lexical candidate list')
    saved = archive['top_oov']
    need(len({entry['word'] for entry in saved}) == len(saved), 'Repeated archived OOV word')
    singleton = rule == 'identity' or all(c in VOWELS for c in skeleton)
    if singleton:
        expected = [] if skeleton in model.lexicon else [skeleton]
        need([entry['word'] for entry in saved] == expected, 'Incorrect singleton candidates')
        found = [model.word(word) for word in expected]
    else:
        need(search['search_complete'] and len(search['found']) >= 2, 'Independent search incomplete')
        need(len(saved) == 2, 'Nontrivial infinite OOV support needs two distinct candidates')
        found = search['found']
    ranks, maximum_log_error = [], 0.
    lows, highs = sorted(entry['lower'] for entry in found), sorted(entry['upper'] for entry in found)
    for rank, entry in enumerate(saved):
        word = entry['word']
        need(word not in model.lexicon and (word == skeleton if rule == 'identity' else image(word) == skeleton),
             'Archived word is lexical or outside channel')
        candidate = model.word(word)
        error = abs(candidate['log_char_probability']-entry['log_char_probability'])
        maximum_log_error = max(error, maximum_log_error)
        need(error <= LOG_AGREEMENT_TOLERANCE, 'Archived word log probability differs')
        checked = rank_interval_check(candidate, lows[rank], highs[rank])
        need(checked['status'] != 'FAIL', 'Archived OOV rank is not within fixed tie tolerance')
        ranks.append(dict(rank=rank+1, word=word, **checked))
    result = dict(status='PASS' if all(r['status'] == 'PASS' for r in ranks) else 'NUMERICALLY_AMBIGUOUS',
                  singleton=singleton, lexical_candidates=len(expected_lexical),
                  archived_words=[entry['word'] for entry in saved],
                  independently_first_words=[entry['word'] for entry in found[:2]],
                  identical_word_order=[entry['word'] for entry in saved] == [entry['word'] for entry in found[:2]],
                  maximum_word_log_probability_difference=maximum_log_error, ranks=ranks)
    if not singleton:
        counters = archive['search_counters']
        need(counters['stop_reason'] == 'TOP2_FRONTIER_CERTIFICATE', 'Unexpected archived stopping condition')
        need(counters['expansion_cap'] == 100000 and counters['heap_cap'] == 1000000 and
             counters['expanded_prefixes'] <= 100000 and counters['heap_peak'] <= 1000000,
             'Archived search budget changed or exceeded')
        need(counters['log_cost_tie_tolerance'] == RANK_TOLERANCE, 'Archived tie tolerance changed')
        need(counters['frontier_cost_lower_bound'] >
             counters['second_word_cost_upper_bound']+RANK_TOLERANCE, 'Archived frontier inequality fails')
        intervals = counters['found_top_cost_intervals']
        need([entry['word'] for entry in intervals[:2]] == result['archived_words'],
             'Archived rank interval word order')
        for entry in intervals:
            raw = model.word(entry['word'])['raw_cost']
            need(entry['lower'] <= entry['upper'] and
                 entry['lower']-LOG_AGREEMENT_TOLERANCE <= raw <= entry['upper']+LOG_AGREEMENT_TOLERANCE,
                 'Archived raw-cost interval disagrees with independent word score')
        need(counters['top1_order_strictly_separated'] == (intervals[0]['upper'] < intervals[1]['lower']),
             'Archived top-one separation flag inconsistent')
        need(counters['top2_set_strictly_separated'] ==
             (len(intervals) < 3 or intervals[1]['upper'] < intervals[2]['lower']),
             'Archived top-two separation flag inconsistent')
        gap_low = intervals[1]['lower']-intervals[0]['upper']
        gap_high = intervals[1]['upper']-intervals[0]['lower']
        need(gap_low <= RANK_TOLERANCE if counters['top2_numerical_tie'] else gap_high > RANK_TOLERANCE,
             'Archived tie flag incompatible with archived outward intervals')
        # No transitive/chain clustering: every boundary comparison refers
        # directly to the second order statistic's interval.
        boundary = [entry for entry in found if entry['lower'] <= highs[1]+TIE_UNITS and
                    entry['upper'] >= lows[1]-TIE_UNITS]
        top_tie = highs[1]-lows[0] <= TIE_UNITS
        result.update(independent_top2_numerical_tie=top_tie,
                      archived_top2_numerical_tie=counters['top2_numerical_tie'],
                      second_rank_boundary_words=[entry['word'] for entry in boundary],
                      boundary_identity_ambiguous=len(boundary) > 1,
                      full_search={key: value for key, value in search.items() if key != 'found'},
                      found_terminal_words=len(found))
    return result


def audit(report):
    need(not (ROOT/RESULTS/'computation-failure.json').exists(), 'V2 attempt did not complete')
    integrity = Integrity()
    original = integrity.verify(DATA+'/input-freeze.json')
    prepared = integrity.verify(DATA+'/prepared-manifest.json')
    version = integrity.verify(INPUT)
    terminal = integrity.verify(RESULTS+'/result-freeze.json')
    need(prepared['input_hashes'] == original['hashes'] == version['parent020_input_hashes'], 'Parent hash map')
    need(prepared['input_freeze_sha256'] == version['parent020_input_freeze_sha256'] ==
         sha(DATA+'/input-freeze.json'), 'Parent hash identity')
    need(terminal['hashes'][INPUT] == sha(INPUT), 'Result binds different v2 freeze')
    need(all(version['hashes'].get(path) == value for path, value in
             itertools.chain(original['hashes'].items(), prepared['hashes'].items())), 'Changed parent/preparation')
    report['integrity'] = dict(input_hashes=len(version['hashes']), original_hashes=len(original['hashes']),
                              prepared_outputs=len(prepared['hashes']), result_hashes=len(terminal['hashes']),
                              unique_checked_files=len(integrity.files), hash_checks=HASH_COUNTS['hash_checks'])
    totals = Counter()
    report['languages'] = {}
    for language in ('la', 'de', 'it'):
        model = SearchModel(IndependentModel(language, integrity))
        path = RESULTS+'/'+language+'-completion-certificates.json'
        need(path in terminal['hashes'], 'Completion certificates absent from result freeze')
        archived = read(path)
        observed = set()
        for mechanism, rule in itertools.product(('ordered', 'iid', 'copy'), ('identity', 'omit-aeiou')):
            testpath = f'{DATA}/tests/{language}/{mechanism}-{rule}.json'
            need(testpath in prepared['hashes'], 'Observation file absent from prepared freeze')
            for row in read(testpath)['rows']:
                for token in row['tokens']:
                    totals['observation_slots'] += 1
                    if token['units'] is not None:
                        totals['known_observation_slots'] += 1
                        observed.add(rule+'|'+''.join(token['units']))
        need(observed == set(archived), 'Certificate inventory differs from all observed slots')
        results = {}
        report['languages'][language] = results
        for number, key in enumerate(sorted(observed), 1):
            rule, skeleton = key.split('|', 1)
            singleton = rule == 'identity' or all(c in VOWELS for c in skeleton)
            search = None if singleton else SearchGraph(model, skeleton).enumerate()
            if search is not None and not search['search_complete']:
                results[key] = dict(status='SEARCH_INCOMPLETE', search=search)
                raise RuntimeError('Independent audit budget exhausted: '+language+'/'+key)
            results[key] = compare_archive(model, skeleton, rule, archived[key], search)
            totals['analyses'] += 1
            totals['singletons' if singleton else 'full_prefix_searches'] += 1
            totals['word_order_differences'] += not results[key]['identical_word_order']
            totals['rank_interval_ambiguities'] += results[key]['status'] != 'PASS'
            if search is not None:
                totals['expanded_prefixes'] += search['expanded_prefixes']
                totals['lexicon_eos_rejected'] += search['lexicon_eos_rejected']
            if number % 100 == 0:
                print(json.dumps(dict(language=language, checked=number, total=len(observed))), file=sys.stderr, flush=True)
        print(json.dumps(dict(language=language, status='AUDITED', analyses=len(results))), file=sys.stderr, flush=True)
    report['totals'] = dict(totals)
    report['status'] = 'PASS' if totals['rank_interval_ambiguities'] == 0 else 'NUMERICALLY_AMBIGUOUS'


def self_check():
    """Pure toy search versus finite brute force with an explicit tail bound."""
    from fractions import Fraction
    from types import SimpleNamespace
    checks = 0
    for value in (0., 1e-300, 1e-20, .1, 1., 17.3):
        lower, upper = enclosure(value)
        exact = Fraction(*value.as_integer_ratio())
        need(Fraction(lower, SCALE) <= exact <= Fraction(upper, SCALE), 'Dyadic enclosure')
        checks += 1
    for probabilities in ((.25, .25, .25, .25), (.1, .2, .3, .4)):
        alphabet = ' abx'
        tables = [np.log(np.tile(probabilities, 4**(n-1))).reshape((4,)*n) for n in range(1, 5)]
        fake = SimpleNamespace(alphabet=alphabet, words={'b', 'ab', 'ba'}, index={c:i for i,c in enumerate(alphabet)},
                               tables=tables, log_nonempty=math.log(1-probabilities[0]))
        model = SearchModel(fake)
        result = SearchGraph(model, 'b').enumerate()
        need(result['search_complete'], 'Toy search incomplete')
        brute = []
        for length in range(1, 11):
            for b_position in range(length):
                word = 'a'*b_position+'b'+'a'*(length-b_position-1)
                if word not in fake.words:
                    brute.append(model.word(word))
        brute.sort(key=lambda entry: (entry['raw_cost'], entry['word']))
        # Every missing word is >=11 letters. Each letter probability <=.4;
        # including EOS only raises its positive cost. This is a valid tail bound.
        need(11*-math.log(max(probabilities[1:])) > brute[1]['raw_cost']+RANK_TOLERANCE, 'Unproven brute tail')
        for rank in (0, 1):
            need(abs(result['found'][rank]['raw_cost']-brute[rank]['raw_cost']) < 1e-14, 'Toy top rank')
            checks += 1
        need(result['lexicon_eos_rejected'] >= 3, 'Lexicon EOS exclusion not exercised')
        need(len([entry for entry in result['found'] if
                  abs(entry['raw_cost']-brute[1]['raw_cost']) <= RANK_TOLERANCE]) >= 3,
             'Boundary ties not exhausted')
        need(not SearchGraph(model, 'b').enumerate(expansion_cap=0)['search_complete'], 'Budget failure accepted')
        checks += 3
    sample = dict(lower=100, upper=101)
    need(rank_interval_check(sample, 100, 101)['status'] == 'PASS', 'Rank enclosure')
    need(rank_interval_check(dict(lower=SCALE, upper=SCALE), 0, 0)['status'] == 'FAIL', 'Bad rank accepted')
    print(json.dumps(dict(status='PASS', scope='SYNTHETIC_ONLY', checks=checks+2,
                          dyadic_scale=SCALE, production_imports=False)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-check', action='store_true')
    parser.add_argument('--output', default=OUTPUT)
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    report = dict(status='RUNNING_NOT_COMPLETE', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  audit_sha256=sha('src/audit_omission_search_v2.py'),
                  dependency=dict(file='src/audit_omission_recovery_v2.py',
                                  sha256=sha('src/audit_omission_recovery_v2.py'),
                                  used='IndependentModel loader/normalizer and Integrity hash verifier only'),
                  method='Lexicon-free finite graph, exact-integer reverse Dijkstra, distinct-prefix A*; lexical EOS rejection',
                  rank_tolerance=RANK_TOLERANCE, log_probability_agreement_tolerance=LOG_AGREEMENT_TOLERANCE,
                  expansion_cap=EXPANSION_CAP, heap_cap=HEAP_CAP, dyadic_scale=SCALE,
                  scope='Independent represented-table global top-two/tie-frontier search; no mass or HMM claim')
    exit_code = 0
    try:
        audit(report)
        if report['status'] != 'PASS':
            exit_code = 1
    except Exception as error:
        report.update(status='FAIL_OR_INCOMPLETE', exception=type(error).__name__, reason=str(error))
        traceback.print_exc()
        exit_code = 1
    report['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    destination = ROOT/args.output
    need(not destination.exists(), 'Refusing to overwrite a previous audit output')
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({key: report[key] for key in ('status', 'totals', 'reason') if key in report}), flush=True)
    sys.exit(exit_code)


if __name__ == '__main__':
    main()
