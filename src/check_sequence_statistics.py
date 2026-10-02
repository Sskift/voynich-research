"""Pure synthetic independent checks for experiment019 statistics.

This file never reads prepared data, source text, manuscript rows, or result files.
The reference uses tuple-valued categories, direct probability arithmetic, explicit
segment/pair enumeration, and ordinary Levenshtein dynamic programming. Production
functions are called only on the synthetic objects constructed below.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
import json
import math

import numpy as np
import sequence_dependence as production


METRICS = ('word_gain', 'edge_gain', 'local_exact1')
SPLITS = ('train', 'development', 'reserve')
EPS = 1e-12
CHECKS = Counter()
MAX_ERROR = 0.0


def equal(actual, expected, label='value'):
    global MAX_ERROR
    if isinstance(expected, dict):
        assert set(actual) == set(expected), (label, set(actual), set(expected))
        for key in expected:
            equal(actual[key], expected[key], f'{label}.{key}')
    elif isinstance(expected, (tuple, list)):
        assert len(actual) == len(expected), (label, len(actual), len(expected))
        for index, (a, e) in enumerate(zip(actual, expected)):
            equal(a, e, f'{label}[{index}]')
    elif isinstance(expected, float):
        error = abs(float(actual)-expected)
        MAX_ERROR = max(MAX_ERROR, error)
        assert math.isfinite(float(actual)) and error <= 1e-10, (label, actual, expected)
    else:
        assert actual == expected, (label, actual, expected)


def row(page, locus, split, words, resets=(), paragraph=False):
    return dict(page=page, locus=locus, split=split, paragraph_first_line=paragraph,
                tokens=[dict(units=None if w is None else list(w), reset_before=i in resets)
                        for i, w in enumerate(words)])


def reference_layout(data, split):
    rows = [r for r in data['rows'] if r['split'] == split]
    pages = list(dict.fromkeys(r['page'] for r in rows))
    slots, row_ranges = [], []
    for row_number, r in enumerate(rows):
        indices = []
        for k, token in enumerate(r['tokens']):
            size = len(r['tokens'])
            role = 'only' if size == 1 else 'first' if k == 0 else 'last' if k == size-1 else 'middle'
            indices.append(len(slots))
            slots.append(dict(word=None if token['units'] is None else tuple(token['units']),
                              reset=token['reset_before'], page=r['page'], row=row_number,
                              role=role, paragraph=r['paragraph_first_line']))
        row_ranges.append(indices)
    segments, current = [], []
    for i, slot in enumerate(slots):
        interrupted = (slot['reset'] or (i > 0 and slot['page'] != slots[i-1]['page'])
                       or slot['word'] is None)
        if interrupted and current:
            segments.append(current)
            current = []
        if slot['word'] is not None:
            current.append(i)
    if current:
        segments.append(current)
    starts, ends = set(), set()
    local = []
    for segment in segments:
        starts.add(segment[0]); ends.add(segment[-1])
        for position, right in enumerate(segment):
            for left in segment[max(0, position-8):position]:
                local.append((left, right))
    pairs = []
    for indices in row_ranges:
        for a, b in zip(indices, indices[1:]):
            if slots[a]['word'] is not None and slots[b]['word'] is not None and not slots[b]['reset']:
                pairs.append((a, b))
    edge_pairs = [(a, b) for a, b in pairs if min(len(slots[a]['word']), len(slots[b]['word'])) >= 2]
    strata = defaultdict(list)
    for i, slot in enumerate(slots):
        if slot['word'] is not None:
            key = (slot['page'], slot['paragraph'], slot['role'], len(slot['word']), i in starts, i in ends)
            strata[key].append(i)
    groups = [strata[k] for k in sorted(strata)]
    words = [s['word'] for s in slots]
    return dict(slots=slots, pages=pages, words=words, starts=starts, ends=ends,
                pairs=(pairs, edge_pairs), local=local, groups=groups,
                mobile=sum(len(g) for g in groups if len({words[i] for i in g}) > 1),
                known=sum(w is not None for w in words), unknown=sum(w is None for w in words))


def categories(left, right, edge):
    return (left[-2:], right[:2]) if edge else (left, right)


def reference_models(train):
    models = []
    for edge in (False, True):
        pairs = [categories(train['words'][a], train['words'][b], edge) for a, b in train['pairs'][edge]]
        models.append(dict(edge=edge, targets=Counter(v for _, v in pairs),
                           contexts=Counter(u for u, _ in pairs), pairs=Counter(pairs), N=len(pairs),
                           tau=10 if edge else 100))
    return models


def probability_gain(model, left, right):
    u, v = categories(left, right, model['edge'])
    V = len(model['targets'])
    p0 = (model['targets'][v] + .5) / (model['N'] + .5*(V+1))
    if not model['contexts'][u]:
        return 0.0
    p1 = (model['pairs'][u, v] + model['tau']*p0) / (model['contexts'][u] + model['tau'])
    return math.log(p1/p0)


@lru_cache(maxsize=None)
def levenshtein(left, right):
    last = list(range(len(right)+1))
    for i, a in enumerate(left, 1):
        current = [i]
        for j, b in enumerate(right, 1):
            current.append(min(last[j]+1, current[-1]+1, last[j-1]+int(a != b)))
        last = current
    return last[-1]


def reference_evaluate(layout, models, words=None):
    words = layout['words'] if words is None else words
    sums = [[0.0 for _ in layout['pages']] for _ in METRICS]
    counts = [[0 for _ in layout['pages']] for _ in METRICS]
    page_ids = {page: i for i, page in enumerate(layout['pages'])}
    coverages = {}
    for k, model in enumerate(models):
        coverage = dict(pairs=len(layout['pairs'][k]), target_seen=0, context_seen=0, both_seen=0, pair_seen=0)
        for a, b in layout['pairs'][k]:
            page = page_ids[layout['slots'][b]['page']]
            left, right = words[a], words[b]
            sums[k][page] += probability_gain(model, left, right)
            counts[k][page] += 1
            u, v = categories(left, right, model['edge'])
            t, c = v in model['targets'], u in model['contexts']
            coverage['target_seen'] += t; coverage['context_seen'] += c
            coverage['both_seen'] += t and c
            coverage['pair_seen'] += (u, v) in model['pairs']
        coverages[METRICS[k]] = coverage
    zero, one = 0, 0
    for a, b in layout['local']:
        page = page_ids[layout['slots'][b]['page']]
        distance = levenshtein(words[a], words[b])
        zero += distance == 0; one += distance == 1
        sums[2][page] += distance == 1; counts[2][page] += 1
    extra = dict(prediction_coverage=coverages, local_pairs=len(layout['local']),
                 local_exact0_pairs=zero, local_exact1_pairs=one, local_at_most1_pairs=zero+one)
    return np.array(sums), np.array(counts, dtype=np.int64), extra


def reference_sample(layout, models, seeds, reps, boots):
    observed, counts, extra = reference_evaluate(layout, models)
    totals = counts.sum(axis=1)
    nulls, page_values = [], np.zeros_like(observed)
    rng = np.random.default_rng(seeds['shuffle'])
    words = layout['words']
    for _ in range(reps):
        shuffled = list(words)
        for group in layout['groups']:
            # Permute local integer indices so no tuple word can be interpreted
            # by numpy as an extra array dimension. This consumes identical RNG.
            order = rng.permutation(len(group))
            for to, source in zip(group, order):
                shuffled[to] = words[group[int(source)]]
        sums, checked_counts, _ = reference_evaluate(layout, models, shuffled)
        assert np.array_equal(checked_counts, counts)
        page_values += sums
        nulls.append([float(sums[k].sum()/n) if n else None for k, n in enumerate(totals)])
    page_means = page_values / reps
    npages = len(layout['pages'])
    draws = np.random.default_rng(seeds['bootstrap']).integers(npages, size=(boots, npages))
    result = {}
    for k, metric in enumerate(METRICS):
        n = int(totals[k])
        score = float(observed[k].sum()/n) if n else None
        values = [v[k] for v in nulls]
        mean = sum(values)/reps if n else None
        sd = math.sqrt(sum((v-mean)**2 for v in values)/(reps-1)) if n else None
        excess = score-mean if n else None
        informative = bool(n and layout['mobile'] and max(values)-min(values) > EPS)
        p = (1+sum(v >= score-EPS for v in values))/(reps+1) if n else None
        denominators = [sum(int(counts[k, p]) for p in draw) for draw in draws]
        zero_draws = denominators.count(0)
        raw_ci = excess_ci = None
        if zero_draws == 0:
            raw = [sum(float(observed[k, p]) for p in draw)/d for draw, d in zip(draws, denominators)]
            exc = [sum(float(observed[k, p]-page_means[k, p]) for p in draw)/d
                   for draw, d in zip(draws, denominators)]
            raw_ci = np.quantile(raw, [.025, .975], method='linear').tolist()
            excess_ci = np.quantile(exc, [.025, .975], method='linear').tolist()
        result[metric] = dict(observed=score, sum=float(observed[k].sum()), comparisons=n,
                             null_mean=mean, null_sd=sd, excess=excess, p_greater=p,
                             informative=informative,
                             positive=bool(informative and p <= .05/3 and (score > EPS if k < 2 else excess > EPS)),
                             threshold=.05/3,
                             mobile_fraction=layout['mobile']/layout['known'] if layout['known'] else None,
                             null_statistics=values, raw_page_bootstrap_ci=raw_ci,
                             excess_page_bootstrap_ci=excess_ci,
                             zero_denominator_bootstrap_draws=zero_draws,
                             pages=[dict(page=page, sum=float(observed[k, i]), count=int(counts[k, i]),
                                         null_mean_sum=float(page_means[k, i])) for i, page in enumerate(layout['pages'])])
    return dict(metrics=result, **extra)


def inspect_layout(corpus, data, split):
    layout = reference_layout(data, split)
    actual = corpus.splits[split]
    equal(actual['page_names'], layout['pages'], 'page_order')
    equal(np.flatnonzero(actual['starts']).tolist(), sorted(layout['starts']), 'segment_starts')
    equal(np.flatnonzero(actual['ends']).tolist(), sorted(layout['ends']), 'segment_ends')
    equal([g.tolist() for g in actual['groups']], layout['groups'], 'strata')
    for i in range(2):
        equal(list(zip(*(a.tolist() for a in actual['pair_indices'][i]))), layout['pairs'][i], 'prediction_pairs')
    equal(list(zip(*(a.tolist() for a in actual['local_indices']))), layout['local'], 'local_pairs')
    for field in ('mobile', 'known', 'unknown'):
        equal(actual[field], layout[field], field)
    CHECKS['layouts'] += 1
    return layout


def inspect_case(data, sample=True):
    corpus = production.Corpus(data)
    layouts = {s: inspect_layout(corpus, data, s) for s in SPLITS}
    models = reference_models(layouts['train'])
    for k, model in enumerate(models):
        summary = corpus.models[k].summary()
        equal(summary['training_pairs'], model['N'])
        equal(summary['training_target_types'], len(model['targets']))
        equal(summary['training_context_types'], len(model['contexts']))
        equal(summary['training_pair_types'], len(model['pairs']))
        # All categories are deliberately queried, including contexts or targets
        # seen only in held-out rows and targets present only as train contexts.
        for left in range(len(corpus.words)):
            right = np.arange(len(corpus.words), dtype=np.int64)
            gains, _ = corpus.models[k].score(np.full(len(right), left, dtype=np.int64), right)
            for j, gain in enumerate(gains):
                equal(float(gain), probability_gain(model, corpus.words[left], corpus.words[j]), 'probability_gain')
                CHECKS['individual_probability_gains'] += 1
    for split in ('development', 'reserve'):
        actual_sums, actual_counts, actual_extra = corpus.evaluate(split)
        sums, counts, extra = reference_evaluate(layouts[split], models)
        equal(actual_sums.tolist(), sums.tolist(), 'page_sums')
        equal(actual_counts.tolist(), counts.tolist(), 'page_counts')
        equal(actual_extra, extra, 'coverage_and_local')
        if sample:
            seeds = dict(shuffle=731 if split == 'development' else 947, bootstrap=501031)
            observed = production.sample_test(corpus, split, seeds)
            expected = reference_sample(layouts[split], models, seeds, production.REPS, production.BOOTS)
            equal({key: observed[key] for key in expected}, expected, f'full_sample_{split}')
            CHECKS['full_sample_tests'] += 1
            CHECKS['independent_null_metrics'] += 3*production.REPS
            CHECKS['independent_bootstrap_metric_draws'] += 3*production.BOOTS
    return corpus


def fixtures():
    # Pages deliberately appear in nonlexicographic order. Consecutive row
    # fragments include known->None, None->known, internal reset, page change,
    # empty rows and a one-word row. Multi-atom labels must remain indivisible.
    return [dict(rows=[
        row('train-z', 't1', 'train', ['aa', 'ab', 'aa', 'bb', 'ba', 'ab'], (0,)),
        row('train-z', 't2', 'train', ['q', 'aa', None, 'zz', 'ab'], (0, 4)),
        row('z', 'd1', 'development', ['aa', 'ab', 'ba', 'bb', 'ab'], (0,), True),
        row('z', 'd2', 'development', ['bb', 'ba', 'aa', 'ab', 'bb']),
        row('z', 'd3', 'development', ['ab', None, None, 'ba', 'aa', 'bb'], (4,)),
        row('z', 'd4', 'development', ['ba']),
        row('a', 'd1', 'development', ['xx', 'aa', 'xy', 'bb', 'a'], (), True),
        row('a', 'd2', 'development', [('A', '@042;'), ('a', '@042;'), 'aa', 'bb', 'xy']),
        row('a', 'd3', 'development', [], ()),
        row('a', 'd4', 'development', ['aa', 'ab', 'ba', 'bb', 'zz'], (0,)),
        row('r', 'r1', 'reserve', ['ba', 'aa', 'ab', 'bb'], (0,)),
        row('r', 'r2', 'reserve', ['ab', 'bb', 'aa', 'ba']),
    ]), dict(rows=[
        row('t', 't1', 'train', ['aa', 'bb'], (0,)),
        row('z', 'd1', 'development', ['aa', 'ab', 'bb'], (0,)),
        row('a', 'd1', 'development', [None], (0,)),
        row('r', 'r1', 'reserve', [None, None], (0,)),
    ]), dict(rows=[])]


def graph_check():
    alphabet = ('a', 'A', '@042;')
    words = [w for length in range(1, 5) for w in product(alphabet, repeat=length)]
    data = dict(rows=[row('toy', 'd1', 'development', words, (0,))])
    corpus = production.Corpus(data)
    layout = inspect_layout(corpus, data, 'development')
    for a, left in enumerate(corpus.words):
        for b, right in enumerate(corpus.words):
            actual = bool(corpus.near[corpus.local_lookup[a], corpus.local_lookup[b]])
            equal(actual, levenshtein(left, right) == 1, 'exact_edit_graph')
            CHECKS['graph_pairs_checked_by_dp'] += 1
    equal(len(layout['local']), sum(min(i, 8) for i in range(len(words))), 'window_eight')
    assert not corpus.near.diagonal().any()


def bijection_check(data):
    encoded = dict(rows=[])
    atoms = sorted({a for r in data['rows'] for t in r['tokens'] if t['units'] is not None for a in t['units']})
    mapping = {a: f'opaque-code-{i}' for i, a in enumerate(reversed(atoms))}
    for r in data['rows']:
        transformed = {key: value for key, value in r.items() if key != 'tokens'}
        transformed['tokens'] = [dict(units=None if t['units'] is None else [mapping[a] for a in t['units']],
                                      reset_before=t['reset_before']) for t in r['tokens']]
        encoded['rows'].append(transformed)
    first, second = production.Corpus(data), production.Corpus(encoded)
    for split in ('development', 'reserve'):
        seeds = dict(shuffle=811, bootstrap=500811)
        equal(production.sample_test(first, split, seeds), production.sample_test(second, split, seeds), 'bijection')
        CHECKS['paired_bijection_tests'] += 1


def main():
    for fixture in fixtures():
        inspect_case(fixture)
    graph_check()
    bijection_check(fixtures()[0])
    print(json.dumps(dict(status='PASS', scope='SYNTHETIC_ONLY_NO_SOURCE_OR_MANUSCRIPT_READS',
                          checks=dict(CHECKS), maximum_absolute_error=MAX_ERROR,
                          permutations_per_test=production.REPS, bootstrap_resamples_per_test=production.BOOTS), indent=2))


if __name__ == '__main__':
    main()
