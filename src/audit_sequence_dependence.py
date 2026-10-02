"""Independent experiment019 statistics audit; never imports either runner.

Reconstructs layout/segments with explicit segment lists, training counts with
tuple-keyed Counters, and local EXACT1 relations using a trie Levenshtein dynamic
program. Replays all frozen index permutations and page bootstrap draws. Source
normalization, mechanism generation and encoding are covered by the separate
input audit. This program requires existing results before evaluating any case.
Only stdout JSON and stderr progress are written; prepared/frozen files are read.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LANGS = ('la', 'de', 'it')
MECHANISMS = ('ordered', 'iid', 'copy')
CHANNELS = ('bijective', 'uniform127', 'position127', 'hidden-space127')
SPLITS = ('train', 'development', 'reserve')
METRICS = ('word_gain', 'edge_gain', 'local_exact1')
REPS, BOOTS, TOL = 999, 1999, 1e-12
CHECKS = Counter()
MAX_ERROR = 0.0


def read(path):
    return json.loads((ROOT/path).read_text())


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def value_sha(value):
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, separators=(',', ':'),
                                     sort_keys=True).encode()).hexdigest()


def require(condition, detail):
    if not condition:
        raise AssertionError(detail)


def compare(actual, expected, path='result'):
    global MAX_ERROR
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and set(actual) == set(expected),
                (path, 'dictionary keys', set(actual), set(expected)))
        for key in expected:
            compare(actual[key], expected[key], f'{path}.{key}')
    elif isinstance(expected, (tuple, list)):
        require(len(actual) == len(expected), (path, 'length', len(actual), len(expected)))
        for i, (left, right) in enumerate(zip(actual, expected)):
            compare(left, right, f'{path}[{i}]')
    elif isinstance(expected, float):
        error = abs(float(actual)-expected)
        MAX_ERROR = max(MAX_ERROR, error)
        # Gains/rates receive a stricter tolerance than unnormalized page sums.
        tolerance = 2e-9 if path.endswith('.sum') or path.endswith('.null_mean_sum') else 1e-10
        require(math.isfinite(float(actual)) and error <= tolerance, (path, actual, expected, error))
        CHECKS['floating_values'] += 1
    else:
        require(actual == expected, (path, actual, expected))
        CHECKS['exact_values'] += 1


def verify_manifest(path):
    manifest = read(path)
    require(bool(manifest.get('hashes')), (path, 'empty hashes'))
    for filename, expected in manifest['hashes'].items():
        require(sha(filename) == expected, (path, filename, 'hash mismatch'))
        CHECKS['frozen_hashes'] += 1
    return manifest


def expected_cases():
    return [f'{lang}/{mechanism}-{channel}' for lang in LANGS
            for mechanism in MECHANISMS for channel in ('plain',)+CHANNELS]


def seed_pair(case, split):
    split_index = ('development', 'reserve').index(split)
    if case == 'native':
        seed = 20261201+900000+10*split_index
    else:
        lang, suffix = case.split('/')
        mechanism, channel = next((m, suffix[len(m)+1:]) for m in MECHANISMS if suffix.startswith(m+'-'))
        channel_index = 0 if channel == 'plain' else CHANNELS.index(channel)
        seed = (20261201+100000*LANGS.index(lang)+30000+1000*MECHANISMS.index(mechanism)
                +100*channel_index+10*split_index)
    return dict(shuffle=seed, bootstrap=seed+500000)


def verify_inputs():
    original = verify_manifest('data/sequence-dependence/input-freeze.json')
    prepared = verify_manifest('data/sequence-dependence/prepared-manifest.json')
    compare(prepared['input_hashes'], original['hashes'], 'input_hashes')
    compare(prepared['input_freeze_sha256'], sha('data/sequence-dependence/input-freeze.json'))
    paths = {f'data/sequence-dependence/{case}.json' for case in expected_cases()+['native']}
    paths |= {'data/sequence-dependence/config.json', 'data/sequence-dependence/template.json'}
    compare(set(prepared['hashes']), paths, 'prepared_inventory')
    config = read('data/sequence-dependence/config.json')
    fixed = dict(languages=list(LANGS), mechanisms=list(MECHANISMS), channels=list(CHANNELS),
                 splits=list(SPLITS), test_splits=['development', 'reserve'], primary_split='development',
                 tau_word=100, tau_edge=10, smoothing=.5, local_window=8, local_primary_distance=1,
                 permutations=REPS, bootstrap_replicates=BOOTS, bootstrap_confidence=.95,
                 familywise_alpha=.05, bonferroni_tests=3,
                 strata=['page', 'paragraph_first_line', 'row_role', 'exact_atom_length', 'segment_start', 'segment_end'],
                 singleton_permutation=True, rng='numpy.random.default_rng')
    for key, expected in fixed.items():
        compare(config[key], expected, f'config.{key}')
    require(config['numpy_version'] == np.__version__, 'Frozen numpy version differs')
    expected_seeds = {case: {s: seed_pair(case, s) for s in SPLITS[1:]} for case in expected_cases()}
    compare(config['seeds']['cases'], expected_seeds, 'case_seeds')
    compare(config['seeds']['native'], {s: seed_pair('native', s) for s in SPLITS[1:]}, 'native_seeds')
    return dict(input_hashes=len(original['hashes']), prepared_hashes=len(prepared['hashes']),
                input_freeze_sha256=sha('data/sequence-dependence/input-freeze.json'),
                prepared_manifest_sha256=sha('data/sequence-dependence/prepared-manifest.json'))


def layout(rows, split):
    selected = [r for r in rows if r['split'] == split]
    pages = list(dict.fromkeys(r['page'] for r in selected))
    page_number = {page: i for i, page in enumerate(pages)}
    slots, row_indices = [], []
    for r in selected:
        require(r['page'] not in {'f49v', 'f58r', 'f58v'}, 'Forbidden manuscript page')
        indices = []
        for j, token in enumerate(r['tokens']):
            word = None if token['units'] is None else tuple(token['units'])
            require(word is None or (word and all(isinstance(a, str) and a for a in word)), 'Invalid observed word')
            count = len(r['tokens'])
            role = 'only' if count == 1 else 'first' if j == 0 else 'last' if j == count-1 else 'middle'
            indices.append(len(slots))
            slots.append(dict(word=word, page=r['page'], reset=token['reset_before'],
                              paragraph=r['paragraph_first_line'], role=role))
        row_indices.append(indices)
    segments, segment = [], []
    for i, slot in enumerate(slots):
        if slot['reset'] or slot['word'] is None or (i and slot['page'] != slots[i-1]['page']):
            if segment:
                segments.append(segment)
                segment = []
        if slot['word'] is not None:
            segment.append(i)
    if segment:
        segments.append(segment)
    first = {s[0] for s in segments}; last = {s[-1] for s in segments}
    local = [(segment[i-lag], index) for segment in segments for i, index in enumerate(segment)
             for lag in range(min(8, i), 0, -1)]
    word_pairs = [(a, b) for indices in row_indices for a, b in zip(indices, indices[1:])
                  if slots[a]['word'] is not None and slots[b]['word'] is not None and not slots[b]['reset']]
    edge_pairs = [(a, b) for a, b in word_pairs if min(len(slots[a]['word']), len(slots[b]['word'])) >= 2]
    grouped = defaultdict(list)
    for i, slot in enumerate(slots):
        if slot['word'] is not None:
            grouped[(slot['page'], slot['paragraph'], slot['role'], len(slot['word']), i in first, i in last)].append(i)
    groups = [np.array(grouped[key], dtype=np.int64) for key in sorted(grouped)]
    words = [s['word'] for s in slots]
    arrays = [np.array(pairs, dtype=np.int64).reshape(-1, 2) for pairs in (word_pairs, edge_pairs, local)]
    return dict(pages=pages, words=words, groups=groups, segments=segments,
                pairs=arrays, page=np.array([page_number[s['page']] for s in slots], dtype=np.int64),
                known=sum(w is not None for w in words), unknown=sum(w is None for w in words),
                mobile=sum(len(g) for g in groups if len({words[int(i)] for i in g}) > 1))


def category(left, right, edge):
    return (left[-2:], right[:2]) if edge else (left, right)


class Model:
    def __init__(self, train, edge, vocabulary):
        self.edge, self.tau = edge, 10 if edge else 100
        pairs = Counter(category(train['words'][a], train['words'][b], edge)
                        for a, b in train['pairs'][int(edge)])
        contexts, targets = Counter(), Counter()
        for (u, v), count in pairs.items():
            contexts[u] += count; targets[v] += count
        n = sum(pairs.values()); require(n == sum(contexts.values()) == sum(targets.values()), 'Count marginal mismatch')
        self.summary = dict(training_pairs=n, training_target_types=len(targets),
                            training_context_types=len(contexts), training_pair_types=len(pairs),
                            unigram_pseudocount=.5, backoff_mass=self.tau,
                            unseen_target='one UNK bucket', unseen_context='exact P0 backoff')
        self.count_hash = value_sha(dict(targets=sorted(targets.items()), contexts=sorted(contexts.items()),
                                         pairs=sorted((u, v, count) for (u, v), count in pairs.items())))
        # Reserve index zero for a single UNK target / unseen-context bucket.
        context_ids = {u: i+1 for i, u in enumerate(sorted(contexts))}
        target_ids = {v: i+1 for i, v in enumerate(sorted(targets))}
        self.width = len(target_ids)+1
        self.left = np.array([context_ids.get(w[-2:] if edge else w, 0) for w in vocabulary], dtype=np.int64)
        self.right = np.array([target_ids.get(w[:2] if edge else w, 0) for w in vocabulary], dtype=np.int64)
        context_count = np.array([0]+[contexts[u] for u in sorted(contexts)], dtype=float)
        target_count = np.array([0]+[targets[v] for v in sorted(targets)], dtype=float)
        self.p0 = (target_count+.5)/(n+.5*(len(targets)+1))
        self.base_gain = np.log(self.tau/(context_count+self.tau))
        entries = sorted((context_ids[u]*self.width+target_ids[v], count) for (u, v), count in pairs.items())
        self.keys = np.array([key for key, _ in entries], dtype=np.int64)
        self.values = np.array([count for _, count in entries], dtype=float)
        CHECKS['training_pair_entries'] += len(pairs)
        CHECKS['training_target_entries'] += len(targets)
        CHECKS['training_context_entries'] += len(contexts)

    def query(self, left_ids, right_ids):
        u, v = self.left[left_ids], self.right[right_ids]
        flat = u*self.width+v
        pair_counts = np.zeros(len(u))
        if len(self.keys):
            index = np.searchsorted(self.keys, flat)
            valid = index < len(self.keys)
            valid[valid] &= self.keys[index[valid]] == flat[valid]
            pair_counts[valid] = self.values[index[valid]]
        # Independent algebraic decomposition: context penalty plus positive
        # observed-pair correction. UNK probability is represented only once.
        gains = self.base_gain[u] + np.log1p(pair_counts/(self.tau*self.p0[v]))
        coverage = dict(pairs=len(u), target_seen=int(np.count_nonzero(v)), context_seen=int(np.count_nonzero(u)),
                        both_seen=int(np.count_nonzero((u != 0) & (v != 0))),
                        pair_seen=int(np.count_nonzero(pair_counts)))
        return gains, coverage


def trie_neighbors(vocabulary):
    """All directed EXACT1 pairs via trie DP, without deletion signatures."""
    trie = {}
    for index, word in enumerate(vocabulary):
        node = trie
        for atom in word:
            node = node.setdefault(atom, {})
        node[None] = index
    relation = []
    width = len(vocabulary)
    for source, word in enumerate(vocabulary):
        stack = [(trie, list(range(len(word)+1)))]
        while stack:
            node, old = stack.pop()
            if None in node and old[-1] == 1:
                relation.append(source*width+node[None])
            for atom, child in node.items():
                if atom is None:
                    continue
                new = [old[0]+1]
                for j, target in enumerate(word, 1):
                    new.append(min(new[-1]+1, old[j]+1, old[j-1]+int(atom != target)))
                if min(new) <= 1:
                    stack.append((child, new))
    CHECKS['directed_exact1_relations'] += len(relation)
    return np.array(sorted(relation), dtype=np.int64)


class Case:
    def __init__(self, data):
        self.layouts = {s: layout(data['rows'], s) for s in SPLITS}
        vocabulary = sorted({w for lay in self.layouts.values() for w in lay['words'] if w is not None})
        ids = {w: i for i, w in enumerate(vocabulary)}
        for lay in self.layouts.values():
            lay['ids'] = np.array([-1 if w is None else ids[w] for w in lay['words']], dtype=np.int64)
        self.models = [Model(self.layouts['train'], edge, vocabulary) for edge in (False, True)]
        # Restrict trie to held-out observed words; no train-only local pairs
        # can occur in observed data or within held-out index permutations.
        visible = sorted({w for s in SPLITS[1:] for w in self.layouts[s]['words'] if w is not None})
        local_ids = {w: i for i, w in enumerate(visible)}
        self.local = np.array([local_ids.get(w, -1) for w in vocabulary], dtype=np.int64)
        self.width = len(visible)
        self.near = trie_neighbors(visible)

    def evaluate(self, split, ids=None):
        lay = self.layouts[split]
        ids = lay['ids'] if ids is None else ids
        sums = np.zeros((3, len(lay['pages'])))
        counts = np.zeros_like(sums, dtype=np.int64)
        coverage = {}
        for k, (indices, model) in enumerate(zip(lay['pairs'], self.models)):
            scores, coverage[METRICS[k]] = model.query(ids[indices[:, 0]], ids[indices[:, 1]])
            for page in range(len(lay['pages'])):
                selected = lay['page'][indices[:, 1]] == page
                sums[k, page] = scores[selected].sum()
                counts[k, page] = np.count_nonzero(selected)
        a, b = lay['pairs'][2].T
        left, right = self.local[ids[a]], self.local[ids[b]]
        require(np.all(left >= 0) and np.all(right >= 0), 'Local pair contains an unknown slot')
        exact1 = np.isin(left*self.width+right, self.near, assume_unique=False, kind='sort')
        exact0 = ids[a] == ids[b]
        for page in range(len(lay['pages'])):
            selected = lay['page'][b] == page
            sums[2, page] = np.count_nonzero(exact1[selected])
            counts[2, page] = np.count_nonzero(selected)
        extras = dict(prediction_coverage=coverage, local_pairs=len(a),
                      local_exact0_pairs=int(np.count_nonzero(exact0)),
                      local_exact1_pairs=int(np.count_nonzero(exact1)),
                      local_at_most1_pairs=int(np.count_nonzero(exact0 | exact1)))
        return sums, counts, extras

    def test(self, name, split):
        lay = self.layouts[split]
        seeds = seed_pair(name, split)
        observed, counts, extras = self.evaluate(split)
        totals = counts.sum(axis=1)
        nulls = np.empty((REPS, 3))
        page_sum = np.zeros_like(observed)
        random = np.random.default_rng(seeds['shuffle'])
        for rep in range(REPS):
            shuffled = lay['ids'].copy()
            for group in lay['groups']:
                order = random.permutation(len(group))
                shuffled[group] = lay['ids'][group[order]]
            sums, permuted_counts, _ = self.evaluate(split, shuffled)
            require(np.array_equal(permuted_counts, counts), 'Permutation changed eligible comparisons')
            page_sum += sums
            nulls[rep] = np.divide(sums.sum(axis=1), totals, out=np.full(3, np.nan), where=totals != 0)
        page_mean = page_sum/REPS
        draw = np.random.default_rng(seeds['bootstrap']).integers(len(lay['pages']), size=(BOOTS, len(lay['pages'])))
        metrics = {}
        for k, metric in enumerate(METRICS):
            n = int(totals[k]); score = float(observed[k].sum()/n) if n else None
            mean = float(nulls[:, k].mean()) if n else None
            sd = float(nulls[:, k].std(ddof=1)) if n else None
            excess = score-mean if n else None
            p = (1+int(np.count_nonzero(nulls[:, k] >= score-TOL)))/(REPS+1) if n else None
            informative = bool(n and lay['mobile'] and float(np.max(nulls[:, k])-np.min(nulls[:, k])) > TOL)
            positive = bool(informative and p <= .05/3 and (score > TOL if k < 2 else excess > TOL))
            denominators = counts[k, draw].sum(axis=1)
            zero_draws = int(np.count_nonzero(denominators == 0))
            raw_ci = excess_ci = None
            if not zero_draws:
                raw = observed[k, draw].sum(axis=1)/denominators
                extra = (observed[k]-page_mean[k])[draw].sum(axis=1)/denominators
                raw_ci = np.quantile(raw, [.025, .975], method='linear').tolist()
                excess_ci = np.quantile(extra, [.025, .975], method='linear').tolist()
            metrics[metric] = dict(observed=score, sum=float(observed[k].sum()), comparisons=n,
                                  null_mean=mean, null_sd=sd, excess=excess, p_greater=p,
                                  informative=informative, positive=positive, threshold=.05/3,
                                  mobile_fraction=lay['mobile']/lay['known'] if lay['known'] else None,
                                  null_statistics=nulls[:, k].tolist() if n else [None]*REPS,
                                  raw_page_bootstrap_ci=raw_ci, excess_page_bootstrap_ci=excess_ci,
                                  zero_denominator_bootstrap_draws=zero_draws,
                                  pages=[dict(page=page, sum=float(observed[k, i]), count=int(counts[k, i]),
                                              null_mean_sum=float(page_mean[k, i])) for i, page in enumerate(lay['pages'])])
        CHECKS['permutation_metric_statistics'] += REPS*3
        CHECKS['bootstrap_metric_resamples'] += BOOTS*3
        CHECKS['tested_splits'] += 1
        return dict(split=split, role='PRIMARY' if split == 'development' else 'SECTION_TRANSFER',
                    known_word_slots=lay['known'], unknown_word_slots=lay['unknown'], mobile_word_slots=lay['mobile'],
                    strata=len(lay['groups']), shuffle_seed=seeds['shuffle'], bootstrap_seed=seeds['bootstrap'],
                    permutations=REPS, bootstrap_resamples=BOOTS, metrics=metrics, **extras,
                    local_exact0_rate=extras['local_exact0_pairs']/extras['local_pairs'] if extras['local_pairs'] else None,
                    local_at_most1_rate=extras['local_at_most1_pairs']/extras['local_pairs'] if extras['local_pairs'] else None)


def audit_case(name):
    result_path = f'results/sequence-dependence/{name}.json'
    # This check occurs before even constructing a model. In particular, an
    # existing prepared native.json never authorizes a new native evaluation.
    require((ROOT/result_path).is_file(), ('Unscored case: audit refused', name))
    saved = read(result_path)
    path = f'data/sequence-dependence/{name}.json'
    compare(saved['source_sha256'], sha(path), name+'.source_sha256')
    compare(saved['input_freeze_sha256'], sha('data/sequence-dependence/input-freeze.json'))
    compare(saved['case'], name)
    compare(saved['status'], 'DESCRIPTIVE_DEPENDENCE_NOT_A_TRANSLATION')
    case = Case(read(path))
    compare(saved['training'], dict(zip(METRICS[:2], [m.summary for m in case.models])), name+'.training')
    for split in SPLITS[1:]:
        compare(saved['tests'][split], case.test(name, split), name+'.'+split)
    CHECKS['cases'] += 1
    print(f'PASS {name}', file=sys.stderr, flush=True)
    return saved, dict(case=name, status='PASS', result_sha256=sha(result_path),
                       independent_training_count_sha256={metric: model.count_hash for metric, model in zip(METRICS[:2], case.models)})


def audit_labels(results):
    saved = read('results/sequence-dependence/calibration.json')
    expected_hashes = {f'results/sequence-dependence/{name}.json': sha(f'results/sequence-dependence/{name}.json') for name in expected_cases()}
    compare(saved['result_hashes'], expected_hashes, 'calibration.result_hashes')
    labels = {}
    invariance = []
    for lang in LANGS:
        for mechanism in MECHANISMS:
            plain, cipher = (results[f'{lang}/{mechanism}-{channel}'] for channel in ('plain', 'bijective'))
            compare(plain['training'], cipher['training'], 'bijection.training')
            compare(plain['tests'], cipher['tests'], 'bijection.tests')
            invariance.append(dict(language=lang, mechanism=mechanism, status='PASS', paired_permutation_statistics=2*3*REPS))
    compare(saved['bijection_invariance'], invariance, 'calibration.bijection_invariance')
    for split in SPLITS[1:]:
        labels[split] = {}
        for channel in CHANNELS:
            fields = {}
            for metric in METRICS[:2]:
                positives = {lang: results[f'{lang}/ordered-{channel}']['tests'][split]['metrics'][metric]['positive'] for lang in LANGS}
                fields[metric] = dict(all_three_ordered_positive=all(positives.values()), by_language=positives)
            matrix = {lang: {mechanism: {key: results[f'{lang}/{mechanism}-{channel}']['tests'][split]['metrics']['local_exact1'][key]
                                       for key in ('informative', 'positive')} for mechanism in MECHANISMS} for lang in LANGS}
            passed = all(matrix[l][m]['informative'] for l in LANGS for m in MECHANISMS)
            passed = passed and all(matrix[l]['copy']['positive'] and not matrix[l]['ordered']['positive']
                                    and not matrix[l]['iid']['positive'] for l in LANGS)
            fields['copy_specific_local_exact1'] = dict(demonstrated_on_these_controls=bool(passed), by_language=matrix)
            labels[split][channel] = fields
    compare(saved['calibration_labels'], labels, 'calibration.labels')
    for field, expected in dict(status='COMPLETE_CONTROLS_BEFORE_NATIVE', encoded_datasets=36,
                                 independent_source_pairs=3, plaintext_invariance_datasets=9).items():
        compare(saved[field], expected, 'calibration.'+field)
    freeze = verify_manifest('results/sequence-dependence/calibration-freeze.json')
    require(set(expected_hashes) <= set(freeze['hashes']), 'Calibration freeze missing case results')
    compare(freeze['status'], 'FROZEN_BEFORE_NATIVE019_SCORES')
    CHECKS['bijection_invariance_cases'] += 9
    CHECKS['calibration_label_groups'] += 2*4*3


def dp(left, right):
    last = list(range(len(right)+1))
    for i, atom in enumerate(left, 1):
        now = [i]
        for j, other in enumerate(right, 1):
            now.append(min(now[-1]+1, last[j]+1, last[j-1]+(atom != other)))
        last = now
    return last[-1]


def self_test():
    words = sorted(w for size in range(1, 5) for w in itertools.product(('a', 'A', '@042;'), repeat=size))
    near = set(trie_neighbors(words).tolist())
    for a, left in enumerate(words):
        for b, right in enumerate(words):
            require((a*len(words)+b in near) == (dp(left, right) == 1), ('Trie DP', left, right))
    rows = []
    for split, page, values, resets in [('train', 't', ['aa', 'ab', 'aa', 'bb'], {0}),
                                        ('development', 'z', ['aa', 'ab', None, 'bb', 'ba'], {0}),
                                        ('development', 'z', ['aa', 'ba', 'ab'], set()),
                                        ('development', 'a', [None], {0})]:
        rows.append(dict(split=split, page=page, locus='toy', paragraph_first_line=False,
                         tokens=[dict(units=list(w) if w is not None else None, reset_before=i in resets)
                                 for i, w in enumerate(values)]))
    case = Case(dict(rows=rows))
    compare(case.layouts['development']['segments'], [[0, 1], [3, 4, 5, 6, 7]])
    compare(case.models[0].summary['training_pairs'], 3)
    observed = case.test('native', 'development')
    require(observed['metrics']['word_gain']['zero_denominator_bootstrap_draws'] > 0, 'Toy zero denominator missing')
    require(observed['metrics']['word_gain']['raw_page_bootstrap_ci'] is None, 'Toy invalid CI')
    empty = case.test('native', 'reserve')
    for metric in METRICS:
        require(empty['metrics'][metric]['observed'] is None and not empty['metrics'][metric]['positive'], 'Toy empty split')
    return dict(status='PASS', scope='SYNTHETIC_ONLY', trie_pairs_verified=len(words)**2,
                checks=dict(CHECKS), maximum_absolute_error=MAX_ERROR)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['self-test', 'calibration', 'native', 'all', 'case'])
    parser.add_argument('--case', help='Already-scored control case; native requires the native/all stage')
    args = parser.parse_args()
    started = time.monotonic()
    if args.stage == 'self-test':
        print(json.dumps(self_test(), indent=2)); return
    integrity = verify_inputs()
    if args.stage == 'case':
        require(args.case in expected_cases(), 'case must name an existing control')
        names = [args.case]
    else:
        names = expected_cases() if args.stage in ('calibration', 'all') else []
        if args.stage in ('native', 'all'):
            verify_manifest('results/sequence-dependence/calibration-freeze.json')
            verify_manifest('results/sequence-dependence/native-freeze.json')
            names.append('native')
    results, cases = {}, []
    for name in names:
        results[name], summary = audit_case(name)
        cases.append(summary)
    if args.stage in ('calibration', 'all'):
        audit_labels(results)
    if 'native' in results:
        compare(results['native']['calibration_sha256'], sha('results/sequence-dependence/calibration.json'))
    print(json.dumps(dict(status='PASS', stage=args.stage, integrity=integrity, checks=dict(CHECKS),
                          maximum_absolute_error=MAX_ERROR, cases=cases,
                          elapsed_seconds=time.monotonic()-started,
                          scope='Independent statistics; source/generator/encoding reconstruction audited separately'), indent=2))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(json.dumps(dict(status='FAIL', error_type=type(error).__name__, error=str(error),
                              checks=dict(CHECKS), maximum_absolute_error=MAX_ERROR), indent=2))
        raise
