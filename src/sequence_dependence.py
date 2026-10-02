"""Experiment019: fixed held-out sequence statistics and conditional permutations.

No cipher values are inferred. Calibration truth is never a scoring feature.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import hashlib
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/sequence-dependence'
OUT = ROOT / 'results/sequence-dependence'
LANGS = ['la', 'de', 'it']
MECHANISMS = ['ordered', 'iid', 'copy']
CHANNELS = ['bijective', 'uniform127', 'position127', 'hidden-space127']
SPLITS = ['train', 'development', 'reserve']
METRICS = ['word_gain', 'edge_gain', 'local_exact1']
REPS, BOOTS, TOL = 999, 1999, 1e-12


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def dump(path, obj):
    path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')


def verify(path):
    obj = read(path)
    assert obj['hashes'], path
    for filename, expected in obj['hashes'].items():
        assert sha(filename) == expected, (path, filename)
    return obj


def verify_prepared():
    original = verify('data/sequence-dependence/input-freeze.json')
    prepared = verify('data/sequence-dependence/prepared-manifest.json')
    assert prepared['input_hashes'] == original['hashes']
    assert prepared['input_freeze_sha256'] == sha('data/sequence-dependence/input-freeze.json')
    expected = {f'data/sequence-dependence/{name}.json' for name in ['template', 'native', 'config']}
    expected |= {f'data/sequence-dependence/{l}/{m}-{c}.json'
                 for l in LANGS for m in MECHANISMS for c in CHANNELS + ['plain']}
    assert expected <= set(prepared['hashes']), expected - set(prepared['hashes'])
    return read('data/sequence-dependence/config.json')


def lookup(sorted_keys, values, requested):
    if not len(sorted_keys):
        return np.zeros(len(requested), dtype=np.int64)
    ix = np.searchsorted(sorted_keys, requested)
    safe = np.minimum(ix, len(sorted_keys)-1)
    return np.where((ix < len(sorted_keys)) & (sorted_keys[safe] == requested), values[safe], 0)


class Corpus:
    """Only the rows/tokens observable interface is admitted here."""
    def __init__(self, data):
        self.words, word_ids = [], {}
        self.rows = {}
        for split in SPLITS:
            output = []
            for row in data['rows']:
                if row['split'] != split:
                    continue
                assert row['page'] not in {'f49v', 'f58r', 'f58v'}
                tokens = []
                for token in row['tokens']:
                    units = token['units']
                    if units is None:
                        wid = -1
                    else:
                        word = tuple(units)
                        assert word and all(isinstance(x, str) and x for x in word)
                        if word not in word_ids:
                            word_ids[word] = len(self.words)
                            self.words.append(word)
                        wid = word_ids[word]
                    tokens.append((wid, bool(token['reset_before'])))
                output.append(dict(page=row['page'], locus=row['locus'],
                                   paragraph_first_line=bool(row['paragraph_first_line']), tokens=tokens))
            self.rows[split] = output
        self.lengths = np.array([len(w) for w in self.words], dtype=np.int64)
        self.suffix, self.prefix = [], []
        suffix_ids, prefix_ids = {}, {}
        for word in self.words:
            for mapping, output, part in [(suffix_ids, self.suffix, word[-2:]),
                                          (prefix_ids, self.prefix, word[:2])]:
                if part not in mapping:
                    mapping[part] = len(mapping)
                output.append(mapping[part])
        self.suffix = np.array(self.suffix, dtype=np.int64)
        self.prefix = np.array(self.prefix, dtype=np.int64)
        self.nsuf, self.npre = len(suffix_ids), len(prefix_ids)
        self.splits = {s: self.layout(self.rows[s]) for s in SPLITS}
        self.models = [Predictor(self, False), Predictor(self, True)]
        self.local_lookup, self.near = self.similarity_graph()

    def layout(self, rows):
        words, resets, pages, roles, para, rid = [], [], [], [], [], []
        page_names = list(dict.fromkeys(r['page'] for r in rows))
        page_index = {p: i for i, p in enumerate(page_names)}
        for rno, row in enumerate(rows):
            size = len(row['tokens'])
            for k, (wid, reset) in enumerate(row['tokens']):
                words.append(wid); resets.append(reset); pages.append(page_index[row['page']])
                roles.append('only' if size == 1 else 'first' if k == 0 else 'last' if k == size-1 else 'middle')
                para.append(row['paragraph_first_line']); rid.append(rno)
        ids = np.array(words, dtype=np.int64)
        page = np.array(pages, dtype=np.int64)
        rowid = np.array(rid, dtype=np.int64)
        starts = np.zeros(len(ids), dtype=bool)
        ends = np.zeros(len(ids), dtype=bool)
        before, after, prior, current = [], [], [], []
        history = []
        for i, wid in enumerate(ids):
            if resets[i] or i == 0 or page[i] != page[i-1]:
                history = []
            if wid < 0:
                if i and ids[i-1] >= 0:
                    ends[i-1] = True
                history = []
                continue
            starts[i] = not history
            if starts[i] and i and ids[i-1] >= 0:
                ends[i-1] = True
            for j in history[-8:]:
                prior.append(j); current.append(i)
            history.append(i)
            history = history[-8:]
            if i and ids[i-1] >= 0 and rowid[i-1] == rowid[i] and not resets[i]:
                before.append(i-1); after.append(i)
        if len(ids) and ids[-1] >= 0:
            ends[-1] = True
        before = np.array(before, dtype=np.int64); after = np.array(after, dtype=np.int64)
        edge = (self.lengths[ids[before]] >= 2) & (self.lengths[ids[after]] >= 2)
        strata = defaultdict(list)
        for i, wid in enumerate(ids):
            if wid >= 0:
                key = (page_names[page[i]], para[i], roles[i], int(self.lengths[wid]), bool(starts[i]), bool(ends[i]))
                strata[key].append(i)
        groups = [np.array(strata[k], dtype=np.int64) for k in sorted(strata)]
        mobile = sum(len(g) for g in groups if len(set(ids[g].tolist())) > 1)
        return dict(ids=ids, page=page, page_names=page_names, rowid=rowid,
                    starts=starts, ends=ends, groups=groups, mobile=mobile,
                    known=int(np.count_nonzero(ids >= 0)), unknown=int(np.count_nonzero(ids < 0)),
                    pair_indices=[(before, after), (before[edge], after[edge])],
                    local_indices=(np.array(prior, dtype=np.int64), np.array(current, dtype=np.int64)))

    def similarity_graph(self):
        visible = sorted({int(x) for s in ['development', 'reserve'] for x in self.splits[s]['ids'] if x >= 0})
        mapper = np.full(len(self.words), -1, dtype=np.int64)
        mapper[visible] = np.arange(len(visible))
        words = [self.words[i] for i in visible]
        index = {w: i for i, w in enumerate(words)}
        near = np.zeros((len(words), len(words)), dtype=bool)
        signatures = defaultdict(list)
        for i, word in enumerate(words):
            for p in range(len(word)):
                deleted = word[:p] + word[p+1:]
                signatures[(len(word), p, deleted)].append(i)
                if deleted in index:
                    j = index[deleted]
                    near[i, j] = near[j, i] = True
        for group in signatures.values():
            if len(group) > 1:
                near[np.ix_(group, group)] = True
        np.fill_diagonal(near, False)
        return mapper, near

    def evaluate(self, split, ids=None):
        layout = self.splits[split]
        ids = layout['ids'] if ids is None else ids
        npages = len(layout['page_names'])
        sums = np.zeros((3, npages)); counts = np.zeros((3, npages), dtype=np.int64)
        coverage = []
        for k, model in enumerate(self.models):
            a, b = layout['pair_indices'][k]
            scores, cover = model.score(ids[a], ids[b])
            counts[k] = np.bincount(layout['page'][b], minlength=npages)
            sums[k] = np.bincount(layout['page'][b], weights=scores, minlength=npages)
            coverage.append(cover)
        a, b = layout['local_indices']
        count0 = np.count_nonzero(ids[a] == ids[b])
        aa, bb = self.local_lookup[ids[a]], self.local_lookup[ids[b]]
        assert np.all(aa >= 0) and np.all(bb >= 0)
        near = self.near[aa, bb]
        counts[2] = np.bincount(layout['page'][b], minlength=npages)
        sums[2] = np.bincount(layout['page'][b], weights=near, minlength=npages)
        extras = dict(prediction_coverage=dict(zip(METRICS[:2], coverage)),
                      local_pairs=len(a), local_exact0_pairs=int(count0),
                      local_exact1_pairs=int(np.count_nonzero(near)),
                      local_at_most1_pairs=int(count0 + np.count_nonzero(near)))
        return sums, counts, extras


class Predictor:
    def __init__(self, corpus, edge):
        self.edge = edge; self.corpus = corpus; self.tau = 10 if edge else 100
        k = int(edge); train = corpus.splits['train']
        left, right = train['pair_indices'][k]
        x, y = self.categories(train['ids'][left], train['ids'][right])
        nc = corpus.nsuf if edge else len(corpus.words)
        nt = corpus.npre if edge else len(corpus.words)
        self.nt = nt
        self.cc = np.bincount(x, minlength=nc)
        self.tc = np.bincount(y, minlength=nt)
        self.keys, self.values = np.unique(x * nt + y, return_counts=True)
        self.N = len(y); self.V = int(np.count_nonzero(self.tc))
        self.p0 = (self.tc + .5) / (self.N + .5*(self.V+1))

    def categories(self, left, right):
        return (self.corpus.suffix[left], self.corpus.prefix[right]) if self.edge else (left, right)

    def score(self, left, right):
        x, y = self.categories(left, right)
        pairs = lookup(self.keys, self.values, x*self.nt+y)
        gains = np.log((pairs/self.p0[y] + self.tau)/(self.cc[x] + self.tau))
        cover = dict(pairs=len(y), target_seen=int(np.count_nonzero(self.tc[y])),
                     context_seen=int(np.count_nonzero(self.cc[x])),
                     both_seen=int(np.count_nonzero((self.tc[y]>0) & (self.cc[x]>0))),
                     pair_seen=int(np.count_nonzero(pairs)))
        return gains, cover

    def summary(self):
        return dict(training_pairs=self.N, training_target_types=self.V,
                    training_context_types=int(np.count_nonzero(self.cc)),
                    training_pair_types=len(self.keys), unigram_pseudocount=.5,
                    backoff_mass=self.tau, unseen_target='one UNK bucket',
                    unseen_context='exact P0 backoff')


def fraction(a, b):
    return float(a/b) if b else None


def sample_test(corpus, split, seeds):
    layout = corpus.splits[split]
    observed, counts, extra = corpus.evaluate(split)
    totals = counts.sum(axis=1)
    scores = np.divide(observed.sum(axis=1), totals, out=np.full(3, np.nan), where=totals > 0)
    nulls = np.zeros((REPS, 3)); page_null_sum = np.zeros_like(observed)
    rng = np.random.default_rng(seeds['shuffle'])
    for b in range(REPS):
        ids = layout['ids'].copy()
        for group in layout['groups']:
            ids[group] = rng.permutation(layout['ids'][group])
        values, checked_counts, _ = corpus.evaluate(split, ids)
        assert np.array_equal(checked_counts, counts)
        page_null_sum += values
        nulls[b] = np.divide(values.sum(axis=1), totals, out=np.full(3, np.nan), where=totals > 0)
    page_null_mean = page_null_sum / REPS
    bootstrap_rng = np.random.default_rng(seeds['bootstrap'])
    npages = len(layout['page_names'])
    samples = bootstrap_rng.integers(npages, size=(BOOTS, npages))
    output = {}
    for k, metric in enumerate(METRICS):
        N = int(totals[k]); observed_value = float(scores[k]) if N else None
        mobile_fraction = fraction(layout['mobile'], layout['known'])
        mean = float(nulls[:, k].mean()) if N else None
        sd = float(nulls[:, k].std(ddof=1)) if N else None
        informative = bool(N and layout['mobile'] and np.ptp(nulls[:, k]) > TOL)
        p = (1+int(np.count_nonzero(nulls[:, k] >= scores[k]-TOL)))/(REPS+1) if N else None
        excess = observed_value - mean if N else None
        passed = bool(informative and p <= .05/3 and (observed_value > TOL if k < 2 else excess > TOL))
        denominators = counts[k][samples].sum(axis=1)
        zero_draws = int(np.count_nonzero(denominators == 0))
        raw_ci = excess_ci = None
        if not zero_draws:
            raw_draws = observed[k][samples].sum(axis=1)/denominators
            excess_draws = (observed[k]-page_null_mean[k])[samples].sum(axis=1)/denominators
            raw_ci = np.quantile(raw_draws, [.025, .975], method='linear').tolist()
            excess_ci = np.quantile(excess_draws, [.025, .975], method='linear').tolist()
        output[metric] = dict(observed=observed_value, sum=float(observed[k].sum()), comparisons=N,
                              null_mean=mean, null_sd=sd, excess=excess, p_greater=p,
                              informative=informative, positive=passed, threshold=.05/3,
                              mobile_fraction=mobile_fraction,
                              null_statistics=nulls[:, k].tolist() if N else [None]*REPS,
                              raw_page_bootstrap_ci=raw_ci, excess_page_bootstrap_ci=excess_ci,
                              zero_denominator_bootstrap_draws=zero_draws,
                              pages=[dict(page=p, sum=float(observed[k, i]), count=int(counts[k, i]),
                                          null_mean_sum=float(page_null_mean[k, i]))
                                     for i, p in enumerate(layout['page_names'])])
    return dict(split=split, role='PRIMARY' if split == 'development' else 'SECTION_TRANSFER',
                known_word_slots=layout['known'], unknown_word_slots=layout['unknown'],
                mobile_word_slots=layout['mobile'], strata=len(layout['groups']),
                shuffle_seed=seeds['shuffle'], bootstrap_seed=seeds['bootstrap'],
                permutations=REPS, bootstrap_resamples=BOOTS,
                metrics=output, **extra,
                local_exact0_rate=fraction(extra['local_exact0_pairs'], extra['local_pairs']),
                local_at_most1_rate=fraction(extra['local_at_most1_pairs'], extra['local_pairs']))


def evaluate_case(name, config):
    if (OUT / f'{name}.json').exists():
        raise FileExistsError(f'Never overwrite an existing case result: {name}')
    path = f'data/sequence-dependence/{name}.json'
    data = read(path); corpus = Corpus(data)
    seeds = config['seeds']['native'] if name == 'native' else config['seeds']['cases'][name]
    result = dict(status='DESCRIPTIVE_DEPENDENCE_NOT_A_TRANSLATION', case=name,
                  source_sha256=sha(path), input_freeze_sha256=sha('data/sequence-dependence/input-freeze.json'),
                  training={m: predictor.summary() for m, predictor in zip(METRICS[:2], corpus.models)},
                  tests={})
    for split in ['development', 'reserve']:
        result['tests'][split] = sample_test(corpus, split, seeds[split])
        summary = {m: {k: result['tests'][split]['metrics'][m][k]
                        for k in ['observed', 'excess', 'p_greater', 'informative', 'positive']} for m in METRICS}
        print(json.dumps(dict(case=name, split=split, metrics=summary)), flush=True)
    dump(f'results/sequence-dependence/{name}.json', result)
    return result


def close(a, b):
    if isinstance(a, dict):
        assert set(a) == set(b), (set(a), set(b))
        for k in a: close(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): close(x, y)
    elif isinstance(a, float):
        assert isinstance(b, (float, int)) and abs(a-b) <= 1e-10, (a, b)
    else:
        assert a == b, (a, b)


def calibrate():
    config = verify_prepared()
    assert not (OUT/'calibration.json').exists(), 'Never overwrite completed calibration'
    results, invariance = {}, []
    for lang in LANGS:
        for mechanism in MECHANISMS:
            plain = evaluate_case(f'{lang}/{mechanism}-plain', config)
            for channel in CHANNELS:
                name = f'{lang}/{mechanism}-{channel}'
                result = evaluate_case(name, config); results[name] = result
                if channel == 'bijective':
                    close(plain['training'], result['training'])
                    close(plain['tests'], result['tests'])
                    invariance.append(dict(language=lang, mechanism=mechanism, status='PASS',
                                           paired_permutation_statistics=2*3*REPS))
    labels = {}
    for split in ['development', 'reserve']:
        labels[split] = {}
        for channel in CHANNELS:
            row = {}
            for metric in METRICS[:2]:
                checks = {lang: results[f'{lang}/ordered-{channel}']['tests'][split]['metrics'][metric]['positive'] for lang in LANGS}
                row[metric] = dict(all_three_ordered_positive=all(checks.values()), by_language=checks)
            matrix = {lang: {m: {k: results[f'{lang}/{m}-{channel}']['tests'][split]['metrics']['local_exact1'][k]
                                  for k in ['informative', 'positive']} for m in MECHANISMS} for lang in LANGS}
            specific = all(v[m]['informative'] for v in matrix.values() for m in MECHANISMS) and all(
                v['copy']['positive'] and not v['ordered']['positive'] and not v['iid']['positive'] for v in matrix.values())
            row['copy_specific_local_exact1'] = dict(demonstrated_on_these_controls=bool(specific), by_language=matrix)
            labels[split][channel] = row
    summary = dict(status='COMPLETE_CONTROLS_BEFORE_NATIVE', encoded_datasets=36,
                   independent_source_pairs=3, plaintext_invariance_datasets=9,
                   bijection_invariance=invariance, calibration_labels=labels,
                   interpretation='Conditional dependence only. Sensitivity/discrimination labels apply to these works/channels, not to all language or all generators.',
                   result_hashes={f'results/sequence-dependence/{l}/{m}-{c}.json':sha(f'results/sequence-dependence/{l}/{m}-{c}.json')
                                  for l in LANGS for m in MECHANISMS for c in CHANNELS+['plain']})
    dump('results/sequence-dependence/calibration.json', summary)
    paths = ['notes/experiment-019-design.md','src/sequence_inputs.py','src/sequence_dependence.py',
             'data/sequence-dependence/input-freeze.json','data/sequence-dependence/prepared-manifest.json',
             'results/sequence-dependence/calibration.json']
    paths += list(summary['result_hashes'])
    dump('results/sequence-dependence/calibration-freeze.json',
         dict(status='FROZEN_BEFORE_NATIVE019_SCORES', hashes={p:sha(p) for p in paths}))
    print('CALIBRATION COMPLETE', json.dumps(labels), flush=True)


def native():
    config = verify_prepared(); verify('results/sequence-dependence/calibration-freeze.json')
    assert not (OUT/'native.json').exists(), 'Never overwrite a scored native profile'
    result = evaluate_case('native', config)
    result['calibration_sha256'] = sha('results/sequence-dependence/calibration.json')
    result['limitation'] = 'All manuscript data are development. None of these statistics assigns meaning, language or copying direction. Failed control sensitivity cannot support an exclusion.'
    dump('results/sequence-dependence/native.json', result)
    print('NATIVE COMPLETE', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['calibrate', 'native'])
    args = parser.parse_args()
    globals()[args.stage]()
