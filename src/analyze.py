"""Reproducible layout constraints and an explicitly nonsemantic suffix screen."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from corpus import read, checks

ROOT = Path(__file__).resolve().parents[1]
REPS = 1999
SEED = 20261002


def bootstrap(pages, sums, counts, rng):
    """Page-cluster bootstrap of the sample mean, retaining unequal cluster sizes."""
    keys = sorted(set(pages))
    groups = {p: i for i, p in enumerate(keys)}
    total, ns = np.zeros(len(keys)), np.zeros(len(keys))
    for p, s, n in zip(pages, sums, counts):
        total[groups[p]] += s
        ns[groups[p]] += n
    ix = rng.integers(0, len(keys), (REPS, len(keys)))
    vals = total[ix].sum(axis=1) / ns[ix].sum(axis=1)
    return [float(x) for x in np.quantile(vals, [.025, .975])]


def by_group(rows, name):
    grouped = defaultdict(list)
    for line, value in rows:
        grouped[line.meta.get(name, '?')].append(value)
    return {g: {'n': len(v), 'mean': float(np.mean(v))} for g, v in sorted(grouped.items())}


def paired_edge(lines, feature, edge, rng):
    rows, edge_values, inner_values, probabilities = [], [], [], []
    for x in lines:
        # Require complete known lines for exact exchangeability and fixed positions.
        if len(x.words) < 3 or any(w is None for w in x.words):
            continue
        values = np.array([int(feature(w)) for w in x.words], dtype=float)
        index = 0 if edge == 'first' else len(values)-1
        observed = values[index]
        inner = (values.sum() - observed) / (len(values)-1)
        rows.append((x, observed-inner))
        edge_values.append(observed)
        inner_values.append(inner)
        probabilities.append((values.mean(), values.sum(), len(values)))
    observed = float(np.mean([v for x, v in rows]))
    null = np.zeros(REPS)
    for prob, count, n in probabilities:
        chosen = rng.binomial(1, prob, REPS)
        null += chosen - (count-chosen)/(n-1)
    null /= len(rows)
    return {'n_lines': len(rows), 'edge_rate': float(np.mean(edge_values)),
            'nonedge_line_mean_rate': float(np.mean(inner_values)),
            'paired_difference': observed,
            'page_bootstrap_95ci': bootstrap([x.page for x, v in rows],
                [v for x, v in rows], [1]*len(rows), rng),
            'permutation_p_greater': float((1+np.count_nonzero(null >= observed))/(REPS+1)),
            'null_sd': float(null.std(ddof=1)),
            'by_currier': by_group(rows, 'L'), 'by_scribe': by_group(rows, 'H')}


def paragraph_test(lines, rng):
    pages = defaultdict(list)
    counts = {'first': [0, 0], 'other': [0, 0]}
    initial = {'first': [0, 0], 'other': [0, 0]}
    for x in lines:
        if x.paragraph is None:
            continue
        words = [w for w in x.words if w is not None]
        if not words:
            continue
        key = 'first' if x.paragraph_start else 'other'
        hits = sum('p' in w or 'f' in w for w in words)
        counts[key][0] += hits
        counts[key][1] += len(words)
        if x.words[0] is not None:
            initial[key][0] += x.words[0][0] in 'ktpf'
            initial[key][1] += 1
        pages[x.page].append((x, hits/len(words)))
    rows, null = [], np.zeros(REPS)
    for p, data in pages.items():
        values = np.array([v for x, v in data])
        flags = np.array([x.paragraph_start for x, v in data])
        n1 = int(flags.sum())
        if n1 == 0 or n1 == len(data):
            continue
        diff = values[flags].mean() - values[~flags].mean()
        rows.append((data[0][0], float(diff)))
        # Random first-line designations within a page; fixed count.
        ix = np.argsort(rng.random((REPS, len(data))), axis=1)[:, :n1]
        selected = values[ix].sum(axis=1)
        null += selected/n1 - (values.sum()-selected)/(len(data)-n1)
    null /= len(rows)
    observed = float(np.mean([v for x, v in rows]))
    def ratios(d):
        return {k: {'hits': a, 'n': b, 'rate': a/b} for k, (a, b) in d.items()}
    return {'pf_token_rates': ratios(counts), 'gallows_line_start_rates': ratios(initial),
            'n_paired_pages': len(rows), 'mean_page_difference_pf': observed,
            'page_bootstrap_95ci': bootstrap([x.page for x, v in rows],
                [v for x, v in rows], [1]*len(rows), rng),
            'permutation_p_greater': float((1+np.count_nonzero(null >= observed))/(REPS+1))}


def repetition(lines, rng, interior=False):
    groups = defaultdict(list)
    bypage = defaultdict(lambda: [0, 0., 0])
    for x in lines:
        if any(w is None for w in x.words):
            continue
        words = x.words[1:-1] if interior else x.words
        if len(words) < 3:
            continue
        # Vocabulary IDs need only be consistent within this line.
        ids = {w: i for i, w in enumerate(set(words))}
        values = [ids[w] for w in words]
        observed = sum(a == b for a, b in zip(words, words[1:]))
        c = Counter(words)
        expected = sum(v*(v-1) for v in c.values())/len(words)
        bypage[x.page][0] += observed
        bypage[x.page][1] += expected
        bypage[x.page][2] += len(words)-1
        groups[len(values)].append(values)
    null = np.zeros(REPS, dtype=np.int64)
    for n, group in sorted(groups.items()):
        arr = np.array(group)
        # Chunk over simulations to keep memory bounded.
        for start in range(0, REPS, 25):
            size = min(25, REPS-start)
            order = np.argsort(rng.random((size, len(arr), n)), axis=2)
            shuffled = np.take_along_axis(np.broadcast_to(arr, order.shape), order, axis=2)
            null[start:start+size] += (shuffled[:, :, 1:] == shuffled[:, :, :-1]).sum(axis=(1, 2))
    observed = sum(v[0] for v in bypage.values())
    expected = sum(v[1] for v in bypage.values())
    pairs = sum(v[2] for v in bypage.values())
    return {'n_lines': sum(len(v) for v in groups.values()), 'n_pairs': pairs,
            'observed_equal_pairs': observed, 'expected_equal_pairs': expected,
            'observed_rate': observed/pairs, 'null_expected_rate': expected/pairs,
            'ratio': observed/expected, 'null_mc_mean': float(null.mean()),
            'null_mc_sd': float(null.std(ddof=1)),
            'permutation_p_greater': float((1+np.count_nonzero(null >= observed))/(REPS+1)),
            'page_bootstrap_excess_rate_95ci': bootstrap(list(bypage),
                [v[0]-v[1] for v in bypage.values()], [v[2] for v in bypage.values()], rng)}


def heldout(page):
    # Both sides and every panel of the same leaf always use the same split.
    leaf = re.match(r'f(\d+)', page)
    key = leaf[1] if leaf else page
    return int(hashlib.sha256(key.encode()).hexdigest(), 16) % 3 == 0


def suffix_screen(lines):
    train = Counter(w for x in lines if not heldout(x.page) for w in x.words[:-1] if w)
    test = [x.words[-1] for x in lines if heldout(x.page) and len(x.words) >= 3
            and x.words[-1] and x.words[-1].endswith('m')]
    # Fixed vocabulary includes an unknown bucket; all methods have same denominator.
    denominator = sum(train.values()) + .5*(len(train)+1)
    result = []
    for suffix in ['KEEP', '', 'r', 'l', 'n', 'y', 'in', 'iin', 'iiin']:
        converted = [w if suffix == 'KEEP' else w[:-1]+suffix for w in test]
        hits = sum(w in train for w in converted)
        result.append({'replacement': suffix, 'n_test': len(test), 'in_train': hits,
                       'coverage': hits/len(test),
                       'mean_log2_score': sum(math.log2((train[w]+.5)/denominator)
                                              for w in converted)/len(test)})
    return {'training_nonfinal_tokens': sum(train.values()),
            'training_types': len(train), 'heldout_final_m_tokens': len(test),
            'candidates': result,
            'test_type_counts': dict(Counter(test).most_common()),
            'warning': 'Vocabulary compatibility only; neither allography nor meaning established.'}


def run(name, spacing):
    rng = np.random.default_rng(SEED)
    loci = read(ROOT/'data/raw'/name, spacing)
    lines = [x for x in loci if x.eligible]
    known = [w for x in loci for w in x.words if w]
    q_followers = Counter()
    for w in known:
        for i, ch in enumerate(w):
            if ch == 'q':
                q_followers[w[i+1] if i+1 < len(w) else '<END>'] += 1
    result = {'source': name, 'uncertain_spaces': spacing, 'seed': SEED, 'randomizations': REPS,
              'audit': {'loci': len(loci), 'pages': len(set(x.page for x in loci)),
                        'paragraphs': sum(x.paragraph_start for x in loci),
                        'raw_token_slots': sum(len(x.words) for x in loci),
                        'accepted_tokens': len(known), 'rejected_tokens': sum(w is None for x in loci for w in x.words),
                        'eligible_P0_lines': len(lines)},
              'q_followers': dict(q_followers.most_common()),
              'paragraph': paragraph_test(lines, rng),
              'line_final_m': paired_edge(lines, lambda w: w.endswith('m'), 'last', rng),
              'line_initial_gallows': paired_edge(lines, lambda w: w[0] in 'ktpf', 'first', rng),
              'repetition_full_line': repetition(lines, rng),
              'repetition_interior': repetition(lines, rng, True),
              'suffix_screen': suffix_screen(lines)}
    stem = name.removesuffix('.txt')+'-'+spacing
    (ROOT/'results'/f'{stem}.json').write_text(json.dumps(result, indent=2)+'\n')
    with (ROOT/'results'/f'{stem}-tokens.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['page','locus','code','currier','scribe','section','paragraph',
                         'paragraph_start','position','slots','raw_token','accepted_token','heldout'])
        for x in loci:
            for i, (raw, word) in enumerate(zip(x.original_words, x.words)):
                writer.writerow([x.page,x.num,x.code,x.meta.get('L',''),x.meta.get('H',''),
                    x.meta.get('I',''),x.paragraph,x.paragraph_start,i,len(x.words),raw,word,heldout(x.page)])
    print(json.dumps({'run':stem,'audit':result['audit'],
        'm_edge':result['line_final_m'], 'pf':result['paragraph'],
        'repetition_interior':result['repetition_interior'],
        'suffix_candidates':result['suffix_screen']['candidates']},ensure_ascii=False),flush=True)
    return result


if __name__ == '__main__':
    checks()
    for name, spacing in [('ZL3b-n.txt','split'), ('IT2a-n.txt','split'), ('ZL3b-n.txt','join')]:
        run(name, spacing)
