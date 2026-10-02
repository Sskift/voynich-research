"""Align pharmaceutical labels and test exact own-page body recurrence."""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import csv
import json

import numpy as np
from corpus import read

ROOT = Path(__file__).resolve().parents[1]
REPS = 9999
SEED = 20261003


def key(x):
    return f'{x.page}.{x.num}'


def labels(loci):
    return [x for x in loci if x.code[1:] in ('Lf','Lc') and x.code[0] != '!']


def images_for(page):
    # Mirror filenames are photograph identifiers, NOT IVTFF folio IDs.
    # f101v spans two photographs; retain both instead of dropping half a page.
    coverage = json.loads((ROOT/'data/pharma-image-coverage.json').read_text())
    return coverage[page]['images']


def build_inventory():
    zl = read(ROOT/'data/raw/ZL3b-n.txt')
    it = read(ROOT/'data/raw/IT2a-n.txt')
    ti = {key(x): x for x in it}
    rows = []
    for x in labels(zl):
        y = ti.get(key(x))
        reason = []
        if not y or y.code[0] == '!' or y.code[1:] != x.code[1:]:
            reason.append('missing_or_different_IT_locus_type')
        if len(x.words) != 1:
            reason.append('ZL_multiple_tokens')
        if any(w is None for w in x.words):
            reason.append('ZL_uncertain_reading')
        if y and (len(y.words) != 1 or any(w is None for w in y.words)):
            reason.append('IT_multiple_or_uncertain_tokens')
        if y and x.words != y.words:
            reason.append('readings_disagree')
        rows.append(dict(locus=key(x),page=x.page,num=x.num,subtype=x.code[1:],
                         zl_raw=x.raw,it_raw=y.raw if y else None,
                         zl_words=x.words,it_words=y.words if y else None,
                         agreed_single_word=not reason,
                         word=x.words[0] if not reason else None,
                         exclusion_reasons=reason,meta=x.meta,images=images_for(x.page)))
    (ROOT/'results/pharma-labels.json').write_text(json.dumps(rows,indent=2)+'\n')
    with (ROOT/'results/pharma-labels.tsv').open('w') as f:
        writer = csv.writer(f,delimiter='\t')
        writer.writerow(['locus','subtype','zl_raw','it_raw','agreed','word','exclusions','images'])
        for x in rows:
            writer.writerow([x['locus'],x['subtype'],x['zl_raw'],x['it_raw'],
                             x['agreed_single_word'],x['word'],','.join(x['exclusion_reasons']),';'.join(x['images'])])
    repeated = defaultdict(list)
    for x in rows:
        if x['agreed_single_word']:
            repeated[x['word']].append(x)
    groups = []
    for word, xs in sorted(repeated.items()):
        if len(xs) >= 2:
            groups.append(dict(word=word,occurrences=[
                {k:x[k] for k in ('locus','page','subtype','images')} for x in xs],
                cross_page_pairs=[[x['locus'],y['locus']] for x,y in combinations(xs,2)
                                  if x['page'] != y['page']]))
    (ROOT/'results/pharma-repeated-labels.json').write_text(json.dumps(groups,indent=2)+'\n')
    summary = dict(total=len(rows), by_subtype=dict(Counter(x['subtype'] for x in rows)),
                   agreed=sum(x['agreed_single_word'] for x in rows),
                   agreement_by_subtype=dict(Counter(x['subtype'] for x in rows if x['agreed_single_word'])),
                   exclusions=dict(Counter(r for x in rows for r in x['exclusion_reasons'])),
                   repeated_agreed_word_types=len(groups),
                   cross_page_repeated_agreed_types=sum(bool(x['cross_page_pairs']) for x in groups))
    return rows, summary


def permutation_test(rows, body, name, rng):
    selected = [x for x in rows if x['agreed_single_word'] and (name=='ALL' or x['subtype']==name)]
    strata = defaultdict(list)
    for x in selected:
        assert body[x['page']], x['page']
        strata[(x['subtype'],x['meta'].get('L','?'),x['meta'].get('H','?'))].append(x)
    null_hit, null_rate = np.zeros(REPS), np.zeros(REPS)
    observed_hit, observed_rate, expected_hit, expected_rate = 0, 0., 0., 0.
    matches = []
    for stratum, xs in sorted(strata.items()):
        counts = np.array([[body[x['page']][y['word']] for y in xs] for x in xs])
        hit = (counts>0).astype(int)
        lengths = np.array([sum(body[x['page']].values()) for x in xs])
        rate = 1000*counts/lengths[:,None]
        observed_hit += int(np.trace(hit))
        observed_rate += float(np.trace(rate))
        expected_hit += float(hit.mean(axis=1).sum())
        expected_rate += float(rate.mean(axis=1).sum())
        n = len(xs)
        for start in range(0,REPS,500):
            end = min(REPS,start+500)
            order = np.argsort(rng.random((end-start,n)),axis=1)
            null_hit[start:end] += hit[np.arange(n)[None,:],order].sum(axis=1)
            null_rate[start:end] += rate[np.arange(n)[None,:],order].sum(axis=1)
        for i,x in enumerate(xs):
            if counts[i,i]:
                matches.append(dict(locus=x['locus'],word=x['word'],body_count=int(counts[i,i]),
                                    body_tokens=int(lengths[i])))
    return dict(n_labels=len(selected), observed_present=observed_hit,
                exact_null_expected_present=expected_hit,
                null_mc_mean=float(null_hit.mean()), null_mc_sd=float(null_hit.std(ddof=1)),
                p_greater_presence=float((1+(null_hit>=observed_hit).sum())/(REPS+1)),
                null_presence_95_interval=[float(x) for x in np.quantile(null_hit,[.025,.975])],
                observed_mean_frequency_per_1000=observed_rate/len(selected),
                null_expected_mean_frequency_per_1000=expected_rate/len(selected),
                p_greater_frequency=float((1+(null_rate>=observed_rate-1e-10).sum())/(REPS+1)),
                matches=matches,
                strata={str(k):len(v) for k,v in strata.items()})


if __name__ == '__main__':
    rows, summary = build_inventory()
    all_results = dict(inventory=summary, randomizations=REPS, seed=SEED, runs=[])
    for filename,spacing in [('ZL3b-n.txt','split'),('IT2a-n.txt','split'),('ZL3b-n.txt','join')]:
        loci = read(ROOT/'data/raw'/filename,spacing)
        body = defaultdict(Counter)
        word_positions = defaultdict(list)
        for x in loci:
            if x.code[1]=='P' and x.code[0]!='!':
                body[x.page].update(w for w in x.words if w)
                for i,w in enumerate(x.words):
                    if w:
                        word_positions[(x.page,w)].append(f'{key(x)}:word{i+1}')
        rng = np.random.default_rng(SEED)
        run = dict(source=filename,spacing=spacing,tests={})
        for name in ['Lf','Lc','ALL']:
            result=permutation_test(rows,body,name,rng)
            for match in result['matches']:
                page=match['locus'].split('.')[0]
                match['body_positions']=word_positions[(page,match['word'])]
            run['tests'][name]=result
        all_results['runs'].append(run)
    (ROOT/'results/pharma-label-page-association.json').write_text(json.dumps(all_results,indent=2)+'\n')
    print(json.dumps(all_results,indent=2))
