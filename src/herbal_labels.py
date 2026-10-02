"""Test preregistered picture correspondences without guessing plant names."""
from collections import defaultdict
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import json

from corpus import read

ROOT = Path(__file__).resolve().parents[1]


def distance(a, b):
    previous = list(range(len(b)+1))
    for i, aa in enumerate(a, 1):
        current = [i]
        for j, bb in enumerate(b, 1):
            current.append(min(current[-1]+1, previous[j]+1, previous[j-1]+(aa != bb)))
        previous = current
    return previous[-1]


def nearest(word, body):
    scores = {w: Fraction(distance(word,w),max(len(word),len(w))) for w in body}
    minimum = min(scores.values())
    return minimum, sorted(w for w,s in scores.items() if s == minimum)


def run(filename, spacing, selected, pairs, inventory):
    loci = read(ROOT/'data/raw'/filename, spacing)
    body, positions, all_positions, body_positions = defaultdict(set), defaultdict(list), defaultdict(list), defaultdict(list)
    first, transcribed, metadata = {}, {}, {}
    for x in loci:
        key = f'{x.page}.{x.num}'
        transcribed[key] = x
        metadata[x.page] = x.meta
        if x.code[0] == '!':
            continue
        for n,w in enumerate(x.words,1):
            if not w:
                continue
            position = f'{key}:word{n}'
            all_positions[w].append(position)
            if x.code[1] == 'P':
                body[x.page].add(w)
                positions[x.page,w].append(position)
                body_positions[w].append(position)
        if x.code[1] == 'P' and x.page not in first:
            first[x.page] = x.words[0]
    words = [inventory[x['locus']]['word'] for x in selected]
    targets = [x['target'] for x in selected]
    scores, hits = [], []
    for word in words:
        scores.append([nearest(word,body[t])[0] for t in targets])
        hits.append([int(word in body[t]) for t in targets])
    n = len(selected)
    null_h, null_d, stratified_d = [], [], []
    groups = [(metadata[t].get('L'), metadata[t].get('H')) for t in targets]
    for order in permutations(range(n)):
        null_h.append(sum(hits[i][order[i]] for i in range(n)))
        null_d.append(sum((scores[i][order[i]] for i in range(n)),Fraction(0))/n)
        if all(groups[i] == groups[order[i]] for i in range(n)):
            stratified_d.append(null_d[-1])
    observed_h = sum(hits[i][i] for i in range(n))
    observed_d = sum((scores[i][i] for i in range(n)),Fraction(0))/n
    result = dict(source=filename, spacing=spacing, n_pairs=n,
                  observed_exact_present=observed_h,
                  exact_null_expected_present=sum(null_h)/len(null_h),
                  p_greater_exact=sum(h>=observed_h for h in null_h)/len(null_h),
                  observed_mean_nearest_normalized_distance=float(observed_d),
                  null_mean_nearest_normalized_distance=float(sum(null_d)/len(null_d)),
                  p_lower_distance=sum(d<=observed_d for d in null_d)/len(null_d),
                  permutations=len(null_h), pairs=[], excluded_reading_diagnostics=[])
    result['posthoc_currier_scribe_sensitivity'] = dict(
        target_groups=groups, permutations=len(stratified_d),
        null_mean_distance=float(sum(stratified_d)/len(stratified_d)),
        p_lower_distance=sum(d<=observed_d for d in stratified_d)/len(stratified_d))
    for pair,word in zip(selected,words):
        target=pair['target']
        d, closest = nearest(word,body[target])
        result['pairs'].append(dict(**pair, word=word,
            exact_target_positions=positions[target,word], first_token=first[target],
            is_first_token=word==first[target], nearest_distance=float(d),
            closest=[dict(word=w, positions=positions[target,w]) for w in closest],
            corpus_all_exact_positions=all_positions[word],
            corpus_body_exact_positions=body_positions[word]))
    for pair in pairs:
        if not pair['locus'] or inventory[pair['locus']]['agreed_single_word']:
            continue
        reading = transcribed[pair['locus']].words
        if len(reading)==1 and reading[0]:
            word=reading[0]
            d, closest=nearest(word,body[pair['target']])
            result['excluded_reading_diagnostics'].append(dict(**pair,word=word,
                exact_target_positions=positions[pair['target'],word],
                corpus_all_exact_positions=all_positions[word],
                nearest_distance=float(d),closest=closest))
    return result


if __name__ == '__main__':
    assert distance('kitten','sitting')==3
    assert distance('chokam','chokam')==0
    pairs=json.loads((ROOT/'data/herbal-image-pairs.json').read_text())
    inventory={x['locus']:x for x in json.loads((ROOT/'results/pharma-labels.json').read_text())}
    selected=[x for x in pairs if x['locus'] and inventory[x['locus']]['agreed_single_word']]
    assert len(pairs)==9 and len(selected)==4
    result=dict(selection=pairs, primary_pairs=selected,
                runs=[run(f,s,selected,pairs,inventory) for f,s in
                      [('ZL3b-n.txt','split'),('IT2a-n.txt','split'),('ZL3b-n.txt','join')]])
    (ROOT/'results/herbal-label-association.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in result['runs']:
        print(r['source'],r['spacing'],{k:r[k] for k in ['n_pairs','observed_exact_present','exact_null_expected_present','p_greater_exact','observed_mean_nearest_normalized_distance','null_mean_nearest_normalized_distance','p_lower_distance']})
        for x in r['pairs']:
            print(x['word'],'->',x['target'],'nearest',x['closest'],'all occurrences',x['corpus_all_exact_positions'])
        print('Excluded sensitivity:',r['excluded_reading_diagnostics'])
