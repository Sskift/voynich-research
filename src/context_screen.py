"""Exploratory check: can preceding-word context distinguish terminal-m readings?"""
import json
import math
from collections import Counter, defaultdict

from analyze import ROOT, heldout
from corpus import read


def run(name, spacing):
    lines = [x for x in read(ROOT/'data/raw'/name, spacing) if x.eligible]
    unigram, bigram, preceding = Counter(), Counter(), Counter()
    for x in lines:
        if heldout(x.page):
            continue
        for i, word in enumerate(x.words[:-1]):
            if word is None:
                continue
            unigram[word] += 1
            if i and x.words[i-1]:
                previous = x.words[i-1]
                bigram[previous, word] += 1
                preceding[previous] += 1
    denom = sum(unigram.values())+.5*(len(unigram)+1)
    controls = {}
    for ending in ['m', 'r', 'l', 'n', 'y']:
        test = [(x, x.words[-2], x.words[-1]) for x in lines
                if heldout(x.page) and len(x.words) >= 3 and x.words[-2]
                and x.words[-1] and x.words[-1].endswith(ending)]
        rows = []
        for scope in ['ALL', 'A', 'B']:
            scoped = [(x, p, w) for x, p, w in test if scope == 'ALL' or x.meta.get('L') == scope]
            if not scoped:
                continue
            for replacement in ['KEEP', 'r', 'l', 'n', 'y', 'in', 'iin']:
                for tau in [10, 50, 200]:
                    sums = [0., 0., 0.]
                    for x, prev, word in scoped:
                        target = word if replacement == 'KEEP' else word[:-1]+replacement
                        pu = (unigram[target]+.5)/denom
                        pc = (bigram[prev, target]+tau*pu)/(preceding[prev]+tau)
                        sums[0] += math.log2(pu)
                        sums[1] += math.log2(pc)
                        sums[2] += int((prev, target) in bigram)
                    u, c, seen = (v/len(scoped) for v in sums)
                    rows.append(dict(scope=scope, replacement=replacement, tau=tau,
                                     n=len(scoped), unigram_log2=u, conditional_log2=c,
                                     context_lift_bits=c-u, pair_seen_rate=seen))
        controls[ending] = rows
    types = Counter(w for x in lines if heldout(x.page) for w in x.words[-1:]
                    if w and w.endswith('m'))
    examples = []
    for word, n in types.most_common(20):
        examples.append(dict(word=word, count=n, candidates={s: {
            'form': word[:-1]+s, 'train_nonfinal_count': unigram[word[:-1]+s]}
            for s in ['r','l','iin']}))
    result = dict(source=name, spacing=spacing, controls=controls, examples=examples,
                  caveat='Exploratory after viewing experiment 001; not a semantic decipherment.')
    stem = name.removesuffix('.txt')+'-'+spacing
    (ROOT/'results'/f'{stem}-context.json').write_text(json.dumps(result, indent=2)+'\n')
    print(stem)
    for row in controls['m']:
        if row['scope'] == 'ALL' and row['tau'] == 50:
            print(row)


if __name__ == '__main__':
    for name, spacing in [('ZL3b-n.txt','split'), ('IT2a-n.txt','split'), ('ZL3b-n.txt','join')]:
        run(name, spacing)
