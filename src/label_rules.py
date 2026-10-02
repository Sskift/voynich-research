"""Finite, fixed boundary rewrite search; external herbal text is not scored here."""
from collections import Counter, defaultdict
from pathlib import Path
import hashlib,json,re,string
import numpy as np
from corpus import read

ROOT=Path(__file__).resolve().parents[1]
SEED=20261005


def transform(word,rule):
    if rule['side']=='identity': return word
    a,b=rule['old'],rule['new']
    if rule['side']=='prefix':
        if not word.startswith(a): return word
        result=b+word[len(a):]
    else:
        if not word.endswith(a): return word
        result=(word[:-len(a)] if a else word)+b
    return result if len(result)>=3 else word


def rule_id(side,a,b): return f'{side}:{a or "∅"}->{b or "∅"}'


def candidates(words):
    unique=sorted(set(words))
    rules=[dict(id='identity',side='identity',old='',new='',changed_types=0)]
    for side in ['prefix','suffix']:
        for a in ['']+list(string.ascii_lowercase):
            outputs=set(string.ascii_lowercase)
            if a:
                outputs|={''}|{a+c for c in string.ascii_lowercase}|{c+a for c in string.ascii_lowercase}
            if side=='suffix' and a=='m': outputs.add('iin')
            outputs.discard(a)
            for b in sorted(outputs):
                rule=dict(id=rule_id(side,a,b),side=side,old=a,new=b)
                changed=sum(transform(w,rule)!=w for w in unique)
                if changed>=8:
                    rule['changed_types']=changed
                    rules.append(rule)
    # Predeclared tie-break: edit count, replacement length, lexical rule ID.
    return [rules[0]]+sorted(rules[1:],key=lambda r:(3 if r['new']=='iin' else 1,len(r['new']),r['id']))


def feature_values(w):
    return [int(w.startswith(c)) for c in string.ascii_lowercase]+[int(w.endswith(c)) for c in string.ascii_lowercase]+[
        int(w.startswith('qo')),int(any(c in w for c in 'pf')),int('ch' in w or 'sh' in w)]


def profiles(rows,loci):
    pools=defaultdict(list)
    for x in loci:
        if x.code[1]=='P' and x.code[0]!='!':
            for w in x.words:
                if w: pools[x.page,len(w)].append(w)
    result=[]
    names=['initial_'+c for c in string.ascii_lowercase]+['final_'+c for c in string.ascii_lowercase]+['initial_qo','contains_pf','contains_ch_or_sh']
    for mode in ['tokens','types']:
        a,b,unmatched=[],[],[]
        for x in rows:
            pool=pools[x['page'],len(x['word'])]
            if not pool:
                unmatched.append(x['locus']);continue
            if mode=='types': pool=sorted(set(pool))
            a.append(feature_values(x['word']))
            b.append(np.mean([feature_values(w) for w in pool],axis=0))
        av,bv=np.mean(a,axis=0),np.mean(b,axis=0)
        result.append(dict(control=mode,n_matched=len(a),unmatched=unmatched,
            features=[dict(feature=f,label_rate=float(x),matched_body_rate=float(y),difference=float(x-y)) for f,x,y in zip(names,av,bv)]))
    return result


def association(matrix,indices,pages,rows):
    """Mean observed minus null gain; allows multiple Currier/hand groups."""
    n=len(indices)
    score=np.zeros(matrix.shape[0])
    strata=defaultdict(list)
    for i in indices: strata[rows[i]['meta'].get('L'),rows[i]['meta'].get('H')].append(i)
    for ids in strata.values():
        ii=np.array(ids)
        observed=matrix[:,ii,pages[ii]].sum(axis=1)
        expected=matrix[:,ii[:,None],pages[ii][None,:]].mean(axis=2).sum(axis=1)
        score+=observed-expected
    return score/n


def fit():
    allrows=json.loads((ROOT/'results/pharma-labels.json').read_text())
    labels=[x for x in allrows if x['subtype']=='Lf' and x['agreed_single_word']]
    weak=json.loads((ROOT/'data/herbal-weak-image-pairs.json').read_text())
    strong=json.loads((ROOT/'data/herbal-image-pairs.json').read_text())
    bykey={x['locus']:x for x in allrows}
    excluded={x['word'] for x in weak if x['primary']}
    excluded|={bykey[x['locus']]['word'] for x in strong if x['locus'] and bykey[x['locus']]['agreed_single_word']}
    rows=[x for x in labels if x['word'] not in excluded]
    # Retain only pharmaceutical pages in the fitting body dictionary.
    pharma={x['page'] for x in allrows}
    loci=read(ROOT/'data/raw/ZL3b-n.txt')
    body=defaultdict(set)
    for x in loci:
        if x.page in pharma and x.code[1]=='P' and x.code[0]!='!':
            body[x.page].update(w for w in x.words if w)
    assert not ({x['target'] for x in weak}&set(body))
    pnames=sorted(body)
    pindex={p:i for i,p in enumerate(pnames)}
    pages=np.array([pindex[x['page']] for x in rows])
    words=[x['word'] for x in rows]
    rules=candidates(words)
    hits=np.array([[[transform(w,r) in body[p] for p in pnames] for w in words] for r in rules],dtype=np.int8)
    gain=hits-hits[0:1]
    indices=list(range(len(rows)))
    scores=association(gain,indices,pages,rows)
    winner=int(np.argmax(scores))
    # Stabilize ties to mathematical equality; floating accumulation is not a tie-break.
    winner=int(np.flatnonzero(scores>=scores.max()-1e-12)[0])
    groups=defaultdict(list)
    for i,x in enumerate(rows): groups[x['meta'].get('L'),x['meta'].get('H')].append(i)
    rng=np.random.default_rng(SEED)
    observed_gain=gain[:,np.arange(len(rows)),pages].mean(axis=1)
    expectations=observed_gain-scores
    maxima=[]
    for _ in range(1999):
        shuffled=pages.copy()
        for ids in groups.values(): shuffled[ids]=rng.permutation(shuffled[ids])
        random_scores=gain[:,np.arange(len(rows)),shuffled].mean(axis=1)-expectations
        maxima.append(float(random_scores.max()))
    order=sorted(range(len(rules)),key=lambda i:(-round(float(scores[i]),12),i))
    cv=[]
    leaves=sorted({re.match(r'f(\d+)',x['page'])[1] for x in rows})
    for leaf in leaves:
        test=[i for i,x in enumerate(rows) if re.match(r'f(\d+)',x['page'])[1]==leaf]
        train=[i for i in indices if i not in test]
        local=association(gain,train,pages,rows)
        eligible=[i for i,r in enumerate(rules) if i==0 or sum(transform(w,r)!=w for w in {words[j] for j in train})>=8]
        top=max(float(local[i]) for i in eligible)
        chosen=next(i for i in eligible if local[i]>=top-1e-12)
        # Test relative to all pharmaceutical pages, not only other panels of one leaf.
        observed=int(gain[chosen,test,pages[test]].sum())
        expected=float(gain[chosen,test][:,pages].mean(axis=1).sum())
        cv.append(dict(leaf=leaf,n_test=len(test),rule=rules[chosen]['id'],
                       training_centered_gain=float(local[chosen]),
                       observed_test_gain=observed,null_expected_test_gain=expected,
                       test_centered_gain_per_label=(observed-expected)/len(test)))
    result=dict(stage='FROZEN_BEFORE_EXTERNAL_TARGET_SCORING',seed=SEED,
        n_training_labels=len(rows),n_training_types=len(set(words)),
        training_loci=[x['locus'] for x in rows],excluded_word_types=sorted(excluded),
        fitting_body_pages=pnames,n_candidates=len(rules),selection_source='ZL3b-n.txt split',
        selected_rule=rules[winner],selected_centered_gain_per_label=float(scores[winner]),
        training_max_statistic_p=(1+sum(v>=scores[winner]-1e-12 for v in maxima))/2000,
        training_null_max_mean=float(np.mean(maxima)),
        top_rules=[dict(**rules[i],centered_gain_per_label=float(scores[i]),
                        uncentered_gain=int(gain[i,np.arange(len(rows)),pages].sum()),
                        own_page_hits=int(hits[i,np.arange(len(rows)),pages].sum())) for i in order[:25]],
        leave_leaf_out=cv,
        artifact_inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'data/herbal-weak-image-pairs.json',ROOT/'results/pharma-labels.json',ROOT/'data/raw/ZL3b-n.txt']})
    (ROOT/'results/label-rule-fit.json').write_text(json.dumps(result,indent=2)+'\n')
    # Profile computation never uses herbal target text for a comparison.
    pr=dict(runs=[])
    for source,spacing in [('ZL3b-n.txt','split'),('IT2a-n.txt','split'),('ZL3b-n.txt','join')]:
        pr['runs'].append(dict(source=source,spacing=spacing,profiles=profiles(labels,read(ROOT/'data/raw'/source,spacing))))
    (ROOT/'results/label-body-profiles.json').write_text(json.dumps(pr,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': fit()
