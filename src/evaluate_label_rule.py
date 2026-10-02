"""Score a previously frozen rule on the prospectively listed weak picture pairs."""
from collections import defaultdict
from pathlib import Path
import hashlib,json
import numpy as np
from corpus import read
from label_rules import transform

ROOT=Path(__file__).resolve().parents[1]


def evaluate(source,spacing,selected,rule):
    body=defaultdict(set);positions=defaultdict(list);meta={}
    targets={x['target'] for x in selected}
    for x in read(ROOT/'data/raw'/source,spacing):
        if x.page not in targets: continue
        meta[x.page]=x.meta
        if x.code[1]!='P' or x.code[0]=='!': continue
        for j,w in enumerate(x.words,1):
            if w:
                body[x.page].add(w);positions[x.page,w].append(f'{x.page}.{x.num}:word{j}')
    words=[x['word'] for x in selected]
    target_list=[x['target'] for x in selected]
    rewritten=[transform(w,rule) for w in words]
    old=np.array([[w in body[p] for p in target_list] for w in words],dtype=int)
    new=np.array([[w in body[p] for p in target_list] for w in rewritten],dtype=int)
    difference=new-old
    groups=defaultdict(list)
    for i,p in enumerate(target_list): groups[meta[p].get('L'),meta[p].get('H')].append(i)
    expected=sum(float(difference[np.ix_(ids,ids)].mean(axis=1).sum()) for ids in groups.values())
    original_expected=sum(float(old[np.ix_(ids,ids)].mean(axis=1).sum()) for ids in groups.values())
    observed=int(np.trace(difference))
    rng=np.random.default_rng(20261005)
    null=[];identity_null=[]
    for _ in range(9999):
        order=np.arange(len(words))
        for ids in groups.values(): order[ids]=rng.permutation(order[ids])
        null.append(int(difference[np.arange(len(words)),order].sum()))
        identity_null.append(int(old[np.arange(len(words)),order].sum()))
    result=dict(source=source,spacing=spacing,n_pairs=len(words),n_changed=sum(x!=y for x,y in zip(words,rewritten)),
        identity_hits=int(np.trace(old)),rewritten_hits=int(np.trace(new)),observed_gain=observed,
        exact_null_expected_gain=expected,centered_gain=observed-expected,
        p_greater_gain=(1+sum(x>=observed for x in null))/10000,
        null_gain_variance=float(np.var(null)),null_uninformative=len(set(null))==1,
        identity_exact_null_expected_hits=original_expected,
        identity_p_greater=(1+sum(x>=int(np.trace(old)) for x in identity_null))/10000,
        identity_null_uninformative=len(set(identity_null))==1,
        strata={str(k):len(v) for k,v in groups.items()},pairs=[])
    for pair,w,r in zip(selected,words,rewritten):
        result['pairs'].append(dict(fragment=pair['fragment'],page=pair['page'],target=pair['target'],
            source=pair['source'],word=w,rewritten=r,changed=w!=r,
            original_positions=positions[pair['target'],w],rewritten_positions=positions[pair['target'],r]))
    return result


if __name__=='__main__':
    fit_path=ROOT/'results/label-rule-fit.json'
    fit=json.loads(fit_path.read_text())
    assert fit['stage']=='FROZEN_BEFORE_EXTERNAL_TARGET_SCORING'
    pairs=json.loads((ROOT/'data/herbal-weak-image-pairs.json').read_text())
    selected=[x for x in pairs if x['primary']]
    assert len(pairs)==22 and len(selected)==12
    assert {x['word'] for x in selected}<=set(fit['excluded_word_types'])
    assert not ({x['target'] for x in selected}&set(fit['fitting_body_pages']))
    result=dict(frozen_fit_sha256=hashlib.sha256(fit_path.read_bytes()).hexdigest(),
        rule=fit['selected_rule'],n_proposals=len(pairs),n_eligible=len(selected),
        runs=[evaluate(f,s,selected,fit['selected_rule']) for f,s in
              [('ZL3b-n.txt','split'),('IT2a-n.txt','split'),('ZL3b-n.txt','join')]])
    (ROOT/'results/label-rule-external-evaluation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
