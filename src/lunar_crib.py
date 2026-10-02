"""Frozen whole-circle crib screening; no phonetic values are assigned."""
from pathlib import Path
from collections import defaultdict
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPS=9999
SEED=20261008
COMPOUNDS=('cth','ckh','cph','cfh','ch','sh')


def units(text,mode):
    if mode=='eva_letters':return tuple(text)
    result=[];i=0
    while i<len(text):
        match=next((c for c in COMPOUNDS if text.startswith(c,i)),text[i])
        result.append(match);i+=len(match)
    return tuple(result)


def distance(a,b):
    row=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        new=[i]
        for j,y in enumerate(b,1):new.append(min(new[-1]+1,row[j]+1,row[j-1]+(x!=y)))
        row=new
    return row[-1]


def similarities(words):
    n=len(words);matrix=np.ones((n,n))
    for i in range(n):
        for j in range(i):matrix[i,j]=matrix[j,i]=1-distance(words[i],words[j])/max(len(words[i]),len(words[j]))
    return matrix


def normalize_rows(a):
    centered=a-a.mean(axis=1,keepdims=True)
    norm=np.sqrt((centered**2).sum(axis=1,keepdims=True))
    assert np.all(norm>0)
    return centered/norm


def summarize(a):
    return dict(observed=float(a[0]),null_mean=float(a[1:].mean()),null_sd=float(a[1:].std()),
        p_greater=float((1+np.count_nonzero(a[1:]>=a[0]-1e-12))/(REPS+1)),
        null_interval_95=[float(x) for x in np.quantile(a[1:],[.025,.975])])


def repeat_classes(rows,spacing):
    groups=defaultdict(list)
    for r in rows:
        word=r['zl_text'] if spacing=='retained' else r['zl_text'].replace('.','')
        groups[word].append(r['locus'])
    return {w:keys for w,keys in groups.items() if len(keys)>1}


def pattern(word):
    lookup={};out=[]
    for c in word:
        if c not in lookup:lookup[c]=len(lookup)
        out.append(lookup[c])
    return tuple(out)


def execute():
    paths=['notes/experiment-008-design.md','data/lunar-crib-names.json',
        'results/lunar-crib-inventory.json','data/raw/ZL3b-n.txt','data/raw/IT2a-n.txt']
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    inv=json.loads((ROOT/paths[2]).read_text());crib=json.loads((ROOT/paths[1]).read_text())
    assert crib['status']=='FROZEN_BEFORE_CRIB_SCORING'
    rows=[r for r in inv['labels'] if r['accepted']];n=len(rows)
    assert n==21 and len(inv['labels'])==28
    slots=np.array([r['slot'] for r in rows]);words=[r['zl_text'].replace('.','') for r in rows]
    left,right=np.triu_indices(n,1)
    phases=[dict(direction=d,start=s) for d in [1,-1] for s in range(28)]
    indices=np.array([[(p['start']+p['direction']*int(i))%28 for i in slots] for p in phases])
    targets={};structures={};models=[]
    for name,values in crib['representations'].items():
        assert len(values)==len(set(values))==28
        m=similarities([tuple(w) for w in values]);target=m[indices[:,left],indices[:,right]]
        targets[name]=normalize_rows(target)
    for mode in ['eva_letters','bench_units']:
        strings=[units(w,mode) for w in words];structures[mode]=similarities(strings)
        for name,values in crib['representations'].items():
            exact=[]
            for k,p in enumerate(phases):
                aligned=[tuple(values[i]) for i in indices[k]]
                length_count=sum(len(a)==len(b) for a,b in zip(strings,aligned))
                pattern_count=sum(pattern(a)==pattern(b) for a,b in zip(strings,aligned))
                exact.append(dict(**p,length_compatible_labels=length_count,within_word_pattern_compatible_labels=pattern_count))
            models.append(dict(id=mode+':'+name,unit_mode=mode,crib=name,strict_substitution_phase_diagnostics=exact))
    repeats={s:repeat_classes(rows,s) for s in ['retained','removed']}
    invariants=dict(repeated_agreed_labels=repeats,all_crib_representations_have28_distinct_items=True,
        exact_position_independent_full_name_decoder_possible=False,
        maximal_distinct_name_matches_by_label_equality={s:n-sum(len(v)-1 for v in groups.values()) for s,groups in repeats.items()},
        qualification='Assumes each full observed label alone deterministically decodes to its full distinct mansion name. Repetition observation preceded the design; this is an exact conditional contradiction, not a p value. Does not reject contextual codes, ambiguous abbreviation, non-name descriptions or all lunar interpretations.')
    runs=[]
    for null_index,null in enumerate(['all_slots','same_initial']):
        rng=np.random.default_rng(SEED+null_index)
        orders=np.tile(np.arange(n),(REPS+1,1))
        groups=[list(range(n))] if null=='all_slots' else [[i for i,w in enumerate(words) if w[0]==c] for c in sorted({w[0] for w in words})]
        for rep in range(1,REPS+1):
            for g in groups:orders[rep,g]=rng.permutation(g)
        family=np.full(REPS+1,-np.inf);model_runs=[]
        for model in models:
            m=structures[model['unit_mode']];target=targets[model['crib']]
            maxima=np.empty(REPS+1);observed=None
            for start in range(0,REPS+1,1000):
                o=orders[start:start+1000];a=m[o[:,left],o[:,right]]
                scores=normalize_rows(a)@target.T
                maxima[start:start+len(o)]=scores.max(axis=1)
                if start==0:observed=scores[0]
            family=np.maximum(family,maxima)
            best=int(np.argmax(observed));same=np.flatnonzero(np.abs(observed-observed[best])<1e-12)
            model_runs.append(dict(id=model['id'],score=summarize(maxima),best_phase=phases[best],
                tied_best_phases=[phases[i] for i in same],all_observed_phase_correlations=[float(x) for x in observed],
                phase_alignment_note='A score-maximizing alignment only, not decoded plaintext.',
                best_alignment=[dict(locus=r['locus'],voynich=r['zl_text'],crib_item=int(indices[best,i]+1),crib_text=crib['representations'][model['crib']][indices[best,i]]) for i,r in enumerate(rows)]))
        runs.append(dict(null=null,randomizations=REPS,seed=SEED+null_index,
            n_accepted_labels=n,n_pairwise_similarities=len(left),independent_sample_size_note='The210 pair scores share21 labels; inference permutes labels, not pair scores.',
            groups_of_exchangeable_loci=[[rows[i]['locus'] for i in g] for g in groups],
            global_six_model_maximum=summarize(family),models=model_runs))
        print(null,runs[-1]['global_six_model_maximum'],flush=True)
        for r in model_runs:print(r['id'],r['best_phase'],r['score'],flush=True)
    result=dict(input_hashes=hashes,exact_invariants=invariants,models=models,runs=runs,
        no_holdout=True,no_decoding_key_assigned=True)
    (ROOT/'results/lunar-crib-screen.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')


if __name__=='__main__':execute()
