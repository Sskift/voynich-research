"""Predeclared ring-order tests with shared permutations and reserved confirmation."""
from collections import defaultdict,Counter
from itertools import combinations
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from corpus import read
from zodiac_inventory import known_text
from herbal_labels import distance

ROOT=Path(__file__).resolve().parents[1]
REPS=9999
SEED=20261007


def length_bin(text):
    n=len(text.replace('.',''))
    return 0 if n<=4 else 1 if n<=6 else 2 if n<=8 else 3


def summarize(values):
    obs=float(values[0]);null=values[1:]
    return dict(observed=obs,null_mean=float(np.mean(null)),null_sd=float(np.std(null,ddof=1)),
        p_greater=float((1+np.count_nonzero(null>=obs-1e-12))/(len(null)+1)),
        null_interval_95=[float(x) for x in np.quantile(null,[.025,.975])],
        null_uninformative=bool(np.ptp(null)<1e-12))


def holm(ps):
    order=sorted(range(len(ps)),key=lambda i:ps[i]);out=[0.]*len(ps);previous=0.
    for rank,i in enumerate(order):
        previous=max(previous,min(1.,ps[i]*(len(ps)-rank)));out[i]=previous
    return out


def permutations_for(groups,texts,lookup,mode,rng):
    orders={};mobility=[]
    for group in groups:
        base=np.array([lookup[texts[k]] if k and k in texts else 0 for k in group['slots']],dtype=np.int32)
        samples=np.tile(base,(REPS+1,1))
        strata=defaultdict(list)
        for col,k in enumerate(group['slots']):
            if base[col]:strata[length_bin(texts[k]) if mode=='ring_length_bin' else 0].append(col)
        movable=0
        for cols in strata.values():
            if len(cols)<2:continue
            movable+=len(cols)
            perm=np.argsort(rng.random((REPS,len(cols))),axis=1)
            samples[1:,cols]=base[np.array(cols)[perm]]
        orders[group['id']]=samples
        mobility.append(dict(group=group['id'],n_accepted=int(np.count_nonzero(base)),
            exchangeable_slots=movable,stratum_sizes=[len(v) for v in strata.values()]))
    return orders,mobility


def pair_comparison(a,b,orders,S,E,texts):
    aa,bb=orders[a['id']],orders[b['id']]
    n=aa.shape[1];assert bb.shape[1]==n
    descriptions=[dict(direction=d,shift=s) for d in [1,-1] for s in range(n)]
    phases=np.array([[(p['shift']+p['direction']*i)%n for i in range(n)] for p in descriptions])
    counts=((aa[0]>0)[None,:] & (bb[0,phases]>0)).sum(axis=1)
    keep=counts>=3
    if not np.any(keep):return None,None,None
    phases=phases[keep];counts=counts[keep];descriptions=[p for p,k in zip(descriptions,keep) if k]
    similarities=S[aa[:,None,:],bb[:,phases]].sum(axis=2)/counts[None,:]
    exact=E[aa[:,None,:],bb[:,phases]].sum(axis=2)/counts[None,:]
    maxima=np.max(similarities,axis=1);exact_max=np.max(exact,axis=1)
    best=int(np.argmax(similarities[0]));phase=phases[best]
    matched=[]
    for i,j in enumerate(phase):
        if aa[0,i] and bb[0,j]:
            left,right=a['slots'][i],b['slots'][int(j)]
            matched.append(dict(left_locus=left,right_locus=right,left=texts[left],right=texts[right],similarity=float(S[aa[0,i],bb[0,j]])))
    info=dict(left=a['id'],right=b['id'],slots=n,eligible_phases=len(phases),
        best_phase=descriptions[best],best_phase_support=int(counts[best]),
        tied_best_phases=int(np.count_nonzero(abs(similarities[0]-maxima[0])<1e-12)),
        edit_similarity=summarize(maxima),exact_similarity=summarize(exact_max),best_matches=matched)
    return info,maxima,exact_max


def run(stage,spacing,mode,inventory):
    groups=[g for g in inventory['groups'] if g['stage']==stage and g['kind']=='circle']
    selected={x['locus'] for x in inventory['labels'] if x['primary'] and x['stage']==stage}
    texts={f'{x.page}.{x.num}':known_text(x) for x in read(ROOT/'data/raw/ZL3b-n.txt',spacing) if f'{x.page}.{x.num}' in selected}
    assert all(texts.values())
    vocab=['']+sorted(set(texts.values()));lookup={w:i for i,w in enumerate(vocab)}
    S=np.zeros((len(vocab),len(vocab)),dtype=float)
    for i in range(1,len(vocab)):
        for j in range(i,len(vocab)):
            S[i,j]=S[j,i]=1-distance(vocab[i],vocab[j])/max(len(vocab[i]),len(vocab[j]))
    E=np.eye(len(vocab));E[0,0]=0
    rng=np.random.default_rng(SEED+(stage=='confirmation'))
    orders,mobility=permutations_for(groups,texts,lookup,mode,rng)
    pairs=[];pair_scores=[];exact_scores=[];excluded=[]
    for a,b in combinations(groups,2):
        if not a['layout'] or a['layout']!=b['layout'] or a['band']!=b['band'] or a['page']==b['page']:continue
        info,score,exact=pair_comparison(a,b,orders,S,E,texts)
        if info is None:excluded.append([a['id'],b['id']]);continue
        pairs.append(info);pair_scores.append(score);exact_scores.append(exact)
    total=np.zeros(REPS+1);edge_total=0;ring_details=[]
    for group in groups:
        order=orders[group['id']];base=order[0];n=len(base)
        if np.count_nonzero(base)<3:continue
        edges=[(i,(i+1)%n) for i in range(n) if base[i] and base[(i+1)%n]]
        if not edges:continue
        left,right=np.array(edges).T
        scores=S[order[:,left],order[:,right]].sum(axis=1)
        total+=scores;edge_total+=len(edges)
        ring_details.append(dict(group=group['id'],n_edges=len(edges),similarity=summarize(scores/len(edges))))
    assert pair_scores and edge_total
    result=dict(stage=stage,source='ZL3b-n.txt',spacing=spacing,null=mode,randomizations=REPS,
        n_primary_labels=len(selected),n_circular_groups=len(groups),n_ring_pairs=len(pairs),
        excluded_ring_pairs=excluded,n_neighbour_edges=edge_total,
        primary_order=summarize(np.mean(pair_scores,axis=0)),
        secondary_exact_order=summarize(np.mean(exact_scores,axis=0)),
        neighbour_similarity=summarize(total/edge_total),pairs=pairs,rings=ring_details,mobility=mobility)
    corrected=holm([result['primary_order']['p_greater'],result['neighbour_similarity']['p_greater']])
    result['primary_order']['holm_p']=corrected[0];result['neighbour_similarity']['holm_p']=corrected[1]
    return result


def execute(stage):
    inv_path=ROOT/'results/zodiac-inventory.json';inventory=json.loads(inv_path.read_text())
    snapshots={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [inv_path,ROOT/'data/zodiac-layout.json',ROOT/'notes/experiment-007-design.md']}
    if stage=='confirmation':
        discovery=ROOT/'results/zodiac-order-discovery.json'
        assert discovery.exists(),'Save discovery result before scoring confirmation diagrams'
        snapshots[str(discovery.relative_to(ROOT))]=hashlib.sha256(discovery.read_bytes()).hexdigest()
    runs=[]
    for spacing in ['split','join']:
        for null in ['ring','ring_length_bin']:
            r=run(stage,spacing,null,inventory);runs.append(r)
            print(stage,spacing,null,{k:r[k] for k in ['n_primary_labels','n_ring_pairs','n_neighbour_edges','primary_order','neighbour_similarity']},flush=True)
    # IT split has exactly the same retained label strings by construction.
    keys={x['locus'] for x in inventory['labels'] if x['primary'] and x['stage']==stage}
    arrays=[]
    for f in ['ZL3b-n.txt','IT2a-n.txt']:
        arrays.append({f'{x.page}.{x.num}':known_text(x) for x in read(ROOT/'data/raw'/f) if f'{x.page}.{x.num}' in keys})
    assert arrays[0]==arrays[1]
    out=dict(stage=stage,input_hashes=snapshots,
        IT_split_note='Exactly identical retained strings by the consensus inclusion rule; not independent replication and not re-simulated.',
        runs=runs)
    (ROOT/f'results/zodiac-order-{stage}.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['discovery','confirmation'])
    execute(parser.parse_args().stage)
