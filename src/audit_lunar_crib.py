"""Independent observed-score and input audit for the frozen lunar crib screen."""
from pathlib import Path
from collections import Counter,defaultdict
from functools import lru_cache
import hashlib,json,math,re
from corpus import read

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(None)
def similarity(a,b):
    table=[[0]*(len(b)+1) for _ in range(len(a)+1)]
    for i in range(len(a)+1):table[i][0]=i
    for j in range(len(b)+1):table[0][j]=j
    for i,x in enumerate(a,1):
        for j,y in enumerate(b,1):
            table[i][j]=min(table[i-1][j]+1,table[i][j-1]+1,table[i-1][j-1]+int(x!=y))
    return 1-table[-1][-1]/max(len(a),len(b))


def corr(a,b):
    am=math.fsum(a)/len(a);bm=math.fsum(b)/len(b)
    cov=math.fsum((x-am)*(y-bm) for x,y in zip(a,b))
    return cov/math.sqrt(math.fsum((x-am)**2 for x in a)*math.fsum((y-bm)**2 for y in b))


def signature(word):
    return [word.index(c) for c in word]


def audit():
    result=json.loads((ROOT/'results/lunar-crib-screen.json').read_text())
    inv=json.loads((ROOT/'results/lunar-crib-inventory.json').read_text())
    crib=json.loads((ROOT/'data/lunar-crib-names.json').read_text())['representations']
    for name,sha in result['input_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    z={x.num:x for x in read(ROOT/'data/raw/ZL3b-n.txt') if x.page=='f69v' and x.code=='@Ri'}
    it={x.num:x for x in read(ROOT/'data/raw/IT2a-n.txt') if x.page=='f69v' and x.code=='@Ri'}
    assert set(z)==set(it)==set(range(4,32))
    for r in inv['labels']:
        num=r['slot']+4;x=z[num];y=it[num]
        assert r['locus']==f'f69v.{num}'
        assert r['accepted']==(all(x.words) and x.words==y.words)
        if r['accepted']:assert r['zl_text']==r['it_text']=='.'.join(x.words)
    rows=[r for r in inv['labels'] if r['accepted']];assert len(rows)==21
    reports=[]
    for mode in ['eva_letters','bench_units']:
        text=[r['zl_text'].replace('.','') for r in rows]
        a=[tuple(s) if mode=='eva_letters' else tuple(re.findall(r'cth|ckh|cph|cfh|ch|sh|.',s)) for s in text]
        pairs=[(i,j) for i in range(21) for j in range(i+1,21)]
        va=[similarity(a[i],a[j]) for i,j in pairs]
        for name,values in crib.items():
            assert len(set(values))==28
            scores=[];counts=[]
            for d in [1,-1]:
                for start in range(28):
                    b=[tuple(values[(start+d*r['slot'])%28]) for r in rows]
                    vb=[similarity(b[i],b[j]) for i,j in pairs]
                    scores.append(corr(va,vb))
                    counts.append((sum(len(x)==len(y) for x,y in zip(a,b)),sum(signature(x)==signature(y) for x,y in zip(a,b))))
            ident=mode+':'+name
            diag=next(m for m in result['models'] if m['id']==ident)['strict_substitution_phase_diagnostics']
            assert counts==[(p['length_compatible_labels'],p['within_word_pattern_compatible_labels']) for p in diag]
            for run in result['runs']:
                m=next(m for m in run['models'] if m['id']==ident)
                assert max(abs(a-b) for a,b in zip(scores,m['all_observed_phase_correlations']))<1e-12
                assert abs(max(scores)-m['score']['observed'])<1e-12
                groups=run['groups_of_exchangeable_loci'];flat=[x for g in groups for x in g]
                assert sorted(flat)==sorted(r['locus'] for r in rows)
                if run['null']=='same_initial':
                    initials={r['locus']:r['zl_text'][0] for r in rows}
                    assert all(len({initials[k] for k in g})==1 for g in groups)
            reports.append(dict(model=ident,direct_max_correlation=max(scores),phases_checked=len(scores)))
    for run in result['runs']:
        assert abs(max(m['score']['observed'] for m in run['models'])-run['global_six_model_maximum']['observed'])<1e-12
        assert run['randomizations']==9999
    repeats={}
    for mode in ['retained','removed']:
        groups=defaultdict(list)
        for r in rows:groups[r['zl_text'] if mode=='retained' else r['zl_text'].replace('.','')].append(r['locus'])
        repeats[mode]={k:v for k,v in groups.items() if len(v)>1}
        assert len(groups)==result['exact_invariants']['maximal_distinct_name_matches_by_label_equality'][mode]
    assert repeats==result['exact_invariants']['repeated_agreed_labels']
    out=dict(status='PASS',raw_transcription_and_mask_checked=True,frozen_hashes_checked=True,
        repeated_word_contradiction_checked=True,independent_full_DP_and_scalar_Pearson=True,
        score_scope='All336 observed phase scores recomputed independently; random samples were not independently regenerated. Exchangeability groups and reported observed family maximum checked.',models=reports)
    (ROOT/'results/audit-005.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':audit()
