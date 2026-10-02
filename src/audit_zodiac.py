"""Direct-loop audit of order scores and analytic neighbour null expectations."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib,json,math
from corpus import read
from zodiac_inventory import known_text

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(None)
def sim(a,b):
    # Full 2-D dynamic programming independently of the experiment's distance helper.
    d=[[0]*(len(b)+1) for _ in range(len(a)+1)]
    for i in range(len(a)+1):d[i][0]=i
    for j in range(len(b)+1):d[0][j]=j
    for i in range(1,len(a)+1):
        for j in range(1,len(b)+1):
            d[i][j]=min(d[i-1][j]+1,d[i][j-1]+1,d[i-1][j-1]+(a[i-1]!=b[j-1]))
    return 1-d[-1][-1]/max(len(a),len(b))


def bin_of(s):
    n=len(s.replace('.',''))
    return 0 if n<=4 else 1 if n<=6 else 2 if n<=8 else 3


def audit():
    inv=json.loads((ROOT/'results/zodiac-inventory.json').read_text())
    assert len(inv['labels'])==299
    assert sum(x['primary'] for x in inv['labels'])==186
    assert len({x['locus'] for x in inv['labels']})==299
    locations=[k for g in inv['groups'] for k in g['slots'] if k]
    assert len(locations)==len(set(locations))==299
    gem=next(g for g in inv['groups'] if g['id']=='f72r2:outer')
    assert len(gem['slots'])==16 and gem['slots'][13] is None
    assert gem['slots'][12]=='f72r2.19' and gem['slots'][14]=='f72r2.20'
    reports=[]
    for stage in ['discovery','confirmation']:
        result=json.loads((ROOT/f'results/zodiac-order-{stage}.json').read_text())
        for path,sha in result['input_hashes'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
        keys={x['locus'] for x in inv['labels'] if x['primary'] and x['stage']==stage}
        groups=[g for g in inv['groups'] if g['stage']==stage and g['kind']=='circle']
        for run in result['runs']:
            texts={f'{x.page}.{x.num}':known_text(x) for x in read(ROOT/'data/raw/ZL3b-n.txt',run['spacing']) if f'{x.page}.{x.num}' in keys}
            pair_values=[]
            for a,b in combinations(groups,2):
                if not a['layout'] or a['layout']!=b['layout'] or a['band']!=b['band'] or a['page']==b['page']:continue
                n=len(a['slots']);scores=[]
                for direction in [1,-1]:
                    for shift in range(n):
                        values=[]
                        for i,k in enumerate(a['slots']):
                            other=b['slots'][(shift+direction*i)%n]
                            if k in texts and other in texts:values.append(sim(texts[k],texts[other]))
                        if len(values)>=3:scores.append(sum(values)/len(values))
                if scores:pair_values.append(max(scores))
            observed=0.;expected=0.;edges=0
            for g in groups:
                valid=[i for i,k in enumerate(g['slots']) if k in texts]
                if len(valid)<3:continue
                buckets=defaultdict(list)
                for i in valid:
                    word=texts[g['slots'][i]]
                    buckets[bin_of(word) if run['null']=='ring_length_bin' else 0].append(i)
                for i,k in enumerate(g['slots']):
                    j=(i+1)%len(g['slots']);kk=g['slots'][j]
                    if k not in texts or kk not in texts:continue
                    if g['id']=='f72r2:outer':assert (k,kk)!=('f72r2.19','f72r2.20')
                    a,b=texts[k],texts[kk];observed+=sim(a,b);edges+=1
                    ka=bin_of(a) if run['null']=='ring_length_bin' else 0
                    kb=bin_of(b) if run['null']=='ring_length_bin' else 0
                    values=[sim(texts[g['slots'][u]],texts[g['slots'][v]]) for u in buckets[ka] for v in buckets[kb] if u!=v]
                    expected+=sum(values)/len(values)
            assert edges==run['n_neighbour_edges']
            assert len(pair_values)==run['n_ring_pairs']
            assert abs(sum(pair_values)/len(pair_values)-run['primary_order']['observed'])<1e-12
            assert abs(observed/edges-run['neighbour_similarity']['observed'])<1e-12
            mc=run['neighbour_similarity'];error=abs(expected/edges-mc['null_mean'])
            assert error<6*mc['null_sd']/math.sqrt(run['randomizations'])+1e-12
            reports.append(dict(stage=stage,spacing=run['spacing'],null=run['null'],
                direct_ring_pairs=len(pair_values),direct_edges=edges,
                direct_order_score=sum(pair_values)/len(pair_values),
                direct_neighbour_score=observed/edges,
                analytic_neighbour_null_mean=expected/edges,mc_mean_error=error))
    out=dict(status='PASS',all_loci_and_missing_slots_accounted=True,
        confirmation_references_unchanged_saved_discovery=True,
        independent_2d_edit_distance_and_loop_order_scores=True,
        neighbour_mc_means_match_analytic_expectations=True,runs=reports)
    (ROOT/'results/audit-004.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':audit()
