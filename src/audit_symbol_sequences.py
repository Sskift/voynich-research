"""Independent direct-membership recount, including all saved permutation p values."""
from pathlib import Path
from collections import Counter
import hashlib,json,re
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def stripped(s):
    s=s.replace('<->','.').replace('<~>','.')
    s=re.sub(r'<[^>]*>','',s)
    return re.sub(r'\s+','',s)


def extract(path,page,num):
    prefix=f'<{page}.{num},'
    lines=[s for s in path.read_text().splitlines() if s.startswith(prefix)]
    assert len(lines)==1
    return lines[0].split('>',1)[1].strip()


def audit():
    inv=json.loads((ROOT/'results/f57v-symbol-inventory.json').read_text())
    out=json.loads((ROOT/'results/f66r-symbol-prediction.json').read_text())
    for obj,key in [(inv,'input_hashes'),(out,'frozen_input_hashes')]:
        for name,sha in obj[key].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    paths={'ZL':'ZL3b-n.txt','IT':'IT2a-n.txt','GC':'GC2a-n.txt','RF':'RF1b-e.txt'}
    for source in paths:
        raw=extract(ROOT/'data/raw'/paths[source],'f57v',3)
        assert raw==inv['sources'][source]['raw_locus']
        s=stripped(raw)
        groups=re.split(r'[.,]',s)
        if source=='GC':groups=re.findall(r'@\d{3};|.',s.replace('.','').replace(',',''))
        if source=='IT':
            assert len(groups)==69 and groups[-2:]==['r','n'];groups=groups[:67]+[None]
        assert len(groups)==68
        records=inv['sources'][source]['comparison_cells']
        assert groups==[r['raw'] for r in records]
        if source in inv['templates']:
            rebuilt=[]
            for pos in range(17):
                options=set()
                for cycle in range(4):
                    cell=groups[cycle*17+pos]
                    if cell is None or '?' in cell:continue
                    # These particular source cells have only a whole-cell alternative.
                    options.update(cell[1:-1].split(':') if cell.startswith('[') and cell.endswith(']') else [cell])
                rebuilt.append(sorted(options))
            assert rebuilt==inv['templates'][source]
    assert inv['templates']['ZL'][8]==['f','p']
    assert inv['templates']['ZL'][3]==inv['templates']['ZL'][11]==['r']
    assert not inv['allograph_merges_applied_to_corpus']
    reports=[]
    for run in out['runs']:
        source=run['source'];known=[]
        for row in run['rows']:
            num=int(row['locus'].split('.')[1]);assert 16<=num<=49
            raw=extract(ROOT/'data/raw'/paths[source],'f66r',num)
            assert raw==row['raw']
            value=stripped(raw).translate(str.maketrans('','','.,{}'))
            if any(c in value for c in '?[]'):
                assert row['value'] is None
            else:
                assert row['value']==value;known.append((row['position'],value))
        assert len(known)==run['n_known']
        template=[set(x) for x in inv['templates'][source]]
        phases=[(p['start'],p['direction']) for p in run['all_phases']]
        predicted=[[template[(s+d*pos)%17] for pos,_ in known] for s,d in phases]
        values=[v for _,v in known]
        def scores(v):return [sum(x in allowed for x,allowed in zip(v,pattern)) for pattern in predicted]
        observed=scores(values)
        assert observed==[p['matches'] for p in run['all_phases']]
        maximum=max(observed);assert maximum==run['observed_max_matches']
        rng=np.random.default_rng(run['seed']);maxima=[];original=np.arange(len(values))
        counter=Counter(values)
        for rep in range(run['randomizations']):
            indices=rng.permutation(original);permuted=[values[i] for i in indices]
            if rep==0:assert Counter(permuted)==counter
            maxima.append(max(scores(permuted)))
        direct_p=(1+sum(m>=maximum for m in maxima))/(1+len(maxima))
        mean=sum(maxima)/len(maxima)/len(values)
        assert abs(direct_p-run['primary_order']['p_greater'])<1e-12
        assert abs(mean-run['primary_order']['null_mean'])<1e-12
        reports.append(dict(source=source,all_phases_checked=len(phases),known_targets=len(known),
            observed_matches=maximum,direct_p=direct_p,permutations_recounted=len(maxima)))
    result=dict(status='PASS',frozen_hashes_checked=True,raw_source_units_and_templates_checked=True,
        target_loci_and_uncertainty_masks_checked=True,all_observed_and_randomized_scores_recomputed=True,
        note='Direct Python membership counts, independent of the experiment lookup-table score matrix; same specified RNG seed.',runs=reports)
    (ROOT/'results/audit-006.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':audit()
