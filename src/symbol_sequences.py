"""Case/escape-preserving inventories and a frozen cross-page sequence test."""
from pathlib import Path
from collections import defaultdict
import argparse,csv,hashlib,json,re
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FILES={'ZL':'ZL3b-n.txt','IT':'IT2a-n.txt','GC':'GC2a-n.txt','RF':'RF1b-e.txt'}


def raw_loci(source,page):
    result={};lines=(ROOT/'data/raw'/FILES[source]).read_text().splitlines();i=0
    while i<len(lines):
        m=re.fullmatch(r'<'+re.escape(page)+r'\.(\d+),([^>]+)>\s*(.*)',lines[i]);i+=1
        if not m:continue
        num,code,raw=m.groups()
        while raw.rstrip().endswith('/'):
            assert lines[i].startswith('/')
            raw=raw.rstrip()[:-1]+lines[i][1:].lstrip();i+=1
        result[int(num)]=dict(code=code,raw=raw)
    return result


def clean(raw):
    raw=raw.replace('<->','.').replace('<~>','.')
    return re.sub(r'\s+','',re.sub(r'<[^>]*>','',raw))


def layout_groups(raw):
    s=clean(raw);groups=[];current='';depth=0
    for c in s:
        if c in '[{':depth+=1
        if c in ']}':depth-=1
        if c in '.,' and depth==0:groups.append(current);current=''
        else:current+=c
    groups.append(current);assert all(groups) and depth==0
    return groups


def alternatives(text):
    m=re.search(r'\[([^\[\]]+)\]',text)
    if m:
        return sorted({s for v in m[1].split(':') for s in alternatives(text[:m.start()]+v+text[m.end():])})
    if '?' in text:return []
    return [text.replace('{','').replace('}','')]


def native_atoms(text):
    s=clean(text).replace('.','').replace(',','').replace('{','').replace('}','')
    atoms=re.findall(r'@\d{3};|\[[^\]]*\]|.',s)
    assert ''.join(atoms)==s
    return atoms


def snapshot(paths):
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def inventory():
    sources={};templates={}
    for source in FILES:
        raw=raw_loci(source,'f57v')[3]['raw'];original=layout_groups(raw)
        groups=native_atoms(raw) if source=='GC' else list(original)
        if source=='IT':
            assert len(groups)==69 and groups[-2:]==['r','n']
            # This preserves an unresolved correspondence span; it is NOT a glyph merge.
            cells=groups[:67]+[None]
        else:
            assert len(groups)==68,(source,len(groups));cells=groups
        rows=[]
        for i,text in enumerate(cells):
            rows.append(dict(cycle=i//17+1,position=i%17+1,raw=text,
                alternatives=alternatives(text) if text is not None else [],
                unresolved_source_fragments=groups[-2:] if text is None else [],
                encoded_atoms=native_atoms(text) if text is not None else []))
        bypos=[[rows[cycle*17+pos]['alternatives'] for cycle in range(4)] for pos in range(17)]
        template=[sorted({v for seq in sequences for v in seq}) for sequences in bypos]
        changes=[dict(position=pos+1,cycle_readings=sequences) for pos,sequences in enumerate(bypos) if any(seq!=sequences[0] for seq in sequences[1:])]
        sources[source]=dict(raw_locus=raw,apparent_spacing_group_count=len(original),
            native_encoded_atom_count=len(native_atoms(raw)),comparison_cells=rows,
            variation_positions=changes,provenance_note='Dependent reference combining ZL and GC, descriptive only.' if source=='RF' else 'Native source, differences retained; no semantic equivalence claimed.')
        if source!='RF':templates[source]=template
    paths=['notes/experiment-009-design.md']+['data/raw/'+f for f in FILES.values()]
    out=dict(status='FROZEN_BEFORE_F66R_COLUMN_INSPECTION_OR_SCORING',input_hashes=snapshot(paths),
        sources=sources,templates=templates,allograph_merges_applied_to_corpus=False,
        label='Four corresponding17-position spans, not a phonetic alphabet.',
        IT_last_span_note='IT has69 fragments. Final r/n is an unresolved span in position17 of cycle4; neither value is used there.')
    (ROOT/'results/f57v-symbol-inventory.json').write_text(json.dumps(out,indent=2)+'\n')
    with (ROOT/'results/f57v-symbol-inventory.tsv').open('w') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(['source','cycle','position','raw','alternatives','unresolved_fragments'])
        for source,obj in sources.items():
            for r in obj['comparison_cells']:writer.writerow([source,r['cycle'],r['position'],r['raw'],','.join(r['alternatives']),','.join(r['unresolved_source_fragments'])])
    for source,obj in sources.items():
        print(source,'spacing groups',obj['apparent_spacing_group_count'],'encoded atoms',obj['native_encoded_atom_count'],'variations',obj['variation_positions'])


def summarize(scores):
    return dict(observed=float(scores[0]),null_mean=float(scores[1:].mean()),null_sd=float(scores[1:].std()),
        p_greater=float((1+np.count_nonzero(scores[1:]>=scores[0]-1e-12))/len(scores)),
        null_interval_95=[float(v) for v in np.quantile(scores[1:],[.025,.975])])


def target():
    inv_path=ROOT/'results/f57v-symbol-inventory.json';inv=json.loads(inv_path.read_text())
    assert inv['status']=='FROZEN_BEFORE_F66R_COLUMN_INSPECTION_OR_SCORING'
    for p,h in inv['input_hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    layout=json.loads((ROOT/'data/f66r-symbol-layout.json').read_text())
    runs=[];reps=9999;seed=20261009
    for source in ['ZL','IT','GC']:
        loci=raw_loci(source,'f66r');rows=[]
        for i,num in enumerate(layout['loci_top_to_bottom']):
            item=loci.get(num);raw=item['raw'] if item else None
            value=clean(raw).replace('.','').replace(',','').replace('{','').replace('}','') if item else None
            reason=[]
            if not item:reason.append('absent_locus')
            elif item['code'][1:]!='L0':reason.append('unexpected_locus_type')
            if value is not None and any(c in value for c in '?[]'):reason.append('uncertain_target_reading')
            if reason:value=None
            rows.append(dict(locus=f'f66r.{num}',position=i,raw=raw,value=value,exclusions=reason,
                encoded_atoms=native_atoms(value) if value else []))
        template=inv['templates'][source];known=[r for r in rows if r['value'] is not None]
        positions=np.array([r['position'] for r in known]);values=[r['value'] for r in known];n=len(values)
        vocab=sorted(set(values));lookup={s:i for i,s in enumerate(vocab)};original=np.array([lookup[v] for v in values])
        match=np.array([[v in ts for ts in template] for v in vocab],dtype=np.int64)
        phases=[dict(start=s,direction=d) for d in [1,-1] for s in range(17)]
        predictions=np.array([[(p['start']+p['direction']*i)%17 for i in positions] for p in phases])
        rng=np.random.default_rng(seed);orders=np.empty((reps+1,n),dtype=np.int64);orders[0]=original
        for rep in range(1,reps+1):orders[rep]=rng.permutation(original)
        score=np.array([match[orders,pred].sum(axis=1) for pred in predictions]).T
        best=int(np.argmax(score[0]));maxima=score.max(axis=1)/n
        allsymbols=set(v for ts in template for v in ts)
        summary=dict(source=source,primary=source=='ZL',randomizations=reps,seed=seed,n_slots=len(rows),n_known=n,
            order_free_template_coverage=sum(v in allsymbols for v in values),
            primary_order=summarize(maxima),observed_max_matches=int(score[0,best]),best_phase=phases[best],
            all_phases=[dict(**p,matches=int(score[0,i])) for i,p in enumerate(phases)],
            missing_positions=[r['position'] for r in rows if r['value'] is None],rows=rows,
            best_matches=[dict(locus=r['locus'],target=r['value'],template_position=int(predictions[best,i]+1),
                template_options=template[predictions[best,i]],matches=bool(match[original[i],predictions[best,i]])) for i,r in enumerate(known)])
        runs.append(summary);print(source,{k:summary[k] for k in ['n_known','order_free_template_coverage','observed_max_matches','best_phase','primary_order']},flush=True)
    out=dict(frozen_input_hashes=snapshot(['results/f57v-symbol-inventory.json','data/f66r-symbol-layout.json','notes/experiment-009-design.md']),
        question='Same17-position template in the separate34-entry column?',
        independence_note='Reserved score inspection, previously available corpus. Transcriptions share manuscript evidence; RF not a test source.',runs=runs)
    (ROOT/'results/f66r-symbol-prediction.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['inventory','target']);a=p.parse_args()
    inventory() if a.stage=='inventory' else target()
