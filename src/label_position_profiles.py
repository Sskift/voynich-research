"""Exploratory positional controls for the previously observed label differences."""
from collections import defaultdict
from pathlib import Path
import json
import numpy as np
from corpus import read
from label_rules import transform,association

ROOT=Path(__file__).resolve().parents[1]


def features(word): return np.array([word.startswith('o'),word.startswith('qo'),word.endswith('m')],dtype=float)


def profile(source,spacing,labels):
    pools=defaultdict(list)
    for x in read(ROOT/'data/raw'/source,spacing):
        if not x.eligible: continue
        for i,w in enumerate(x.words):
            if not w:continue
            tags=['all']
            if i==0: tags.append('paragraph_first' if x.paragraph_start else 'other_line_first')
            if i==len(x.words)-1: tags.append('line_last')
            if 0<i<len(x.words)-1: tags.append('line_internal')
            for tag in tags: pools[tag,x.page,len(w)].append(w)
    output=[]
    for mode in ['tokens','types']:
        for tag in ['all','paragraph_first','other_line_first','line_internal','line_last']:
            label_values,body_values,all_values=[],[],[]
            excluded=[];included=[]
            for x in labels:
                pool=pools[tag,x['page'],len(x['word'])]
                general=pools['all',x['page'],len(x['word'])]
                if not pool:excluded.append(x['locus']);continue
                if mode=='types':pool=sorted(set(pool));general=sorted(set(general))
                a=features(x['word']);b=np.mean([features(w) for w in pool],axis=0)
                c=np.mean([features(w) for w in general],axis=0)
                label_values.append(a);body_values.append(b);all_values.append(c)
                included.append(dict(locus=x['locus'],word=x['word'],pool_size=len(pool),
                    label_features=a.tolist(),position_expected_features=b.tolist(),all_expected_features=c.tolist()))
            av,bv,cv=map(lambda v:np.mean(v,axis=0),[label_values,body_values,all_values])
            output.append(dict(control=mode,position=tag,n_matched=len(included),excluded=excluded,
                features=[dict(feature=f,label_rate=float(a),position_rate=float(b),
                               all_position_rate_same_support=float(c),label_minus_position=float(a-b))
                          for f,a,b,c in zip(['initial_o','initial_qo','final_m'],av,bv,cv)],
                matched=included))
    return dict(source=source,spacing=spacing,profiles=output)


def diagnostics(rows,fit):
    body=defaultdict(set)
    for x in read(ROOT/'data/raw/ZL3b-n.txt'):
        if x.page in fit['fitting_body_pages'] and x.code[1]=='P' and x.code[0]!='!':body[x.page].update(w for w in x.words if w)
    pages=sorted(body);pindex={p:i for i,p in enumerate(pages)}
    target=np.array([pindex[x['page']] for x in rows])
    specs=[('prefix','o','qo'),('prefix','o',''),('suffix','m','r'),('suffix','m','l'),('suffix','m','iin'),('suffix','y','dy')]
    rules=[dict(side='identity',old='',new='',id='identity')]+[
        dict(side=s,old=a,new=b,id=f'{s}:{a}->{b or "∅"}') for s,a,b in specs]
    hits=np.array([[[transform(x['word'],r) in body[p] for p in pages] for x in rows] for r in rules],dtype=np.int8)
    scores=association(hits-hits[0:1],list(range(len(rows))),target,rows)
    return [dict(**r,changed_types=sum(transform(w,r)!=w for w in {x['word'] for x in rows}),
                 training_own_page_hits=int(hits[i,np.arange(len(rows)),target].sum()),
                 training_centered_gain_per_label=float(scores[i])) for i,r in enumerate(rules)]


if __name__=='__main__':
    labels=[x for x in json.loads((ROOT/'results/pharma-labels.json').read_text()) if x['subtype']=='Lf' and x['agreed_single_word']]
    fit=json.loads((ROOT/'results/label-rule-fit.json').read_text())
    training=[x for x in labels if x['locus'] in fit['training_loci']]
    result=dict(exploratory=True,runs=[profile(s,m,labels) for s,m in [('ZL3b-n.txt','split'),('IT2a-n.txt','split'),('ZL3b-n.txt','join')]],
                fixed_training_rule_diagnostics=diagnostics(training,fit))
    (ROOT/'results/label-position-profiles.json').write_text(json.dumps(result,indent=2)+'\n')
    for run in result['runs']:
        for x in run['profiles']:
            if x['control']=='tokens': print(run['source'],run['spacing'],x['position'],x['n_matched'],x['features'])
    print('Fixed training rule diagnostics:',result['fixed_training_rule_diagnostics'])
