"""Ring-aware label inventory, with original layout slots and uncertainties."""
from pathlib import Path
from collections import Counter
import csv,json,re
from corpus import read

ROOT=Path(__file__).resolve().parents[1]


def known_text(x):
    return '.'.join(x.words) if x.words and all(x.words) else None


def build():
    zl=read(ROOT/'data/raw/ZL3b-n.txt')
    it={f'{x.page}.{x.num}':x for x in read(ROOT/'data/raw/IT2a-n.txt')}
    layout=json.loads((ROOT/'data/zodiac-layout.json').read_text())
    index={};groups=[]
    for page in layout:
        for group in page['groups']:
            slots=[]
            for num in range(group['range'][0],group['range'][1]+1):
                key=f'{page["page"]}.{num}'
                index[key]=dict(group=page['page']+':'+group['id'],slot=len(slots),stage=page['stage'],kind=group['kind'])
                slots.append(key)
                if num==group.get('missing_after'):slots.append(None)
            groups.append(dict(id=page['page']+':'+group['id'],page=page['page'],sign=page['sign'],stage=page['stage'],layout=page['layout'],kind=group['kind'],band=group['id'],slots=slots))
    rows=[]
    for x in zl:
        if x.code[1:]!='Lz' or x.code[0]=='!':continue
        key=f'{x.page}.{x.num}';y=it.get(key)
        assert key in index,key
        reason=[]
        if not y or y.code[1:]!='Lz' or y.code[0]=='!':reason.append('missing_or_wrong_IT_locus')
        a,b=known_text(x),known_text(y) if y else None
        if a is None:reason.append('uncertain_ZL')
        if b is None:reason.append('uncertain_IT')
        if x.words!=(y.words if y else None):reason.append('disagreeing_reading_or_segmentation')
        clock=re.search(r'<!([0-9]{2}:[0-9]{2})',x.raw)
        rows.append(dict(locus=key,page=x.page,num=x.num,zl_raw=x.raw,it_raw=y.raw if y else None,zl_text=a,it_text=b,primary=not reason,exclusions=reason,clock=clock[1] if clock else None,**index[key]))
    assert len(rows)==len(index)==299
    assert len({x['locus'] for x in rows})==299
    bykey={x['locus']:x for x in rows}
    for g in groups:
        g['n_slots']=len(g['slots']);g['n_labels']=sum(k is not None for k in g['slots'])
        g['n_primary']=sum(bykey[k]['primary'] for k in g['slots'] if k)
    summary=dict(n_labels=len(rows),n_primary=sum(x['primary'] for x in rows),n_ZL_known=sum(x['zl_text'] is not None for x in rows),n_groups=len(groups),
        page_counts={p['page']:dict(labels=sum(x['page']==p['page'] for x in rows),primary=sum(x['page']==p['page'] and x['primary'] for x in rows)) for p in layout},
        exclusion_flags=dict(Counter(e for x in rows for e in x['exclusions'])))
    out=dict(summary=summary,groups=groups,labels=rows)
    (ROOT/'results/zodiac-inventory.json').write_text(json.dumps(out,indent=2)+'\n')
    with (ROOT/'results/zodiac-inventory.tsv').open('w') as f:
        writer=csv.writer(f,delimiter='\t');cols=['locus','group','slot','stage','clock','zl_text','it_text','primary','exclusions']
        writer.writerow(cols)
        for x in rows:writer.writerow([','.join(x[k]) if k=='exclusions' else x[k] for k in cols])
    print(json.dumps(summary,indent=2))


if __name__=='__main__':build()
