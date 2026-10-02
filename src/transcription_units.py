"""Experiment013: frozen same-locus transcription units and language constraints."""
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,re
from language_constraints import corpus_result,compare

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/transcription-units';OUT=ROOT/'results/transcription-units'
FILES={'GC':'GC2a-n.txt','ZL':'ZL3b-n.txt','IT':'IT2a-n.txt'}
SPLITS=['train','development','reserve'];MODES=['primary','comma_join']
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def write(p,o):(ROOT/p).write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n')

def parse_source(source):
    lines=(ROOT/'data/raw'/FILES[source]).read_text().splitlines();rows={};pageorder={};i=0
    while i<len(lines):
        m=re.fullmatch(r'<(f[^.>]+)\.(\d+),([^>]+)>\s*(.*)',lines[i]);i+=1
        if not m:continue
        page,num,code,raw=m.groups();locus=f'{page}.{num}'
        while raw.rstrip().endswith('/'):
            assert lines[i].startswith('/');raw=raw.rstrip()[:-1]+lines[i][1:].lstrip();i+=1
        assert locus not in rows
        rows[locus]=dict(page=page,locus=locus,code=code,raw=raw)
        pageorder.setdefault(page,[]).append(locus)
    for page,loci in pageorder.items():
        for j,locus in enumerate(loci):
            nxt=rows[loci[j+1]]['code'] if j+1<len(loci) else None
            row=rows[locus];row['next_code']=nxt
            row['eligible']=row['code'] in {'@P0','+P0','*P0'} and not (nxt and nxt[0] in '=~' and nxt[1]=='P')
    return rows

def split_words(raw,comma_join=False):
    s=raw.replace('<->','.').replace('<~>','.')
    s=re.sub(r'\s+','',re.sub(r'<[^>]*>','',s))
    rawwords=[];word='';brackets=0;braces=0
    for a in s:
        if a=='[':brackets+=1
        elif a==']':brackets-=1
        elif a=='{':braces+=1
        elif a=='}':braces-=1
        assert brackets>=0 and braces>=0
        if a in '.,' and not brackets and not braces:
            if a==',' and comma_join:continue
            if word:rawwords.append(word)
            word=''
        else:word+=a
    if word:rawwords.append(word)
    assert brackets==braces==0
    words=[]
    for i,w in enumerate(rawwords):
        bad=('?' in w or '[' in w or ']' in w)
        item=dict(index=i,raw=w,known=not bad,atoms=None)
        if not bad:
            s=w.replace('{','').replace('}','');atoms=re.findall(r'@\d{3};|.',s)
            assert ''.join(atoms)==s
            for a in atoms:
                if len(a)>1:assert 128<=int(a[1:4])<=255
                else:assert a not in '@;[]{}?.,<>/',(w,a)
            item['atoms']=atoms
        words.append(item)
    return words

def prepare():
    DATA.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    man=read('data/homophonic/manuscript.json');allrows={s:parse_source(s) for s in FILES};rows=[]
    for old in man['rows']:
        entry=dict(locus=old['locus'],page=old['page'],split=old['split'],sources={})
        for s in FILES:
            r=allrows[s].get(old['locus'])
            if r is None:entry['sources'][s]=dict(eligible=False,reason='missing_locus');continue
            e=dict(r);e['reason']=None if r['eligible'] else 'nonbody_or_separate_continuation'
            if r['eligible']:
                e['words']={mode:split_words(r['raw'],mode=='comma_join') for mode in MODES}
            entry['sources'][s]=e
        assert entry['sources']['GC']['eligible'] and entry['sources']['GC']['raw']==old['raw']
        rows.append(entry)
    pairs={};exclusions={}
    for mode in MODES:
        pairs[mode]=[];exclusions[mode]=[]
        for row in rows:
            reason=[]
            for s,r in row['sources'].items():
                if not r['eligible']:reason.append(s+':'+r['reason'])
            sizes={s:len(r['words'][mode]) for s,r in row['sources'].items() if r['eligible']}
            if len(set(sizes.values()))>1:reason.append('unequal_word_counts')
            if reason:
                exclusions[mode].append(dict(locus=row['locus'],split=row['split'],reason=reason,word_counts=sizes));continue
            for i in range(next(iter(sizes.values()))):
                ws={s:r['words'][mode][i] for s,r in row['sources'].items()}
                if not all(w['known'] for w in ws.values()):
                    exclusions[mode].append(dict(locus=row['locus'],split=row['split'],word_index=i,reason=['uncertain_word'],raw={s:w['raw'] for s,w in ws.items()}));continue
                pairs[mode].append(dict(locus=row['locus'],split=row['split'],word_index=i,words=ws))
    selection=dict(status='FROZEN_BEFORE_UNIT_PROFILES',rows=rows,pairs=pairs,paired_exclusions=exclusions,source_files=FILES)
    write('data/transcription-units/selection.json',selection)
    paths=['notes/experiment-013-design.md','notes/experiment-013-parser-rules.md','src/transcription_units.py','src/language_constraints.py','data/homophonic/manuscript.json','data/language-screen/source-freeze.json','results/language-screen/constraints.json','data/transcription-units/selection.json','data/raw/IVTFF_format.txt','data/raw/transcr.html','data/raw/sta-aaa.html']+['data/raw/'+f for f in FILES.values()]
    write('data/transcription-units/selection-freeze.json',dict(status='SELECTION_FROZEN_BEFORE_COMPARATIVE_PROFILES',files={p:sha(p) for p in paths},selected_GC_loci=len(rows),meaning='Source word positions/components only. No phonetic or semantic assignment.'))
    print('FROZEN',len(rows),'GC loci; no comparative profile calculated.')

def collect(selection,source,split,mode,paired):
    if paired:return [p['words'][source]['atoms'] for p in selection['pairs'][mode] if p['split']==split]
    return [w['atoms'] for row in selection['rows'] if row['split']==split and row['sources'][source]['eligible'] for w in row['sources'][source]['words'][mode] if w['known']]

def screen(stats,refs):
    result={}
    for language in ['la','de','it']:
        names=[n for n,m in refs['metadata'].items() if m['language']==language and not m['supplementary']]
        result[language]={}
        for sensitivity in ['primary','apostrophe_split']:
            bs=[(n,b) for n in names for b in refs['sources'][n][sensitivity]['blocks']]
            means=[b['summary']['mean_word_length'] for n,b in bs]
            admitted=[dict(source=n,block=b['index']) for n,b in bs if b['summary']['max_frequency_opaque_upper']+1e-12>=stats['max_frequency']]
            within=min(means)<=stats['mean_word_length']<=max(means)
            result[language][sensitivity]=dict(mean_word_length_range=[min(means),max(means)],native_mean_in_range=within,frequency_admitted=admitted,screen_pass=within and bool(admitted),blocks=[dict(source=n,index=b['index'],**compare(stats,b['summary'])) for n,b in bs])
    return result

def score():
    frozen=read('data/transcription-units/selection-freeze.json')
    for p,h in frozen['files'].items():assert sha(p)==h,p
    selection=read('data/transcription-units/selection.json');refs=read('results/language-screen/constraints.json');full={};paired={};coverage={}
    for mode in MODES:
        full[mode]={};paired[mode]={};coverage[mode]={}
        for split in SPLITS:
            full[mode][split]={};paired[mode][split]={};coverage[mode][split]={}
            subset=[r for r in selection['rows'] if r['split']==split]
            for s in FILES:
                fw=collect(selection,s,split,mode,False);pw=collect(selection,s,split,mode,True)
                f=corpus_result(fw,False);p=corpus_result(pw,False)
                f['language_screen']=screen(f['full'],refs);p['language_screen']=screen(p['full'],refs)
                full[mode][split][s]=f;paired[mode][split][s]=p
                eligible=[r['sources'][s] for r in subset if r['sources'][s]['eligible']]
                coverage[mode][split][s]=dict(scope_loci=len(subset),eligible_loci=len(eligible),missing_or_ineligible_loci=len(subset)-len(eligible),raw_words=sum(len(r['words'][mode]) for r in eligible),known_words=len(fw),unknown_words=sum(not w['known'] for r in eligible for w in r['words'][mode]),paired_known_words=len(pw),paired_word_fraction=len(pw)/len(fw),paired_unit_fraction=p['full']['units']/f['full']['units'])
    for split in SPLITS:
        assert full['primary'][split]['GC']['full']==refs['native'][split]['full']
    triplets={};ratios={}
    for mode in MODES:
        triplets[mode]={};ratios[mode]={}
        for split in SPLITS:
            ps=[p for p in selection['pairs'][mode] if p['split']==split]
            ct=Counter(tuple(''.join(p['words'][s]['atoms']) for s in FILES) for p in ps)
            triplets[mode][split]=[dict(GC=t[0],ZL=t[1],IT=t[2],tokens=n) for t,n in sorted(ct.items(),key=lambda v:(-v[1],v[0]))]
            ratios[mode][split]={}
            for s in ['ZL','IT']:
                base=paired[mode][split]['GC']['full'];dest=paired[mode][split][s]['full']
                ratios[mode][split][s]=dict(total_unit_ratio=dest['units']/base['units'],words=len(ps),shorter=sum(len(p['words'][s]['atoms'])<len(p['words']['GC']['atoms']) for p in ps),equal=sum(len(p['words'][s]['atoms'])==len(p['words']['GC']['atoms']) for p in ps),longer=sum(len(p['words'][s]['atoms'])>len(p['words']['GC']['atoms']) for p in ps))
    gate=full['primary']['train']['ZL']['language_screen'];passing=[la for la,m in gate.items() if m['primary']['screen_pass']]
    out=dict(status='TRANSCRIPTION_UNIT_SCREEN_NOT_A_DECIPHERMENT',selection_freeze_sha256=sha('data/transcription-units/selection-freeze.json'),primary_gate='full primary-space training ZL, reference primary apostrophe condition',primary_passing_languages=passing,reference_sha256=sha('results/language-screen/constraints.json'),full=full,paired=paired,coverage=coverage,paired_ratios=ratios,paired_triplets=triplets,limitations='Positional matches with equal source word counts can retain compensating split/join differences; not a verified character mapping. Three transcriptions concern one manuscript.')
    write('results/transcription-units/profiles.json',out)
    summary=dict(primary_passing_languages=passing,full_train={s:{k:full['primary']['train'][s]['full'][k] for k in ['units','words','mean_word_length','max_frequency','most_frequent','types']} for s in FILES},paired_train={s:{k:paired['primary']['train'][s]['full'][k] for k in ['units','words','mean_word_length','max_frequency','most_frequent','types']} for s in FILES},ratios=ratios['primary']['train'],coverage=coverage['primary']['train'])
    write('results/transcription-units/summary.json',summary);print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','score']);args=p.parse_args()
    prepare() if args.stage=='prepare' else score()
