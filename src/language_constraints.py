"""Frozen experiment012 diagnostics; no language values assigned to Voynich glyphs."""
from collections import Counter
from pathlib import Path
import hashlib,json,random,math,unicodedata as U

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/language-screen'
OUT=ROOT/'results/language-screen'
def read(p):return json.loads((ROOT/p).read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,obj): (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def opaque(s):return any(U.category(c)=='Co' or c.isdigit() or c in '~&⁊' for c in s)

def blocks(words):
    out=[];start=0;size=0
    for i,w in enumerate(words):
        size+=len(w)
        if size>=5000:
            out.append(dict(start_word=start,end_word=i+1,words=words[start:i+1]))
            start=i+1;size=0
    return out,dict(start_word=start,words=len(words)-start,units=size)

def summarize(words,check_opaque=True):
    c=Counter(a for w in words for a in w);lengths=Counter(map(len,words));n=c.total();nw=len(words)
    ranks=sorted(c.values(),reverse=True); cumulative=[];total=0
    for k in ranks:total+=k;cumulative.append(total/n)
    oc={a:k for a,k in c.items() if check_opaque and opaque(a)}
    clear={a:k for a,k in c.items() if a not in oc}
    return dict(units=n,words=nw,types=len(c),counts=dict(sorted(c.items())),word_lengths=dict(sorted(lengths.items())),mean_word_length=n/nw,max_frequency=max(c.values())/n,most_frequent=sorted(a for a,k in c.items() if k==max(c.values())),opaque_counts=oc,opaque_units=sum(oc.values()),max_clear_frequency=max(clear.values(),default=0)/n,max_opaque_frequency=max(oc.values(),default=0)/n,max_frequency_opaque_upper=(max(clear.values(),default=0)+sum(oc.values()))/n,ranked_cumulative=cumulative,entropy_bits=-sum(k/n*math.log2(k/n) for k in c.values()))

def deficit(cipher,reference):
    q=cipher['ranked_cumulative'];p=reference['ranked_cumulative']
    return max([0]+[(q[k] if k<len(q) else 1)-(p[k] if k<len(p) else 1) for k in range(max(len(q),len(p)))])

def ks(a,b):
    ca=cb=0;out=0
    aa={int(k):v for k,v in a['word_lengths'].items()};bb={int(k):v for k,v in b['word_lengths'].items()}
    for k in sorted(set(aa)|set(bb)):
        ca+=aa.get(k,0);cb+=bb.get(k,0);out=max(out,abs(ca/a['words']-cb/b['words']))
    return out

def compare(q,p):
    return dict(cumulative_deficit=deficit(q,p),word_length_KS=ks(q,p),mean_length_ratio=p['mean_word_length']/q['mean_word_length'],top_frequency_admits=q['max_frequency']<=p['max_frequency']+1e-12,top_frequency_admits_with_opaque=q['max_frequency']<=p['max_frequency_opaque_upper']+1e-12)

def corpus_result(words,check_opaque=True,reference=None):
    bs,tail=blocks(words);full=summarize(words,check_opaque)
    result=dict(full=full,discarded_block_tail=tail,blocks=[])
    if reference:result['versus_native_H']=compare(reference,full)
    for i,b in enumerate(bs):
        s=summarize(b['words'],check_opaque);item=dict(index=i,start_word=b['start_word'],end_word=b['end_word'],summary=s)
        if reference:item['versus_native_H']=compare(reference,s)
        result['blocks'].append(item)
    return result

def main():
    OUT.mkdir(exist_ok=True)
    freeze=read('data/language-screen/source-freeze.json')
    for p,h in freeze['files'].items():assert digest(ROOT/p)==h,p
    previous={}
    for name in ['substitution','homophonic']:
        f=read(f'results/{name}/solver-freeze.json')
        # Freeze schemas differ; inspect the common file-hash map explicitly.
        d=f.get('files',f.get('hashes',f.get('frozen_hashes',f.get('inputs'))))
        if d is None:
            d=next(v for v in f.values() if isinstance(v,dict) and v and all(isinstance(h,str) and len(h)==64 for h in v.values()))
        for p,h in d.items():assert digest(ROOT/p)==h,p
        previous[name]=len(d)
    man=read('data/homophonic/manuscript.json');native={};nativewords={};nativeexclusions={}
    for split in ['train','development','reserve']:
        ws=[];excluded=0
        for r in man['rows']:
            if r['split']!=split:continue
            word=[]
            for a in r['atoms']+[' ']:
                if a==' ':
                    if word:
                        if '?' in word:excluded+=1
                        else:ws.append(word)
                    word=[]
                else:word.append(a)
        nativewords[split]=ws;native[split]=corpus_result(ws,False);nativeexclusions[split]=excluded
    h=native['train']['full'];sources={};sourcewords={};metadata={}
    for entry in freeze['sources']:
        obj=read(entry['file']);name=obj['name'];metadata[name]={k:obj[k] for k in ['language','supplementary','diplomatic']};sources[name]={}
        for sensitivity,key in [('primary','words'),('apostrophe_split','words_apostrophe_split')]:
            ws=[w for r in obj['records'] for w in r[key]]
            sources[name][sensitivity]=corpus_result(ws,True,h)
            if sensitivity=='primary':sourcewords[name]=ws
    controls=[];seed=20261012
    for si,name in enumerate(n for n in sources if not metadata[n]['supplementary']):
        bs,_=blocks(sourcewords[name])
        for bi,b in enumerate(bs[:2]):
            plain=b['words'];ps=summarize(plain,False);vocab=sorted(ps['counts'])
            for mode in ['bijection','homophonic']:
                controlseed=seed+100*si+10*bi+(mode=='homophonic');rng=random.Random(controlseed)
                codes=list(range(len(vocab)*(4 if mode=='homophonic' else 1)));rng.shuffle(codes)
                key={};groups={};ci=0
                for a in vocab:
                    count=rng.randint(1,4) if mode=='homophonic' else 1
                    groups[a]=[f'c{n}' for n in codes[ci:ci+count]];ci+=count
                    key.update({symbol:a for symbol in groups[a]})
                cipher=[[rng.choice(groups[a]) for a in w] for w in plain];cs=summarize(cipher,False)
                exact=([len(w) for w in cipher]==[len(w) for w in plain] and [[key[a] for a in w] for w in cipher]==plain and deficit(cs,ps)<1e-12 and cs['max_frequency']<=ps['max_frequency']+1e-12)
                assert exact
                control=dict(source=name,block=bi,mode=mode,seed=controlseed,key=key,cipher_words=cipher,plaintext=ps,ciphertext=cs,own_paired_invariants_pass=exact)
                control['cross_work_reference']=[]
                for other in sources:
                    if other==name or metadata[other]['supplementary'] or metadata[other]['language']!=metadata[name]['language']:continue
                    for ref in sources[other]['primary']['blocks']:
                        d=compare(cs,ref['summary'])
                        control['cross_work_reference'].append(dict(source=other,block=ref['index'],**d,naive_exact_majorization_rejects=d['cumulative_deficit']>1e-12))
                controls.append(control)
    decision={}
    for language in ['la','de','it']:
        decision[language]={}
        names=[n for n in sources if metadata[n]['language']==language and not metadata[n]['supplementary']]
        for sensitivity in ['primary','apostrophe_split']:
            entries=[dict(source=n,**b) for n in names for b in sources[n][sensitivity]['blocks']]
            means=[e['summary']['mean_word_length'] for e in entries]
            admits=[dict(source=e['source'],block=e['index']) for e in entries if e['versus_native_H']['top_frequency_admits_with_opaque']]
            inrange=min(means)<=h['mean_word_length']<=max(means)
            maxupper=max(e['summary']['max_frequency_opaque_upper'] for e in entries)
            decision[language][sensitivity]=dict(works=names,blocks=len(entries),mean_word_length_range=[min(means),max(means)],native_mean_in_range=inrange,top_frequency_admitting_blocks=admits,top_frequency_range=[min(e['summary']['max_frequency'] for e in entries),max(e['summary']['max_frequency'] for e in entries)],max_opaque_upper=maxupper,screen_pass=bool(admits) and inrange,minimum_length_expansion_to_reach_block_range=min(means)/h['mean_word_length'],one_unit_model_minimum_mass_needing_other_treatment=max(0,h['max_frequency']-maxupper))
    # Paired edition changes are descriptive, not a decoded plaintext comparison.
    reading=read('data/language-screen/brant-reading.json');diplomatic=read('data/language-screen/brant-diplomatic.json')
    paired=[]
    for a,b in zip(reading['records'],diplomatic['records'],strict=True):
        assert a['locator']==b['locator']
        na=sum(map(len,a['words']));nb=sum(map(len,b['words']))
        paired.append(dict(chapter=a['locator'],reading_units=na,diplomatic_units=nb,reading_words=len(a['words']),diplomatic_words=len(b['words']),unit_expansion_ratio=na/nb))
    pair=dict(chapters=paired,total_unit_expansion_ratio=sum(p['reading_units'] for p in paired)/sum(p['diplomatic_units'] for p in paired),reading_to_diplomatic_word_ratio=sum(p['reading_words'] for p in paired)/sum(p['diplomatic_words'] for p in paired),note='Same edition and work, not independent samples. Corrections, expansions, glyph normalization and occasionally spacing all contribute; not a pure abbreviation estimate.')
    out=dict(status='CONSTRAINT_SCREEN_NOT_A_DECIPHERMENT',source_freeze_sha256=digest(DATA/'source-freeze.json'),script_sha256=digest(Path(__file__)),previous_solver_freezes_verified=previous,native=native,native_excluded_unreadable_words=nativeexclusions,metadata=metadata,sources=sources,decisions=decision,brant_paired=pair,controls_summary=dict(controls=len(controls),all_paired_invariants_pass=all(c['own_paired_invariants_pass'] for c in controls),cross_work_comparisons=sum(len(c['cross_work_reference']) for c in controls),cross_work_naive_rejections=sum(r['naive_exact_majorization_rejects'] for c in controls for r in c['cross_work_reference'])))
    write('constraints.json',out);write('synthetic-controls.json',controls)
    print(json.dumps(dict(native_H={k:h[k] for k in ['units','words','max_frequency','mean_word_length']},decisions=decision,controls=out['controls_summary'],brant_unit_expansion=pair['total_unit_expansion_ratio']),ensure_ascii=False,indent=2))

if __name__=='__main__':main()
