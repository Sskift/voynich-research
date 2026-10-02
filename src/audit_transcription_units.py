"""Independent raw-source/selection/profile audit for experiment013."""
from pathlib import Path
from collections import Counter
import hashlib,json,re,sys
from audit_language_constraints import check_corpus,independent_compare
R=Path(__file__).resolve().parents[1]
def read(p):return json.loads((R/p).read_text())
def sha(p):return hashlib.sha256((R/p).read_bytes()).hexdigest()

def independently_parse(file):
    out={};page_loci={};pending=None
    for line in (R/file).read_text().splitlines():
        if pending:
            assert line.startswith('/');out[pending]['raw']+=line[1:].lstrip()
            if out[pending]['raw'].rstrip().endswith('/'):out[pending]['raw']=out[pending]['raw'].rstrip()[:-1]
            else:pending=None
            continue
        m=re.match(r'^<(f[^.>]+\.\d+),([^>]+)>\s*(.*)$',line)
        if not m:continue
        loc,code,raw=m.groups();page=loc.rsplit('.',1)[0]
        out[loc]=dict(raw=raw,code=code);page_loci.setdefault(page,[]).append(loc)
        if raw.rstrip().endswith('/'):out[loc]['raw']=raw.rstrip()[:-1];pending=loc
    assert pending is None
    for locs in page_loci.values():
        for i,l in enumerate(locs):
            nxt=out[locs[i+1]]['code'] if i+1<len(locs) else None
            out[l]['next_code']=nxt;out[l]['eligible']=out[l]['code'] in ('@P0','+P0','*P0') and not (nxt and nxt.startswith(('=P','~P')))
    return out

def words(raw,mode):
    raw=re.sub('<(?:-|~)>','.',raw);raw=re.sub('<[^>]+>','',raw);raw=''.join(raw.split())
    # Split only outside bracket/braced runs, independently of production scanner.
    chunks=re.findall(r'\[[^\]]*\]|\{[^}]*\}|@\d{3};|.',raw);groups=[];buf=''
    for c in chunks:
        if c==',' and mode=='comma_join':continue
        if c in ('.',','):
            if buf:groups.append(buf);buf=''
        else:buf+=c
    if buf:groups.append(buf)
    result=[]
    for i,g in enumerate(groups):
        known=not any(c in g for c in '?[]');atoms=None
        if known:
            flat=g.replace('{','').replace('}','');atoms=[];j=0
            while j<len(flat):
                if flat[j]=='@':a=flat[j:j+5];assert re.fullmatch('@[12][0-9][0-9];',a) and 128<=int(a[1:4])<=255
                else:a=flat[j]
                atoms.append(a);j+=len(a)
        result.append(dict(index=i,raw=g,known=known,atoms=atoms))
    return result

def main():
    s=read('data/transcription-units/selection.json');p=read('results/transcription-units/profiles.json');refs=read('results/language-screen/constraints.json');count=0
    for f in ['data/transcription-units/selection-freeze.json','data/transcription-units/score-freeze.json','data/language-screen/source-freeze.json','results/substitution/solver-freeze.json','results/homophonic/solver-freeze.json']:
        o=read(f);d=next(v for v in o.values() if isinstance(v,dict) and v and all(isinstance(x,str) and re.fullmatch('[a-f0-9]{64}',x) for x in v.values()))
        for path,h in d.items():assert sha(path)==h,path;count+=1
    original={source:independently_parse('data/raw/'+file) for source,file in s['source_files'].items()};rows=0;wordslots=0
    for r in s['rows']:
        for source,entry in r['sources'].items():
            native=original[source].get(r['locus'])
            if native is None:assert entry==dict(eligible=False,reason='missing_locus');continue
            assert all(entry[k]==native[k] for k in ['code','raw','next_code','eligible'])
            if entry['eligible']:
                for mode in ['primary','comma_join']:
                    assert words(native['raw'],mode)==entry['words'][mode];wordslots+=len(entry['words'][mode])
            rows+=1
    matched={};exclusions={}
    for mode in ['primary','comma_join']:
        matched[mode]=[];exclusions[mode]=[]
        for r in s['rows']:
            bad=[source+':'+v['reason'] for source,v in r['sources'].items() if not v['eligible']]
            sizes={source:len(v['words'][mode]) for source,v in r['sources'].items() if v['eligible']}
            if len(set(sizes.values()))>1:bad.append('unequal_word_counts')
            if bad:exclusions[mode].append(dict(locus=r['locus'],split=r['split'],reason=bad,word_counts=sizes));continue
            for i in range(next(iter(sizes.values()))):
                triplet={source:v['words'][mode][i] for source,v in r['sources'].items()}
                if all(w['known'] for w in triplet.values()):matched[mode].append(dict(locus=r['locus'],split=r['split'],word_index=i,words=triplet))
                else:exclusions[mode].append(dict(locus=r['locus'],split=r['split'],word_index=i,reason=['uncertain_word'],raw={source:w['raw'] for source,w in triplet.items()}))
        assert matched[mode]==s['pairs'][mode] and exclusions[mode]==s['paired_exclusions'][mode]
    summaries=0;screens=0
    for mode in ['primary','comma_join']:
        for split in ['train','development','reserve']:
            selected=[r for r in s['rows'] if r['split']==split];pair=[q for q in matched[mode] if q['split']==split]
            for source in s['source_files']:
                eligible=[r['sources'][source] for r in selected if r['sources'][source]['eligible']]
                fw=[w['atoms'] for r in eligible for w in r['words'][mode] if w['known']];pw=[q['words'][source]['atoms'] for q in pair]
                c=p['coverage'][mode][split][source]
                assert c['scope_loci']==len(selected) and c['eligible_loci']==len(eligible) and c['missing_or_ineligible_loci']==len(selected)-len(eligible)
                assert c['known_words']==len(fw) and c['paired_known_words']==len(pw)
                assert c['raw_words']==sum(len(r['words'][mode]) for r in eligible) and c['unknown_words']==c['raw_words']-len(fw)
                assert c['paired_word_fraction']==len(pw)/len(fw) and c['paired_unit_fraction']==sum(map(len,pw))/sum(map(len,fw))
                for scope,ww in [('full',fw),('paired',pw)]:
                    obj=p[scope][mode][split][source];check_corpus(ww,obj,False);summaries+=1+len(obj['blocks'])
                    for language,conditions in obj['language_screen'].items():
                        for condition,gate in conditions.items():
                            bs=[(n,b) for n,meta in refs['metadata'].items() if meta['language']==language and not meta['supplementary'] for b in refs['sources'][n][condition]['blocks']]
                            means=[b['summary']['mean_word_length'] for n,b in bs]
                            admitted=[dict(source=n,block=b['index']) for n,b in bs if b['summary']['max_frequency_opaque_upper']+1e-12>=obj['full']['max_frequency']]
                            inside=min(means)<=obj['full']['mean_word_length']<=max(means)
                            assert gate['mean_word_length_range']==[min(means),max(means)] and gate['native_mean_in_range']==inside and gate['frequency_admitted']==admitted and gate['screen_pass']==(inside and bool(admitted))
                            for (n,b),entry in zip(bs,gate['blocks'],strict=True):assert entry['source']==n and entry['index']==b['index'];independent_compare(obj['full'],b['summary'],entry);screens+=1
            ct=Counter(tuple(''.join(q['words'][ss]['atoms']) for ss in s['source_files']) for q in pair)
            assert p['paired_triplets'][mode][split]==[dict(GC=t[0],ZL=t[1],IT=t[2],tokens=n) for t,n in sorted(ct.items(),key=lambda v:(-v[1],v[0]))]
            for source in ['ZL','IT']:
                ratio=p['paired_ratios'][mode][split][source];a=[len(q['words']['GC']['atoms']) for q in pair];b=[len(q['words'][source]['atoms']) for q in pair]
                assert ratio==dict(total_unit_ratio=sum(b)/sum(a),words=len(pair),shorter=sum(v<u for u,v in zip(a,b)),equal=sum(v==u for u,v in zip(a,b)),longer=sum(v>u for u,v in zip(a,b)))
    assert p['primary_passing_languages']==[la for la,c in p['full']['primary']['train']['ZL']['language_screen'].items() if c['primary']['screen_pass']]
    out=dict(status='PASS',frozen_hashes=count,source_loci=rows,source_word_slots=wordslots,profile_full_and_blocks=summaries,reference_comparisons=screens,primary_passing_languages=p['primary_passing_languages'],limitations='Audits positional source selection and statistics, not physical word alignment, phonemes or German identity.',hashes={f:sha(f) for f in ['src/audit_transcription_units.py','results/transcription-units/profiles.json','data/transcription-units/selection-freeze.json','data/transcription-units/score-freeze.json']})
    (R/'results/audit-010-units.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
