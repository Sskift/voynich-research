"""Independent experiment012 audit. Does not import extraction/scoring modules."""
from pathlib import Path
from collections import Counter
from html import unescape
import json,re,hashlib,math,random,unicodedata as U
import xml.etree.ElementTree as E

R=Path(__file__).resolve().parents[1]; D=R/'data/language-screen'; RAW=R/'data/raw/language-screen'
def read(p):return json.loads((R/p).read_text())
def sha(p):return hashlib.sha256((R/p).read_bytes()).hexdigest()
def norm(s):return ' '.join(s.split())
def htmltext(s):return unescape(re.sub('<[^>]+>','',re.sub(r'<br\s*/?>','\n',s,flags=re.I)))
def joinlines(s):return re.sub(r'[-⸗\u00ad]\s*\n\s*(?=\w|[\ue000-\uf8ff])','',s)

def source_texts(name):
    if name=='macer':
        s=(RAW/'macer.html').read_text();out=[]
        for row in re.findall(r'<tr[^>]*>(.*?)</tr>',s,re.S|re.I):
            cells=re.findall(r'<td\b[^>]*>(.*?)</td>',row,re.S|re.I)
            if len(cells)!=2:continue
            # Last verse p in each plant cell has an omitted HTML closing tag.
            for p in re.findall(r'<p\b[^>]*>(.*?)(?:</p>|(?=<p\b)|$)',cells[0],re.S|re.I):
                if re.search(r'<b\b',p,re.I):continue
                t=htmltext(p)
                if norm(t):out.append(norm(re.sub(r'^\s*\d+\s+','',t)))
        assert len(out)==2269
        return out
    if name=='aquinas':
        s=(RAW/'aquinas-q1.html').read_text(encoding='latin1');out=[]
        for p in re.findall(r'<P TITLE="I q\. 1 [^"]*">(.*?)</P>',s,re.S):
            out.append(norm(htmltext(re.sub(r'<SPAN CLASS="ref">.*?</SPAN>','',p,flags=re.S))))
        assert len(out)==66
        return out
    if name=='tepl':
        s=(RAW/'tepl-47465.html').read_text();s=s[s.index('<a id="Das_erst_capitel"'):s.index('</div> <!-- end of Text div')]
        assert len(re.findall('<blockquote>',s))>0
        s=re.sub(r'<blockquote>.*?</blockquote>','',s,flags=re.S)
        out=[]
        for p in re.findall(r'<p\b[^>]*>(.*?)</p>',s,re.S):
            p=re.sub(r'<span class="pagenum">.*?</span>','',p,flags=re.S)
            p=re.sub(r'<a class="fnanchor[^"]*"[^>]*>.*?</a>','',p,flags=re.S)
            t=norm(joinlines(re.sub(r'\[Z\.\s*\d+\]','',htmltext(p))))
            if t:out.append(t)
        assert len(out)==92
        return out
    if name=='boccaccio':
        s=(RAW/'boccaccio-21424.html').read_text();out=[]
        for ident,p in re.findall(r'<p id="id(\d+)"[^>]*>(.*?)</p>',s,re.S):
            n=int(ident)
            if 20<=n<=174 and n not in (84,142):out.append(norm(joinlines(htmltext(p))))
        assert len(out)==95
        return out
    if name=='dante':
        s=(RAW/'dante-1012.txt').read_text();start=s.index('  Nel mezzo del cammin di nostra vita');last='  l’amor che move il sole e l’altre stelle.'
        s=s[start:s.index(last,start)+len(last)]
        out=[line.strip() for line in s.splitlines() if line.strip() and not re.fullmatch(r'\s*(?:(?:Inferno|Purgatorio|Paradiso) • Canto [IVXLCDM]+|INFERNO|PURGATORIO|PARADISO)\s*',line)]
        assert len(out)==14233
        return out
    ns='{http://www.tei-c.org/ns/1.0}';dip=name.endswith('diplomatic')
    tree=E.parse(RAW/('locher-1497.xml' if name.startswith('locher') else 'brant-1494.xml'));out=[]
    def traverse(e):
        tag=e.tag.split('}')[-1]
        if tag in {'note','rdg','rdgGrp'}:return []
        if tag in {'choice','app'}:
            if tag=='app':matches=e.findall('.//'+ns+'lem')
            else:
                names=('orig','sic','abbr') if dip else ('reg','corr','expan')
                matches=[child for label in names for child in list(e) if child.tag==ns+label]
            assert len(matches)==1
            return traverse(matches[0])
        parts=[e.text or '']
        for child in list(e):parts.extend(traverse(child));parts.append(child.tail or '')
        if tag in {'l','p','lg'}:parts.append('\n')
        return parts
    for d in tree.findall('.//'+ns+'body/'+ns+'div'):
        if not d.get('{http://www.w3.org/XML/1998/namespace}id','').startswith('B'):continue
        chunks=[''.join(traverse(a)) for a in d.iter(ns+'ab') if a.get('type')=='mainText']
        if chunks:out.append(norm(joinlines('\n'.join(chunks))))
    assert len(out)==(110 if name.startswith('locher') else 112)
    return out

def independent_words(text,lang,diplomatic,split):
    s=U.normalize('NFC',text.lower())
    if not diplomatic:
        for a,b in [('ﬀ','ff'),('ﬁ','fi'),('ﬂ','fl'),('ﬃ','ffi'),('ﬄ','ffl'),('ﬅ','st'),('ﬆ','st'),('ſ','s'),('ꝛ','r')]:s=s.replace(a,b)
        if lang=='la':
            for a,b in [('j','i'),('v','u'),('æ','ae'),('œ','oe')]:s=s.replace(a,b)
        if lang=='it':s=''.join(c for c in U.normalize('NFD',s) if not U.category(c).startswith('M'))
        s=U.normalize('NFC',s)
    groups=[];w=[]
    for i,c in enumerate(s):
        category=U.category(c)
        if category.startswith('M') and w:w[-1]+=c
        elif category.startswith('L') or category=='Co' or c in '&⁊':w.append(c)
        elif (c.isdigit() or c=='~') and w and i+1<len(s) and U.category(s[i+1]).startswith('L'):w.append(c)
        elif not split and c in "'’‘ʼ":continue
        else:
            if w:groups.append(w)
            w=[]
    if w:groups.append(w)
    return groups

def check_stats(words,s,opaque):
    c=Counter(a for w in words for a in w);lc=Counter(map(len,words));n=sum(c.values());nw=len(words)
    assert dict(c)==s['counts']
    assert {str(k):v for k,v in lc.items()}==s['word_lengths']
    assert (n,nw,len(c))==(s['units'],s['words'],s['types'])
    assert s['most_frequent']==sorted(a for a,v in c.items() if v==max(c.values()))
    oc={a:k for a,k in c.items() if opaque and any(U.category(ch)=='Co' or ch.isdigit() or ch in '~&⁊' for ch in a)}
    assert s['opaque_counts']==oc and s['opaque_units']==sum(oc.values())
    assert s['max_frequency']==max(c.values())/n and s['mean_word_length']==n/nw
    clearmax=max((v for a,v in c.items() if a not in oc),default=0)
    assert s['max_clear_frequency']==clearmax/n
    assert s['max_opaque_frequency']==max(oc.values(),default=0)/n
    assert s['max_frequency_opaque_upper']==(clearmax+sum(oc.values()))/n
    expected=[];acc=0
    for v in sorted(c.values(),reverse=True):acc+=v;expected.append(acc/n)
    assert expected==s['ranked_cumulative']
    assert abs(s['entropy_bits']+sum(v/n*math.log2(v/n) for v in c.values()))<1e-12

def independent_compare(c,p,entry):
    cc=sorted(c['counts'].values(),reverse=True);pc=sorted(p['counts'].values(),reverse=True)
    # Integer cumulative numerators avoid assuming paired invariants from floats.
    csum=psum=0;d=0
    for k in range(max(len(cc),len(pc))):
        if k<len(cc):csum+=cc[k]
        if k<len(pc):psum+=pc[k]
        d=max(d,csum/c['units']-psum/p['units'])
    all_lengths=sorted({int(x) for x in c['word_lengths']}|{int(x) for x in p['word_lengths']})
    ca=pa=0;ks=0
    for k in all_lengths:
        ca+=c['word_lengths'].get(str(k),0);pa+=p['word_lengths'].get(str(k),0)
        ks=max(ks,abs(ca/c['words']-pa/p['words']))
    assert abs(entry['cumulative_deficit']-d)<1e-14
    assert abs(entry['word_length_KS']-ks)<1e-14
    assert entry['mean_length_ratio']==p['mean_word_length']/c['mean_word_length']
    assert entry['top_frequency_admits']==(c['max_frequency']<=p['max_frequency']+1e-12)
    assert entry['top_frequency_admits_with_opaque']==(c['max_frequency']<=p['max_frequency_opaque_upper']+1e-12)

def check_corpus(words,obj,opaque,reference=None):
    check_stats(words,obj['full'],opaque)
    if reference:independent_compare(reference,obj['full'],obj['versus_native_H'])
    start=0
    for i,b in enumerate(obj['blocks']):
        end=b['end_word'];size=sum(map(len,words[start:end]))
        assert b['index']==i and b['start_word']==start
        assert size>=5000 and size-len(words[end-1])<5000
        check_stats(words[start:end],b['summary'],opaque)
        if reference:independent_compare(reference,b['summary'],b['versus_native_H'])
        start=end
    tail=obj['discarded_block_tail'];assert tail==dict(start_word=start,words=len(words)-start,units=sum(map(len,words[start:]))) and tail['units']<5000

def main():
    x=read('results/language-screen/constraints.json');freeze=read('data/language-screen/source-freeze.json')
    for p,h in freeze['files'].items():assert sha(p)==h,p
    assert sha('data/language-screen/source-freeze.json')==x['source_freeze_sha256']
    assert sha('src/language_constraints.py')==x['script_sha256']
    # Original IVTFF lines, independent of saved atoms and of old tokenizer imports.
    lines=(R/'data/raw/GC2a-n.txt').read_text().splitlines();original={};i=0
    while i<len(lines):
        m=re.match(r'<(f[^>,]+),[^>]+>\s*(.*)',lines[i]);i+=1
        if not m:continue
        locus,raw=m.groups()
        while raw.rstrip().endswith('/'):
            assert lines[i].startswith('/');raw=raw.rstrip()[:-1]+lines[i][1:].lstrip();i+=1
        original[locus]=raw
    nativewords={k:[] for k in x['native']};unread=Counter();man=read('data/homophonic/manuscript.json')
    for row in man['rows']:
        assert original[row['locus']]==row['raw']
        s=row['raw'].replace('<->','.').replace('<~>','.')
        s=re.sub(r'\s+','',re.sub(r'<[^>]*>','',s)).replace('{','').replace('}','')
        for word in re.split('[.,]',s):
            if not word:continue
            if '?' in word or '[' in word:unread[row['split']]+=1;continue
            units=[];j=0
            while j<len(word):
                code=re.match(r'@\d{3};',word[j:])
                a=code[0] if code else word[j];units.append(a);j+=len(a)
            nativewords[row['split']].append(units)
    for split,ws in nativewords.items():
        assert unread[split]==x['native_excluded_unreadable_words'][split]
        check_corpus(ws,x['native'][split],False)
    h=x['native']['train']['full'];words={};records_checked=0;source_summaries=0
    for source in freeze['sources']:
        name=source['name'];data=read(source['file']);assert sha(data['raw_file'])==data['raw_sha256']
        texts=source_texts(name);assert texts==[r['text'] for r in data['records']],name
        for is_split,key,mode in [(False,'words','primary'),(True,'words_apostrophe_split','apostrophe_split')]:
            ws=[]
            for r in data['records']:
                independent=independent_words(r['text'],data['language'],data['diplomatic'],is_split)
                assert independent==r[key],(name,r['locator']);ws.extend(independent)
            check_corpus(ws,x['sources'][name][mode],True,h)
            source_summaries+=len(x['sources'][name][mode]['blocks'])+1
            if not is_split:words[name]=ws
        records_checked+=len(texts)
    controls=read('results/language-screen/synthetic-controls.json');cross_total=reject=0;order=list(n for n in x['sources'] if not x['metadata'][n]['supplementary'])
    for c in controls:
        src=c['source'];b=x['sources'][src]['primary']['blocks'][c['block']];plain=words[src][b['start_word']:b['end_word']]
        check_stats(plain,c['plaintext'],False);check_stats(c['cipher_words'],c['ciphertext'],False)
        assert [[c['key'][a] for a in w] for w in c['cipher_words']]==plain
        assert list(map(len,c['cipher_words']))==list(map(len,plain))
        nc=sorted(c['ciphertext']['counts'].values(),reverse=True);np=sorted(c['plaintext']['counts'].values(),reverse=True)
        assert all(sum(nc[:k])<=sum(np[:k]) for k in range(1,max(len(nc),len(np))+1))
        seed=20261012+100*order.index(src)+10*c['block']+(c['mode']=='homophonic');assert seed==c['seed'];rng=random.Random(seed)
        vocab=sorted(c['plaintext']['counts']);codes=list(range(len(vocab)*(4 if c['mode']=='homophonic' else 1)));rng.shuffle(codes)
        groups={};off=0;key={}
        for a in vocab:
            n=rng.randint(1,4) if c['mode']=='homophonic' else 1
            groups[a]=[f'c{i}' for i in codes[off:off+n]];off+=n;key.update({k:a for k in groups[a]})
        assert key==c['key'];assert [[rng.choice(groups[a]) for a in w] for w in plain]==c['cipher_words']
        for entry in c['cross_work_reference']:
            p=x['sources'][entry['source']]['primary']['blocks'][entry['block']]['summary']
            independent_compare(c['ciphertext'],p,entry)
            assert entry['naive_exact_majorization_rejects']==(entry['cumulative_deficit']>1e-12)
            cross_total+=1;reject+=entry['naive_exact_majorization_rejects']
    assert (len(controls),cross_total,reject)==(24,x['controls_summary']['cross_work_comparisons'],x['controls_summary']['cross_work_naive_rejections'])
    for lang,modes in x['decisions'].items():
        for mode,declared in modes.items():
            bs=[(name,b) for name in declared['works'] for b in x['sources'][name][mode]['blocks']]
            lengths=[b['summary']['mean_word_length'] for n,b in bs];freq=[b['summary']['max_frequency'] for n,b in bs];upper=max(b['summary']['max_frequency_opaque_upper'] for n,b in bs)
            admitted=[dict(source=n,block=b['index']) for n,b in bs if b['summary']['max_frequency_opaque_upper']+1e-12>=h['max_frequency']]
            assert declared['top_frequency_admitting_blocks']==admitted
            assert declared['mean_word_length_range']==[min(lengths),max(lengths)]
            assert declared['blocks']==len(bs) and declared['top_frequency_range']==[min(freq),max(freq)]
            within=min(lengths)<=h['mean_word_length']<=max(lengths)
            assert declared['native_mean_in_range']==within and declared['screen_pass']==(within and bool(admitted))
            assert declared['max_opaque_upper']==upper
            assert declared['minimum_length_expansion_to_reach_block_range']==min(lengths)/h['mean_word_length']
            assert declared['one_unit_model_minimum_mass_needing_other_treatment']==max(0,h['max_frequency']-upper)
    paired=x['brant_paired'];aa=read('data/language-screen/brant-reading.json')['records'];bb=read('data/language-screen/brant-diplomatic.json')['records']
    for a,b,p in zip(aa,bb,paired['chapters'],strict=True):
        assert a['locator']==b['locator']==p['chapter']
        for record,prefix in [(a,'reading'),(b,'diplomatic')]:
            assert sum(map(len,record['words']))==p[prefix+'_units'];assert len(record['words'])==p[prefix+'_words']
        assert p['unit_expansion_ratio']==p['reading_units']/p['diplomatic_units']
    assert paired['total_unit_expansion_ratio']==sum(p['reading_units'] for p in paired['chapters'])/sum(p['diplomatic_units'] for p in paired['chapters'])
    assert paired['reading_to_diplomatic_word_ratio']==sum(p['reading_words'] for p in paired['chapters'])/sum(p['diplomatic_words'] for p in paired['chapters'])
    previous={}
    for name in ['substitution','homophonic']:
        f=read(f'results/{name}/solver-freeze.json');d=next(v for v in f.values() if isinstance(v,dict) and v and all(isinstance(h,str) and len(h)==64 for h in v.values()))
        for p,hsh in d.items():assert sha(p)==hsh,p
        previous[name]=len(d)
    out=dict(status='PASS',source_records_independently_reselected=records_checked,source_full_and_block_summaries=source_summaries,native_loci_rechecked=len(man['rows']),synthetic_controls_regenerated=len(controls),cross_work_comparisons=cross_total,cross_work_exact_rejections=reject,paired_brant_chapters=len(paired['chapters']),source_frozen_hashes=len(freeze['files']),old_solver_frozen_hashes=previous,limitations='Checks extraction, units, arithmetic and exact paired invariants; cannot authenticate every historical spelling, estimate whole-language probabilities, prove a cipher family, or decipher the manuscript.',audited_hashes={p:sha(p) for p in ['src/audit_language_constraints.py','data/language-screen/source-freeze.json','results/language-screen/constraints.json','results/language-screen/synthetic-controls.json']})
    (R/'results/audit-009.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n');print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
