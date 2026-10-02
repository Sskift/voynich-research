"""Independent LM, cipher, objective and native-output audit for experiment014."""
from pathlib import Path
from collections import Counter
import hashlib,json,math,re,unicodedata as U
import numpy as np
from audit_transcription_units import independently_parse,words as native_words
R=Path(__file__).resolve().parents[1]
def read(p):return json.loads((R/p).read_text())
def sha(p):return hashlib.sha256((R/p).read_bytes()).hexdigest()
def close(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),(a.keys(),b.keys())
        for k in a:close(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):close(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10),(a,b)
    else:assert a==b,(a,b)
def gram_counts(rows,n,limit):
    c=Counter()
    for row in rows:
        for q in zip(*(row[i:] for i in range(n))):
            if all(v<=limit for v in q):c[q]+=1
    return c

class Auditor:
    def __init__(self,version):
        self.version=version;self.folder='german-components'+('-v2' if version==2 else '');self.data='data/'+self.folder;self.out='results/'+self.folder
        self.model=read(self.data+'/language.json');self.alpha=self.model['alphabet'];self.L=len(self.alpha);self.A=self.L+1;self.objectives=0
        bfile='data/german-components-v2/brant-base.json' if version==2 else 'data/language-screen/brant-reading.json'
        tfile='data/german-components-v2/tepl-base.json' if version==2 else 'data/language-screen/tepl.json'
        brant=read(bfile);tepl=read(tfile)
        assert self.model['source']==bfile and self.model['source_sha256']==sha(bfile)
        def encode_words(ws):
            seq=[self.L]
            for w in ws:
                seq.extend(self.alpha.index(c) if c in self.alpha else self.L+1 for c in w);seq.append(self.L)
            return seq
        rows=[encode_words(r['words']) for r in brant['records']]
        previous=None
        for n in range(1,5):
            ct=gram_counts(rows,n,self.L);array=np.zeros((self.A,)*n)
            for q,k in ct.items():array[q]=k
            if n==1:prob=(array+.1)/(sum(ct.values())+.1*self.A);close(prob.tolist(),self.model['unigram'])
            else:prob=(array+10*np.expand_dims(previous,0))/(array.sum(axis=n-1,keepdims=True)+10)
            assert np.allclose(prob.sum(axis=-1),1,atol=1e-14);previous=prob
        self.logp=np.log(prob).reshape(-1);stored=np.fromfile(R/self.data/'conditional4.bin',dtype=np.float64)
        assert np.allclose(self.logp,stored,rtol=0,atol=1e-13) and sha(self.data+'/conditional4.bin')==self.model['lm_sha256']
        vocab=Counter(''.join(w) for r in brant['records'] for w in r['words'] if all(c in self.alpha for c in w));assert dict(vocab)==self.model['word_counts']
        six=Counter()
        for row in rows:
            s=''.join(' ' if c==self.L else '□' if c>self.L else self.alpha[c] for c in row)
            six.update(s[i:i+6] for i in range(len(s)-5) if '□' not in s[i:i+6])
        assert dict(six)==self.model['sixgram_counts']
        assert sum(c<self.L for r in rows for c in r)==self.model['known_letter_units']
        assert sum(c>self.L for r in rows for c in r)==self.model['opaque_unit_masks']
        ws=[w for r in tepl['records'] for w in r['words']];cut=math.floor(.6*len(ws));plain=read(self.data+'/plain-control.json')
        assert plain['source']==tfile and plain['source_sha256']==sha(tfile) and plain['split_word_index']==cut and plain['total_word_tokens']==len(ws)
        for split,slice_ in [('train',ws[:cut]),('test',ws[cut:])]:assert [encode_words(slice_[i:i+10]) for i in range(0,len(slice_),10)]==plain['plain_'+split]
        self.plain=plain;self.vocab=vocab;self.six=six
    def probabilities(self,rows,key,G):
        c=Counter(v for r in rows for v in r if v<G);totals=Counter()
        for v,n in c.items():totals[key[v]]+=n
        assert len(c)==G
        return [c[v]/totals[key[v]] for v in range(G)]
    def lm_sum(self,rows,key,G):
        c=gram_counts(rows,4,G);score=0
        for q,n in c.items():
            offset=0
            for v in q:offset=offset*self.A+key[v]
            score+=n*self.logp[offset]
        return score,sum(c.values())
    def metrics(self,rows,key,prob,G):
        total,mass=self.lm_sum(rows,key,G);ec=0;hits=sn=wh=wn=0
        for r in rows:
            ec+=sum(math.log(prob[c]) for c in r if c<G)
            text=''.join('□' if c>G else ' ' if c==G else self.alpha[key[c]] for c in r)
            for i in range(len(text)-5):
                gram=text[i:i+6]
                if '□' not in gram:sn+=1;hits+=gram in self.six
            for word in text.split():
                if len(word)>=3 and '□' not in word:wn+=1;wh+=word in self.vocab
        return dict(conditional4=float(total/mass),emission_per_window=ec/mass,combined=float((total+ec)/mass),windows=mass,sixgram_seen=hits/sn if sn else None,sixgram_hits=hits,sixgram_windows=sn,word_seen=wh/wn if wn else None,word_hits=wh,complete_words_ge3=wn)
    def fit(self,name,rows,G):
        f=read(self.out+'/'+name+'.frozen.json');ct=Counter(v for r in rows for v in r if v<G);assert f['train_counts']==[ct[i] for i in range(G)]
        assert f['model_sha256']==sha(self.data+'/conditional4.bin')
        rawkeys=[]
        for line in (R/self.out/(name+'.keys.txt')).read_text().splitlines():
            x=line.split();rawkeys.append(dict(restart=int(x[0]),objective=float(x[1]),key=list(map(int,x[2:]))))
        assert f['restarts']==rawkeys and len(rawkeys)==16
        for k in rawkeys:
            assert len(k['key'])==G+1 and k['key'][-1]==self.L and all(0<=v<self.L for v in k['key'][:-1])
            score,mass=self.lm_sum(rows,k['key'],G);prob=self.probabilities(rows,k['key'],G);objective=(score+sum(n*math.log(prob[c]) for c,n in ct.items()))/mass
            close(float(objective),k['objective']);assert f['train_windows']==mass;self.objectives+=1
        assert f['selected']==max(rawkeys,key=lambda k:k['objective']);close(f['emission_probabilities'],self.probabilities(rows,f['selected']['key'],G))
        expected=gram_counts(rows,4,G);actual={}
        for line in (R/self.out/(name+'.grams.txt')).read_text().splitlines():
            a=list(map(int,line.split()));actual[tuple(a[:4])]=a[4]
        assert actual==expected
        assert list(map(int,(R/self.out/(name+'.counts.txt')).read_text().split()))==[G]+f['train_counts']
        return f
    def accuracy(self,d,key):
        G=d['glyph_count'];alln=pred=right=wn=wr=0
        for row,tr in zip(d['test'],d['plain_test'],strict=True):
            group=[]
            for c,t in zip(row,tr,strict=True):
                alln+=t<self.L
                if c<G:pred+=1;right+=key[c]==t
                if c==G:
                    if group and all(c<G for c,t in group):wn+=1;wr+=all(key[c]==t for c,t in group)
                    group=[]
                else:group.append((c,t))
        return dict(letter_accuracy=right/pred,correct_letters=right,predicted_letters=pred,total_letters=alln,predicted_fraction=pred/alln,word_accuracy=wr/wn,correct_words=wr,complete_words=wn,key_accuracy=sum(key[i]==d['true_key'][i] for i in range(G))/G)
    def controls(self):
        records=read(self.out+f'/calibration-v{self.version}.json');cfg=read(self.data+f'/config-v{self.version}.json');L=self.L
        for trial,entry in enumerate(records['controls']):
            d=read(self.data+'/'+entry['name']+'.json');seed=cfg['seed']+trial;assert seed==d['seed'];rng=np.random.default_rng(seed)
            weights=np.array(self.model['unigram'][:L]);weights/=weights.sum();fraction=(d['N_requested']-L)*weights;extras=np.floor(fraction).astype(int)
            for i in np.argsort(-(fraction-extras),kind='stable')[:d['N_requested']-L-int(extras.sum())]:extras[i]+=1
            sizes=extras+1;assert sizes.tolist()==d['homophones_per_letter'];permutation=list(map(int,rng.permutation(d['N_requested'])));groups=[];truth={};i=0
            for letter,count in enumerate(sizes):
                group=permutation[i:i+count];i+=count;groups.append(group)
                for sym in group:truth[sym]=letter
            raw={}
            for split in ['train','test']:
                raw[split]=[]
                for row in self.plain['plain_'+split]:
                    coded=[]
                    for letter in row:
                        if letter==L:coded.append(' ');continue
                        if rng.random()<d['mask_probability']:coded.append('?');continue
                        group=groups[letter]
                        if len(group)==1 or d['primary_probability'] is not None and rng.random()<d['primary_probability']:sym=group[0]
                        else:sym=int(rng.choice(group if d['primary_probability'] is None else group[1:]))
                        coded.append(sym)
                    raw[split].append(coded)
            counts=Counter(sym for row in raw['train'] for sym in row if isinstance(sym,int));labels=sorted(counts,key=lambda k:(-counts[k],k));G=len(labels);lookup={sym:i for i,sym in enumerate(labels)}
            assert labels==d['cipher_labels'] and G==d['glyph_count'] and d['true_key']==[truth[c] for c in labels]+[L]
            for split in ['train','test']:
                assert [[G if c==' ' else lookup.get(c,G+1) for c in row] for row in raw[split]]==d[split]
                assert d['plain_'+split]==self.plain['plain_'+split]
                assert sum(isinstance(c,int) and c not in lookup for row in raw[split] for c in row)==d[split+'_unseen_positions']
            fit=self.fit(entry['name'],d['train'],G);key=fit['selected']['key'];close(self.accuracy(d,key),entry['accuracy']);close(self.metrics(d['test'],key,fit['emission_probabilities'],G),entry['test_metrics'])
            close(self.metrics(d['test'],d['true_key'],self.probabilities(d['train'],d['true_key'],G),G),entry['true_key_metrics'])
            for restart,acc in zip(fit['restarts'],entry['restart_accuracy'],strict=True):close(dict(restart=restart['restart'],objective=restart['objective'],**self.accuracy(d,restart['key'])),acc)
            assert sha(self.data+'/'+entry['name']+'.json')==entry['data_sha256'] and sha(self.out+'/'+entry['name']+'.frozen.json')==entry['fit_sha256']
        assert records['status']==('PASS' if all(e['accuracy']['letter_accuracy']>=.95 for e in records['controls']) else 'FAIL')
        return records['status']

def main():
    verified=0
    for path in ['data/transcription-units/selection-freeze.json','data/transcription-units/score-freeze.json','data/language-screen/source-freeze.json','results/substitution/solver-freeze.json','results/homophonic/solver-freeze.json','results/german-components/failed-v1-freeze.json','data/german-components-v2/normalization-freeze.json','results/german-components-v2/solver-freeze.json']:
        o=read(path);d=next(v for v in o.values() if isinstance(v,dict) and v and all(isinstance(h,str) and re.fullmatch('[a-f0-9]{64}',h) for h in v.values()))
        for p,h in d.items():assert sha(p)==h,p;verified+=1
    for src,dest in [('brant-reading','brant-base'),('tepl','tepl-base')]:
        old=read('data/language-screen/'+src+'.json');new=read('data/german-components-v2/'+dest+'.json')
        for a,b in zip(old['records'],new['records'],strict=True):
            assert a['locator']==b['locator'] and a['raw']==b['raw'] and a['text']==b['text']
            for key in ['words','words_apostrophe_split']:
                expected=[]
                for word in a[key]:
                    ww=[]
                    for c in word:
                        c={'\ue42c':'a','\ue644':'o','\ue72b':'u','\ue72d':'u'}.get(c,c);ww.append(''.join(x for x in U.normalize('NFD',c) if not U.category(x).startswith('M')))
                    expected.append(ww)
                assert expected==b[key] and [len(w) for w in a[key]]==[len(w) for w in b[key]]
    assert 'E42C' in (R/'data/raw/titus-a.html').read_text() and all(c in (R/'data/raw/titus-u.html').read_text() for c in ['E72B','E72D'])
    assert 'E644' in (R/'data/raw/ocrd-character-table.html').read_text()
    v1=Auditor(1);assert v1.controls()=='FAIL';v2=Auditor(2);assert v2.controls()=='PASS';G=63;L=v2.L
    man=read(v2.data+'/manuscript.json');oldman=read('data/german-components/manuscript.json');assert man==oldman
    original=independently_parse('data/raw/ZL3b-n.txt');lookup={g:i for i,g in enumerate(man['glyphs'])};splits={k:[] for k in man['coverage']};scope=read('data/transcription-units/selection.json')['rows'];expected=[(r['locus'],r['split']) for r in scope if r['sources']['ZL']['eligible']]+[(loc,'score_reserve_f1r') for loc,r in original.items() if loc.startswith('f1r.') and r['eligible']]
    assert [(r['locus'],r['split']) for r in man['rows']]==expected
    for row in man['rows']:
        assert original[row['locus']]['raw']==row['raw'] and original[row['locus']]['eligible'];ws=native_words(row['raw'],'primary');assert ws==row['words'];encoded=[G]
        for w in ws:encoded.extend([lookup.get(c,G+1) for c in w['atoms']] if w['known'] else [G+1]);encoded.append(G)
        assert encoded==row['encoded'];splits[row['split']].append(encoded)
    for split,c in man['coverage'].items():
        rows=[r for r in man['rows'] if r['split']==split];expected=dict(loci=len(rows),known_input_atoms=sum(len(w['atoms']) for r in rows for w in r['words'] if w['known']),unreadable_whole_words=sum(not w['known'] for r in rows for w in r['words']),training_seen_atoms=sum(c<G for r in rows for c in r['encoded']),unseen_known_atoms=sum(a not in lookup for r in rows for w in r['words'] if w['known'] for a in w['atoms']))
        assert c==expected
    out=read(v2.out+'/attack.json');assert out['complete'];fit=v2.fit('voynich',splits['train'],G);key=fit['selected']['key'];prob=fit['emission_probabilities'];assert key==out['key'];close(prob,out['emission_probabilities']);target='score_reserve_f1r'
    for split,rows in splits.items():close(v2.metrics(rows,key,prob,G),out['metrics'][split])
    for restart,record in zip(fit['restarts'],out['restart_reserve_metrics'],strict=True):close(dict(restart=restart['restart'],objective=restart['objective'],key=restart['key'],metrics=v2.metrics(splits[target],restart['key'],v2.probabilities(splits['train'],restart['key'],G),G)),record)
    cfg=read(v2.data+'/config-v2.json');rng=np.random.default_rng(cfg['seed']+2000);randomscores=[]
    for _ in range(999):
        perm=rng.permutation(L);candidate=[int(perm[k]) for k in key[:-1]]+[L];score,n=v2.lm_sum(splits[target],candidate,G);randomscores.append(float(score/n))
    close(randomscores,out['random_label_scores']);assert out['random_label_p']==(1+sum(s>=out['metrics'][target]['conditional4'] for s in randomscores))/1000
    for n,record in enumerate(out['nulls']):
        assert record['trial']==n and record['shuffle_seed']==cfg['seed']+3000+n;rng=np.random.default_rng(record['shuffle_seed'])
        def permute(rows):
            result=[]
            for row in rows:
                positions=np.flatnonzero(np.array(row)<G);letters=rng.permutation(np.array(row)[positions]);rr=list(row)
                for i,c in zip(positions,letters):rr[int(i)]=int(c)
                result.append(rr)
            return result
        train=permute(splits['train']);test=permute(splits[target]);f=v2.fit('null-'+str(n),train,G);close(v2.metrics(test,f['selected']['key'],f['emission_probabilities'],G),record['reserve_metrics']);close(f['selected']['objective'],record['train_objective']);assert sha(v2.out+'/null-'+str(n)+'.frozen.json')==record['fit_sha256']
    assert len(out['nulls'])==9 and out['null_conditional_p']==(1+sum(n['reserve_metrics']['conditional4']>=out['metrics'][target]['conditional4'] for n in out['nulls']))/10
    dec=read(v2.out+'/voynich-decoded.json');assert dec['status']=='NOT_A_TRANSLATION' and dec['glyph_to_letter']==dict(zip(man['glyphs'],[v2.alpha[k] for k in key[:-1]]))
    for row,d in zip(man['rows'],dec['rows'],strict=True):assert d==dict(locus=row['locus'],split=row['split'],source=row['raw'],output=''.join(' ' if c==G else '□' if c>G else v2.alpha[key[c]] for c in row['encoded']))
    report=dict(status='PASS',frozen_hashes_verified=verified,conditional_probability_entries=[v1.A**4,v2.A**4],synthetic_cipher_conditions_regenerated=6,known_cipher_restart_accuracies_checked=96,all_train_objectives_checked=v1.objectives+v2.objectives,native_loci_rebuilt=len(man['rows']),random_label_scores=999,shuffled_ciphers_regenerated=9,v1_calibration='FAIL',v2_calibration='PASS',limitations='Validates normalization, models, recovery metrics, raw input and search scores, not global optimum, historical spelling, a German identity or a decipherment. V2 controls were development after V1.',hashes={p:sha(p) for p in ['src/audit_german_components.py',v2.out+'/solver-freeze.json',v2.out+'/attack.json',v2.out+'/voynich-decoded.json','notes/experiment-014-character-source-correction.md','data/raw/ocrd-character-table.html']})
    (R/'results/audit-010-decoder.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
