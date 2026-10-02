"""Independent015 audit: source counts, model table, ciphers, all restart scores."""
from pathlib import Path
from collections import Counter
from itertools import product
import hashlib,json,math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; D=ROOT/'data/variable-units-v1'; O=ROOT/'results/variable-units-v1'
ALPH='abcdefghiklmnopqrstuxyz'; L=23; A=24
UNITS=list(ALPH)+[a+b for a in ALPH for b in ALPH]
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-10),(a,b)
def verify(a,b):
    assert a.keys()==b.keys(),(a.keys(),b.keys())
    for k,v in a.items():
        if isinstance(v,(float,int)):close(v,b[k])
        else:assert v==b[k],(k,v,b[k])

def main():
    assert read('results/variable-units-v1/calibration.json')['status'] in {'PASS','FAIL'}
    freezes=['results/substitution/solver-freeze.json','results/homophonic/solver-freeze.json','data/language-screen/source-freeze.json','data/german-components-v2/normalization-freeze.json','results/german-components-v2/solver-freeze.json','data/variable-units-v1/calibration-freeze.json']
    frozen=0
    for f in freezes:
        ob=read(f); hs=ob.get('hashes',ob.get('input_hashes',ob.get('files')))
        assert hs is not None,f
        for p,h in hs.items():assert sha(p)==h,p;frozen+=1
    words=[]
    for n in range(1,8):words.extend((ROOT/f'data/latin/caesar-gallic-{n}.txt').read_text().split())
    lang=read('data/variable-units-v1/language.json');assert Counter(words)==Counter(lang['word_counts']);assert len(words)==lang['word_tokens']
    ct=[Counter() for _ in range(4)]
    for w in words:
        seq=[L,L,L]+[ALPH.index(c) for c in w]+[L]
        for i in range(3,len(seq)):
            for n in range(1,5):ct[n-1][tuple(seq[i-n+1:i+1])]+=1
    levels=[]
    for n in range(1,5):
        arr=np.zeros((A,)*n)
        for gram,count in ct[n-1].items():arr[gram]=count
        if n==1:p=(arr+.1)/(sum(ct[0].values())+.1*A)
        else:p=(arr+10*levels[-1][None,...])/(arr.sum(axis=-1)[...,None]+10)
        levels.append(p)
    lm=np.fromfile(D/'word-conditional4.bin').reshape((A,)*4)
    assert np.allclose(np.log(levels[-1]),lm,rtol=0,atol=1e-14)
    length=Counter(map(len,words));prior=[0]+[(length[n]+20*.8**(n-1))/(len(words)+100) for n in range(1,129)]
    assert np.allclose(prior,lang['length_prior'],rtol=1e-15,atol=1e-20)
    # Independent forward implementation: shift one possible next letter at a time.
    state=np.zeros((A,)*3);state[L,L,L]=1;z=[]
    for n in range(129):
        z.append(float((state*levels[-1][:,:,:,L]).sum()))
        nxt=np.zeros_like(state)
        for d in range(L):nxt[:,:,d]=(state*levels[-1][:,:,:,d]).sum(axis=0)
        state=nxt
    assert np.allclose(z,lang['word_probability_mass_by_length'],rtol=1e-12,atol=0)
    terms=np.fromfile(D/'length-terms.bin').reshape(65,129);termchecks=0
    for g in range(1,65):
        for n in range(g,2*g+1):
            close(terms[g,n],math.log(prior[n])-math.log(sum(prior[g:2*g+1]))-math.log(z[n])-math.log(math.comb(g,n-g)));termchecks+=1
    # Enumerate the generative path space instead of decoding cipher words.
    toy=[]
    for g in range(1,4):
        pi=lambda n:.2*.8**(n-1);den=sum(pi(n) for n in range(g,2*g+1))
        for units,phi in [(['a','b','aa','ab','ba','bb'],[1]*6),(['a','b','aa','ab','ba','bb','a'],[.3,1,1,1,1,1,.7]),(['a','b'],[1,1])]:
            generative=Counter();fail=0.
            for lens in product([1,2],repeat=g):
                n=sum(lens);cp=pi(n)/den/math.comb(g,n-g)
                for letters in product('ab',repeat=n):
                    mass=cp*math.prod(.6 if c=='a' else .4 for c in letters);us=[];offset=0
                    for size in lens:us.append(''.join(letters[offset:offset+size]));offset+=size
                    groups=[[i for i,u in enumerate(units) if u==v] for v in us]
                    if any(not gr for gr in groups):fail+=mass;continue
                    for cipher in product(*groups):generative[cipher]+=mass*math.prod(phi[c] for c in cipher)
            close(sum(generative.values())+fail,1)
            for cipher,p in generative.items():
                w=''.join(units[c] for c in cipher);n=len(w)
                q=pi(n)/den*math.prod(.6 if c=='a' else .4 for c in w)/math.comb(g,n-g)*math.prod(phi[c] for c in cipher)
                close(p,q)
            toy.append(dict(g=g,code_count=len(units),cipher_words=len(generative),failure=fail))
    sourcepairs=Counter(w[i:i+2] for w in words for i in range(len(w)-1))
    def seg(w,pairs):
        result=[]
        while w:
            u=w[:2] if w[:2] in pairs else w[:1];result.append(u);w=w[len(u):]
        return result
    def plainrows(rows):
        result=[]
        for row in rows:
            text=''.join(' ' if x==L else ALPH[x] for x in row);result.extend(text.split())
        return result
    def logword(w):
        seq=[L]*3+[ALPH.index(c) for c in w]+[L]
        return sum(float(lm[tuple(seq[i-3:i+1])]) for i in range(3,len(seq)))
    def score(rows,key,cc):
        totals=Counter()
        for i,u in enumerate(key):totals[u]+=cc[i]
        phi=[cc[i]/totals[u] for i,u in enumerate(key)]; langpart=lenpart=emitpart=mass=outchars=wn=hits=0
        for cipher,num in Counter(tuple(r['cipher']) for r in rows if min(r['cipher'])>=0).items():
            w=''.join(UNITS[key[c]] for c in cipher);n=len(w);g=len(cipher)
            langpart+=num*logword(w);lenpart+=num*terms[g,n];emitpart+=num*sum(math.log(phi[c]) for c in cipher)
            mass+=num*g;outchars+=num*n;wn+=num;hits+=num*(w in lang['word_counts'])
        return dict(objective=(langpart+lenpart+emitpart)/mass,word_language_per_input_unit=langpart/mass,length_normalization_per_input_unit=lenpart/mass,emission_per_input_unit=emitpart/mass,input_units=mass,output_letters=outchars,expansion_ratio=outchars/mass,complete_words=wn,word_hits=hits,source_word_fraction=hits/wn)
    def distance(a,b):
        mat=[[0]*(len(b)+1) for _ in range(len(a)+1)]
        for i in range(len(a)+1):mat[i][0]=i
        for j in range(len(b)+1):mat[0][j]=j
        for i,x in enumerate(a,1):
            for j,y in enumerate(b,1):mat[i][j]=min(mat[i-1][j]+1,mat[i][j-1]+1,mat[i-1][j-1]+int(x!=y))
        return mat[-1][-1]
    def accuracy(d,key):
        ww=[r for r in d['test'] if min(r['cipher'])>=0];wc=Counter((tuple(r['cipher']),r['plain'],tuple(r['units'])) for r in ww)
        ed=wr=tl=pl=ur=ul=un=0
        for (cipher,truth,units),num in wc.items():
            outputs=[UNITS[key[c]] for c in cipher];pred=''.join(outputs)
            ed+=num*distance(pred,truth);wr+=num*(pred==truth);tl+=num*len(truth);pl+=num*len(pred);un+=num*len(cipher)
            ur+=num*sum(a==b for a,b in zip(outputs,units));ul+=num*sum(len(a)==len(b) for a,b in zip(outputs,units))
        allletters=sum(len(r['plain']) for r in d['test'])
        return dict(edit_accuracy=1-ed/tl,edit_distance=ed,word_accuracy=wr/len(ww),correct_words=wr,complete_words=len(ww),total_words=len(d['test']),word_coverage=len(ww)/len(d['test']),letter_coverage=tl/allletters,true_letters=tl,total_true_letters=allletters,output_letters=pl,length_ratio=pl/tl,unit_string_accuracy=ur/un,unit_length_accuracy=ul/un,key_string_accuracy=sum(a==b for a,b in zip(key,d['true_key']))/len(key))
    calibration=read('results/variable-units-v1/calibration.json');rcount=0;objectives=0;regenerated=0;output=[]
    for trial,(name,source,npairs,extras) in enumerate([('letter23','cicero',0,0),('pair127','cicero',104,0),('homophone127','celsus',40,64)]):
        d=read(f'data/variable-units-v1/{name}.json');fit=read(f'results/variable-units-v1/{name}.fit.json');rec=calibration['controls'][trial]
        pairs=sorted(sourcepairs,key=lambda u:(-sourcepairs[u],u))[:npairs];assert pairs==d['pairs']
        us=list(ALPH)+pairs;sf=Counter(u for w in words for u in seg(w,set(pairs)));freq=np.array([sf[u] for u in us],dtype=float);freq/=sum(freq)
        ex=np.floor(extras*freq).astype(int)
        for i in sorted(range(len(us)),key=lambda i:(-(extras*freq[i]-ex[i]),us[i]))[:extras-sum(ex)]:ex[i]+=1
        rng=np.random.default_rng(20261025+trial);labels=list(map(int,rng.permutation(len(us)+extras)));groups={};off=0
        for i,u in enumerate(us):groups[u]=labels[off:off+1+ex[i]];off+=1+int(ex[i])
        assert groups==d['groups'];plain=read(f'data/latin/control-{source}.json');lookup={c:i for i,c in enumerate(d['cipher_labels'])}
        for split in ['train','test']:
            ww=plainrows(plain['plain_'+split]);assert len(ww)==len(d[split])
            for w,r in zip(ww,d[split]):
                units=seg(w,set(pairs));encoded=[]
                for u in units:
                    gr=groups[u];encoded.append(gr[0] if len(gr)==1 or rng.random()<.9 else int(rng.choice(gr[1:])))
                assert r==dict(plain=w,units=units,labels=encoded,cipher=[lookup.get(c,-1) for c in encoded]);regenerated+=1
        true=[next(UNITS.index(u) for u,gr in groups.items() if c in gr) for c in d['cipher_labels']];assert true==d['true_key']
        ct=Counter(c for r in d['train'] for c in r['cipher']);cc=[ct[c] for c in range(d['glyph_count'])];assert cc==fit['train_counts']
        assert fit['selected']==max(fit['restarts'],key=lambda r:r['objective'])
        for r,acc in zip(fit['restarts'],rec['restart_accuracy']):
            close(score(d['train'],r['key'],cc)['objective'],r['objective']);objectives+=1
            ac=accuracy(d,r['key']);verify(dict(restart=r['restart'],objective=r['objective'],**ac),acc);rcount+=1
        key=fit['selected']['key'];verify(accuracy(d,key),rec['accuracy'])
        for split in ['train','test']:
            verify(score(d[split],key,cc),rec[split+'_metrics']);verify(score(d[split],true,cc),rec['oracle_'+split+'_metrics'])
        a=rec['accuracy'];passed=a['edit_accuracy']>=.95 and a['word_accuracy']>=.85 and .95<=a['length_ratio']<=1.05;assert passed==rec['pass']
        output.append(dict(name=name,passed=passed,selected_objective=rec['train_metrics']['objective'],oracle_objective=rec['oracle_train_metrics']['objective']))
    assert (calibration['status']=='PASS')==all(r['passed'] for r in output)
    out=dict(status='PASS',meaning='Computational audit, not calibration success or decipherment.',frozen_hashes=frozen,lm_entries=lm.size,length_terms=termchecks,toy_generative_checks=toy,regenerated_cipher_words=regenerated,recomputed_restart_objectives=objectives,recomputed_restart_accuracies=rcount,controls=output)
    (ROOT/'results/audit-011-variable-units.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
