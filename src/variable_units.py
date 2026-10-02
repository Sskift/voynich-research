"""Experiment015 v1: fixed one/two-letter code units; calibrated before manuscript use."""
from pathlib import Path
from collections import Counter
from itertools import product
import argparse, hashlib, json, math, subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/variable-units-v1'; OUT=ROOT/'results/variable-units-v1'
ALPH='abcdefghiklmnopqrstuxyz'; L=len(ALPH); A=L+1
UNITS=list(ALPH)+[a+b for a in ALPH for b in ALPH]; UID={s:i for i,s in enumerate(UNITS)}
MAXG=64; MAXN=128
def read(p): return json.loads((ROOT/p).read_text())
def dump(p,x):
    p=ROOT/p; p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def digest(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def snap(ps): return {p:digest(p) for p in ps}
def verify(p):
    for f,h in read(p)['hashes'].items(): assert digest(f)==h,f
def plain_words(rows):
    out=[]
    for row in rows:
        w=[]
        for c in row:
            if c==L:
                if w: out.append(''.join(w)); w=[]
            else:
                assert 0<=c<L; w.append(ALPH[c])
        assert not w
    return out
def segment(w,pairs):
    out=[]; i=0
    while i<len(w):
        u=w[i:i+2] if w[i:i+2] in pairs else w[i]
        out.append(u); i+=len(u)
    return out

def prepare():
    DATA.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
    assert not (DATA/'calibration-freeze.json').exists(), 'already frozen'
    paths=[f'data/latin/caesar-gallic-{n}.txt' for n in range(1,8)]
    words=[w for p in paths for w in (ROOT/p).read_text().split()]
    assert all(all(c in ALPH for c in w) for w in words)
    cs=[Counter() for _ in range(4)]
    for w in words:
        hist=[L,L,L]
        for v in [*(ALPH.index(c) for c in w),L]:
            for order in range(1,5): cs[order-1][tuple((hist+[v])[-order:])]+=1
            hist=(hist+[v])[-3:]
    ps=[]
    for order,ct in enumerate(cs,1):
        ar=np.zeros((A,)*order)
        for q,v in ct.items(): ar[q]=v
        p=(ar+.1)/(ar.sum()+A*.1) if order==1 else (ar+10*ps[-1][None,...])/(ar.sum(axis=-1,keepdims=True)+10)
        assert np.allclose(p.sum(axis=-1),1,atol=1e-14); ps.append(p)
    prob=ps[-1]; np.log(prob).tofile(DATA/'word-conditional4.bin')
    f=np.zeros((A,A,A)); f[L,L,L]=1.; z=[]
    for n in range(MAXN+1):
        z.append(float(np.sum(f*prob[:,:,:,L])))
        nf=np.zeros_like(f); nf[:,:,:L]=np.einsum('abc,abcd->bcd',f,prob[:,:,:,:L],optimize=True); f=nf
    lengths=Counter(map(len,words)); assert max(lengths)<=MAXN
    pi=[0.]+[(lengths[n]+100*.2*.8**(n-1))/(len(words)+100) for n in range(1,MAXN+1)]
    terms=np.full((MAXG+1,MAXN+1),np.nan)
    for g in range(1,MAXG+1):
        den=sum(pi[g:2*g+1])
        for n in range(g,2*g+1): terms[g,n]=math.log(pi[n]/den)-math.log(z[n])-math.log(math.comb(g,n-g))
    terms.tofile(DATA/'length-terms.bin')
    mono=Counter(''.join(words)); pairs=Counter(w[i:i+2] for w in words for i in range(len(w)-1))
    dump('data/variable-units-v1/language.json',dict(alphabet=ALPH,units=UNITS,source_hashes=snap(paths),word_tokens=len(words),word_counts=dict(Counter(words)),length_counts=dict(lengths),length_prior=pi,word_probability_mass_by_length=z,monograph_counts=dict(mono),digraph_counts=dict(pairs),context_mass=10,unigram_pseudocount=.1,prior_pseudowords=100,prior_geometric_stop=.2,max_input_word_units=MAXG,lm_sha256=digest('data/variable-units-v1/word-conditional4.bin'),terms_sha256=digest('data/variable-units-v1/length-terms.bin')))
    config=dict(seed=20261025,restarts=16,cycles=3,steps_per_cycle=25000,polish_passes=20,temperature=[.02,.00005],criterion=dict(edit_accuracy=.95,word_accuracy=.85,length_ratio_min=.95,length_ratio_max=1.05))
    dump('data/variable-units-v1/config.json',config)
    for trial,(name,source,npairs,extras) in enumerate([('letter23','cicero',0,0),('pair127','cicero',104,0),('homophone127','celsus',40,64)]):
        chosen=sorted(pairs,key=lambda u:(-pairs[u],u))[:npairs]; units=list(ALPH)+chosen; ps=set(chosen)
        freq=Counter(u for w in words for u in segment(w,ps)); n=len(units)
        fractions=np.array([freq[u] for u in units],dtype=float); fractions/=fractions.sum()
        extra=np.floor(extras*fractions).astype(int)
        ranking=sorted(range(n),key=lambda i:(-(extras*fractions[i]-extra[i]),units[i]))
        for i in ranking[:extras-int(extra.sum())]: extra[i]+=1
        rng=np.random.default_rng(config['seed']+trial); labels=list(map(int,rng.permutation(n+extras)))
        groups={}; truth={}; off=0
        for i,u in enumerate(units):
            groups[u]=labels[off:off+1+extra[i]]; off+=1+int(extra[i])
            for c in groups[u]: truth[c]=UID[u]
        plain=read(f'data/latin/control-{source}.json'); encoded={}
        for split in ['train','test']:
            ww=plain_words(plain['plain_'+split]); rows=[]
            for w in ww:
                us=segment(w,ps); cc=[]
                for u in us:
                    group=groups[u]
                    c=group[0] if len(group)==1 or rng.random()<.90 else int(rng.choice(group[1:]))
                    cc.append(c)
                rows.append(dict(plain=w,units=us,labels=cc))
            encoded[split]=rows
        ct=Counter(c for r in encoded['train'] for c in r['labels']); seen=sorted(ct,key=lambda c:(-ct[c],c)); lookup={c:i for i,c in enumerate(seen)}
        for rows in encoded.values():
            for r in rows: r['cipher']=[lookup.get(c,-1) for c in r['labels']]
        assert all(len(r['cipher'])<=MAXG for rows in encoded.values() for r in rows)
        control=dict(name=name,source=f'data/latin/control-{source}.json',source_sha256=digest(f'data/latin/control-{source}.json'),seed=config['seed']+trial,requested_codes=n+extras,glyph_count=len(seen),cipher_labels=seen,true_key=[truth[c] for c in seen],groups=groups,pairs=chosen,primary_probability=.90,**encoded)
        dump(f'data/variable-units-v1/{name}.json',control)
    print('PREPARED',len(words),'source words; 3 controls',flush=True)

def word_log(w,lm):
    a=b=c=L; total=0.
    for d in [*(ALPH.index(x) for x in w),L]:
        total+=lm[a,b,c,d]; a,b,c=b,c,d
    return float(total)

def normalization():
    lang=read('data/variable-units-v1/language.json'); lm=np.fromfile(DATA/'word-conditional4.bin').reshape((A,)*4)
    checks=[]
    for n in range(1,4):
        total=sum(math.exp(word_log(''.join(w),lm)) for w in product(ALPH,repeat=n)); expected=lang['word_probability_mass_by_length'][n]
        assert math.isclose(total,expected,rel_tol=1e-12,abs_tol=1e-14),(n,total,expected)
        checks.append(dict(test='Latin brute letter enumeration',length=n,enumerated_words=L**n,total=total,expected=expected))
    allunits=['a','b','aa','ab','ba','bb']; pi=lambda n:.2*.8**(n-1)
    books=[('complete',allunits,[1.]*6),('homophonic',allunits+['a'],[.3,1,1,1,1,1,.7]),('monograph',['a','b'],[1.,1.]),('digraph',allunits[2:],[1.]*4)]
    for name,units,phis in books:
        for g in range(1,4):
            den=sum(pi(n) for n in range(g,2*g+1)); total=0.
            for cs in product(range(len(units)),repeat=g):
                w=''.join(units[c] for c in cs); n=len(w)
                pword=math.prod(.6 if c=='a' else .4 for c in w)
                total+=pi(n)/den*pword/math.comb(g,n-g)*math.prod(phis[c] for c in cs)
            expected=pi(g)/den if name=='monograph' else pi(2*g)/den if name=='digraph' else 1.
            assert math.isclose(total,expected,rel_tol=1e-12,abs_tol=1e-14),(name,g,total,expected)
            checks.append(dict(test='exact toy codebook mass',codebook=name,cipher_word_length=g,total_valid=total,expected=expected,failure_mass=1-total))
    # Opposing length-only cases are genuine model probabilities, not averages.
    dump('results/variable-units-v1/normalization.json',dict(status='PASS',checks=checks,model_hash=digest('data/variable-units-v1/word-conditional4.bin'),limitation='Normalization is not a recovery or length-unbiasedness proof.'))
    print('NORMALIZATION PASS',len(checks),'checks',flush=True)

def freeze():
    assert read('results/variable-units-v1/normalization.json')['status']=='PASS'
    paths=['notes/experiment-015-design.md','notes/experiment-015-v1-protocol.md','src/variable_units.py','src/variable_unit_search.cpp','bin/variable_unit_search','results/variable-units-v1/normalization.json']
    paths += [str(p.relative_to(ROOT)) for p in sorted(DATA.iterdir()) if p.is_file()]
    assert 'data/variable-units-v1/calibration-freeze.json' not in paths
    dump('data/variable-units-v1/calibration-freeze.json',dict(status='FROZEN_BEFORE_CALIBRATION_SCORES',hashes=snap(paths)))
    print('FROZEN',len(paths),flush=True)

def counts(words,G): return np.bincount([c for w in words for c in w],minlength=G)
def emission(cc,key):
    sums=np.bincount(key,weights=cc,minlength=len(UNITS))
    return np.array([cc[c]/sums[u] for c,u in enumerate(key)])
def solve(name,words,G,seed):
    config=read('data/variable-units-v1/config.json'); lang=read('data/variable-units-v1/language.json')
    assert all(w and all(0<=c<G for c in w) and len(w)<=MAXG for w in words)
    ct=Counter(tuple(w) for w in words); cc=counts(words,G); assert min(cc)>0
    base=OUT/name
    base.with_suffix('.words.txt').write_text(''.join(f'{len(w)} {n} '+' '.join(map(str,w))+'\n' for w,n in sorted(ct.items())))
    base.with_suffix('.counts.txt').write_text(f'{G}\n'+' '.join(map(str,cc))+'\n')
    mono=np.array([lang['monograph_counts'].get(u,0) for u in UNITS[:L]],float); mono/=mono.sum()
    pairs=np.array([lang['digraph_counts'].get(u,0) for u in UNITS[L:]],float); pairs/=pairs.sum()
    rng=np.random.default_rng(seed+100000); initials=[]
    for r in range(config['restarts']):
        rho=0 if r==0 else [.25,.50,.75][(r-1)%3]
        remain=np.r_[(1-rho)*mono,rho*pairs]*cc.sum(); key=[0]*G
        for c in np.argsort(-cc,kind='stable'):
            u=int(np.argmax(remain)); key[int(c)]=u; remain[u]-=cc[c]
        if r>1:
            perm=rng.permutation(L); key=[UID[''.join(ALPH[perm[ALPH.index(ch)]] for ch in UNITS[u])] for u in key]
        initials.append(key)
    base.with_suffix('.initial.txt').write_text(''.join(' '.join(map(str,k))+'\n' for k in initials))
    cmd=[str(ROOT/'bin/variable_unit_search'),str(DATA/'word-conditional4.bin'),str(DATA/'length-terms.bin'),str(base.with_suffix('.words.txt')),str(base.with_suffix('.counts.txt')),str(base.with_suffix('.initial.txt')),str(seed),str(config['restarts']),str(config['steps_per_cycle']),str(config['cycles']),str(config['polish_passes']),str(base.with_suffix('.keys.txt'))]
    with base.with_suffix('.log').open('w') as f: subprocess.run(cmd,check=True,stdout=f,stderr=f)
    restarts=[]
    for line in base.with_suffix('.keys.txt').read_text().splitlines():
        xs=line.split(); restarts.append(dict(restart=int(xs[0]),objective=float(xs[1]),key=list(map(int,xs[2:]))))
    assert len(restarts)==config['restarts']
    selected=max(restarts,key=lambda r:r['objective']); phi=emission(cc,selected['key'])
    out=dict(name=name,seed=seed,restarts=restarts,selected=selected,train_counts=cc.tolist(),input_units=int(cc.sum()),word_tokens=len(words),word_types=len(ct),emission_probabilities=phi.tolist())
    dump(f'results/variable-units-v1/{name}.fit.json',out)
    print('FIT',name,selected['objective'],flush=True); return out

def edit(a,b):
    row=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        nxt=[i]
        for j,y in enumerate(b,1): nxt.append(min(nxt[-1]+1,row[j]+1,row[j-1]+(x!=y)))
        row=nxt
    return row[-1]
def accuracy(data,key):
    G=data['glyph_count']; allwords=len(data['test']); allletters=sum(len(r['plain']) for r in data['test'])
    selected=[r for r in data['test'] if all(c>=0 for c in r['cipher'])]
    grouped=Counter((tuple(r['cipher']),r['plain'],tuple(r['units'])) for r in selected)
    edits=exact=letters=predletters=rightunits=totalunits=rightlengths=0
    for (cw,truth,truthunits),n in grouped.items():
        pu=[UNITS[key[c]] for c in cw]; pred=''.join(pu)
        edits+=n*edit(pred,truth); exact+=n*(pred==truth); letters+=n*len(truth); predletters+=n*len(pred)
        rightunits+=n*sum(a==b for a,b in zip(pu,truthunits)); totalunits+=n*len(cw)
        rightlengths+=n*sum(len(a)==len(b) for a,b in zip(pu,truthunits))
    return dict(edit_accuracy=1-edits/letters,edit_distance=edits,word_accuracy=exact/len(selected),correct_words=exact,complete_words=len(selected),total_words=allwords,word_coverage=len(selected)/allwords,letter_coverage=letters/allletters,true_letters=letters,total_true_letters=allletters,output_letters=predletters,length_ratio=predletters/letters,unit_string_accuracy=rightunits/totalunits,unit_length_accuracy=rightlengths/totalunits,key_string_accuracy=sum(a==b for a,b in zip(key,data['true_key']))/G)
def metrics(words,key,phi):
    lm=np.fromfile(DATA/'word-conditional4.bin').reshape((A,)*4); terms=np.fromfile(DATA/'length-terms.bin').reshape(MAXG+1,MAXN+1)
    vocab=read('data/variable-units-v1/language.json')['word_counts']; ct=Counter(tuple(w) for w in words if w and min(w)>=0)
    language=length=emiss=units=outletters=tokens=hits=0
    for w,n in ct.items():
        decoded=''.join(UNITS[key[c]] for c in w); m=len(decoded)
        language+=n*word_log(decoded,lm); length+=n*terms[len(w),m]
        emiss+=n*sum(math.log(phi[c]) for c in w)
        units+=n*len(w); outletters+=n*m; tokens+=n; hits+=n*(decoded in vocab)
    return dict(objective=float((language+length+emiss)/units),word_language_per_input_unit=language/units,length_normalization_per_input_unit=float(length/units),emission_per_input_unit=emiss/units,input_units=units,output_letters=outletters,expansion_ratio=outletters/units,complete_words=tokens,word_hits=hits,source_word_fraction=hits/tokens)

def controls():
    verify('data/variable-units-v1/calibration-freeze.json'); records=[]
    for trial,name in enumerate(['letter23','pair127','homophone127']):
        data=read(f'data/variable-units-v1/{name}.json'); G=data['glyph_count']; words=[r['cipher'] for r in data['train']]
        fit=solve(name,words,G,20261125+trial); key=fit['selected']['key']; phi=fit['emission_probabilities']
        test=[r['cipher'] for r in data['test']]; oracle_phi=emission(np.array(fit['train_counts']),data['true_key'])
        tr=metrics(words,key,phi); assert abs(tr['objective']-fit['selected']['objective'])<1e-8
        rec=dict(name=name,selected_restart=fit['selected']['restart'],accuracy=accuracy(data,key),train_metrics=tr,test_metrics=metrics(test,key,phi),oracle_train_metrics=metrics(words,data['true_key'],oracle_phi),oracle_test_metrics=metrics(test,data['true_key'],oracle_phi),restart_accuracy=[dict(restart=r['restart'],objective=r['objective'],**accuracy(data,r['key'])) for r in fit['restarts']])
        ac=rec['accuracy']; rec['pass']=ac['edit_accuracy']>=.95 and ac['word_accuracy']>=.85 and .95<=ac['length_ratio']<=1.05
        records.append(rec); dump('results/variable-units-v1/calibration.json',dict(status='IN_PROGRESS',controls=records))
        print('CONTROL',name,rec['accuracy'],'PASS',rec['pass'],flush=True)
    passed=all(r['pass'] for r in records)
    dump('results/variable-units-v1/calibration.json',dict(status='PASS' if passed else 'FAIL',controls=records,criterion=read('data/variable-units-v1/config.json')['criterion'],note='Cicero is reused for first two conditions; all known works previously inspected in earlier experiments.'))
    print('CALIBRATION', 'PASS' if passed else 'FAIL',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('command',choices=['prepare','normalization','freeze','controls']); a=p.parse_args(); globals()[a.command]()
