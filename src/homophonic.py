"""Experiment011. One many-to-one key, fitted cipher emissions, reserved scores."""
from pathlib import Path
from collections import Counter
import argparse, json, math, re, subprocess
import numpy as np
from substitution import ROOT,ALPH,read,dump,digest,snap
from symbol_sequences import raw_loci,clean

DATA=ROOT/'data/homophonic';OUT=ROOT/'results/homophonic'
RESERVE=['f88r','f88v','f89r1','f89r2','f89v2','f89v1','f99r','f99v','f100r','f100v','f101r','f101v','f102r1','f102r2','f102v2','f102v1']

def windows(rows,n,G):
    counts=Counter()
    for row in rows:
        for i in range(len(row)-n+1):
            q=tuple(row[i:i+n])
            if max(q)<=G:counts[q]+=1
    return counts

def encode(rows,glyphs):
    lookup={g:i for i,g in enumerate(glyphs)};G=len(glyphs)
    return [[G if a==' ' else lookup.get(a,G+1) for a in r['atoms']] for r in rows]

def make_model():
    ns=[Counter() for _ in range(4)];words=Counter();six=Counter()
    files=read('data/latin-model.json')['training_files']
    for file,h in files.items():
        assert digest(file)==h
        for text in (ROOT/file).read_text().splitlines():
            words.update(text.split());s=' '+text+' '
            row=[23 if c==' ' else ALPH.index(c) for c in s]
            for n in range(1,5):ns[n-1].update(tuple(row[i:i+n]) for i in range(len(row)-n+1))
            six.update(s[i:i+6] for i in range(len(s)-5))
    probs=[];counts=[]
    for n,ct in enumerate(ns,1):
        arr=np.zeros((24,)*n)
        for q,v in ct.items():arr[q]=v
        counts.append(arr)
        if n==1:p=(arr+.1)/(arr.sum()+2.4)
        else:p=(arr+10*probs[-1][None,...])/(arr.sum(axis=-1,keepdims=True)+10)
        assert np.allclose(p.sum(axis=-1),1)
        probs.append(p)
    np.log(probs[-1]).tofile(DATA/'conditional4.bin')
    dump('data/homophonic/language.json',dict(alphabet=ALPH,unigram=probs[0].tolist(),training_files=files,unigram_smoothing=.1,context_backoff_mass=10,conditional4_sha256=digest('data/homophonic/conditional4.bin'),word_counts=dict(words),sixgram_counts=dict(six)))

def prepare():
    DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    prev=read('data/substitution-manuscript.json');counts=prev['training_glyph_counts']
    glyphs=sorted(counts,key=lambda x:(-counts[x],x));rows=[]
    for r in prev['rows']:
        rows.append(dict(page=r['page'],locus=r['locus'],raw=r['raw'],atoms=r['atoms'],split='train' if r['split']=='train' else 'development'))
    metadata={}
    for line in (ROOT/'data/raw/GC2a-n.txt').read_text().splitlines():
        m=re.fullmatch(r'<(f[^.>]+)>\s*(.*)',line)
        if m:metadata[m[1]]=dict(re.findall(r'\$([A-Z])=(\S)',m[2]))
    for page in RESERVE:
        assert all(metadata[page].get(k)==v for k,v in {'L':'A','H':'1','I':'P'}.items())
        items=list(raw_loci('GC',page).items())
        for j,(num,r) in enumerate(items):
            if r['code'] not in {'@P0','+P0','*P0'}:continue
            if j+1<len(items) and items[j+1][1]['code'][0] in '=~' and items[j+1][1]['code'][1]=='P':continue
            s=clean(r['raw']).replace('{','').replace('}','');atoms=[]
            for word in re.split('[.,]',s):
                if not word:continue
                tokens=['?'] if '?' in word or '[' in word else re.findall(r'@\d{3};|.',word)
                if atoms:atoms.append(' ')
                atoms.extend(tokens)
            rows.append(dict(page=page,locus=f'{page}.{num}',raw=r['raw'],atoms=[' ']+atoms+[' '],split='reserve'))
    encoded=encode(rows,glyphs)
    for r,e in zip(rows,encoded):r['encoded']=e
    G=len(glyphs);coverage={}
    for split in ['train','development','reserve']:
        rs=[r for r in rows if r['split']==split];rr=[r['encoded'] for r in rs]
        coverage[split]=dict(pages=sorted({r['page'] for r in rs}),loci=len(rs),known=sum(c<G for r in rr for c in r),masked=sum(c==G+1 for r in rr for c in r),windows=sum(windows(rr,4,G).values()))
    dump('data/homophonic/manuscript.json',dict(status='FROZEN_BEFORE_HOMOPHONIC_SCORES',glyphs=glyphs,space=G,unknown=G+1,rows=rows,coverage=coverage,inputs=snap(['notes/experiment-011-design.md','data/raw/GC2a-n.txt','data/substitution-manuscript.json'])))
    make_model()
    dump('data/homophonic/config-v1.json',dict(version='v1',seed=20261011,restarts=16,steps_per_cycle=25000,cycles=3,temperature=[.02,.00005],emission_weight=1))
    print('PREPARED',G,coverage,flush=True)

def cipher_counts(rows,G):return np.array([sum(row.count(i) for row in rows) for i in range(G)],dtype=int)

def emission(counts,key):
    totals=np.bincount(key[:len(counts)],weights=counts,minlength=23)
    p=np.array([counts[i]/totals[key[i]] for i in range(len(counts))])
    assert np.isfinite(p).all() and (p>0).all()
    return p

def solve(name,rows,G,config,seed):
    counts=cipher_counts(rows,G);assert min(counts)>0
    ct=windows(rows,4,G);base=OUT/name
    base.with_suffix('.grams.txt').write_text(''.join(' '.join(map(str,q))+f' {n}\n' for q,n in sorted(ct.items())))
    base.with_suffix('.counts.txt').write_text(str(G)+'\n'+' '.join(map(str,counts))+'\n')
    lp=np.array(read('data/homophonic/language.json')['unigram'][:23]);lp/=lp.sum();remaining=lp*counts.sum();initial=[0]*G
    for i in np.argsort(-counts,kind='stable'):
        k=int(remaining.argmax());initial[int(i)]=k;remaining[k]-=counts[i]
    initial.append(23);base.with_suffix('.initial.txt').write_text(' '.join(map(str,initial))+'\n')
    command=[str(ROOT/'bin/homophonic_search'),str(DATA/'conditional4.bin'),str(base.with_suffix('.grams.txt')),str(base.with_suffix('.counts.txt')),str(base.with_suffix('.initial.txt')),str(seed),str(config['restarts']),str(config['steps_per_cycle']),str(config['cycles']),str(base.with_suffix('.keys.txt'))]
    with base.with_suffix('.log').open('w') as f:subprocess.run(command,stderr=f,check=True)
    results=[]
    for line in base.with_suffix('.keys.txt').read_text().splitlines():
        c=line.split();results.append(dict(restart=int(c[0]),objective=float(c[1]),key=list(map(int,c[2:]))))
    best=max(results,key=lambda x:x['objective']);prob=emission(counts,best['key'])
    out=dict(status='KEY_AND_EMISSIONS_FROZEN_BEFORE_TEST_SCORE',name=name,config=config,seed=seed,glyph_count=G,train_counts=counts.tolist(),windows=sum(ct.values()),restarts=results,selected=best,emission_probabilities=prob.tolist())
    dump(f'results/homophonic/{name}.frozen.json',out);print('SOLVED',name,best['objective'],flush=True)
    return out

def metrics(rows,key,prob,G):
    model=np.fromfile(DATA/'conditional4.bin',dtype=np.float64);ct=windows(rows,4,G);mass=sum(ct.values());lm=0
    for q,n in ct.items():
        d=[key[c] for c in q];lm+=n*model[((d[0]*24+d[1])*24+d[2])*24+d[3]]
    ec=sum(math.log(prob[c]) for row in rows for c in row if c<G)
    language=read('data/homophonic/language.json');six=language['sixgram_counts'];vocab=language['word_counts']
    hits=total=wordhits=wordtotal=0
    for row in rows:
        s=''.join('□' if c==G+1 else ' ' if c==G else ALPH[key[c]] for c in row)
        for i in range(len(s)-5):
            q=s[i:i+6]
            if '□' not in q:total+=1;hits+=q in six
        for w in s.split():
            if len(w)>=3 and '□' not in w:wordtotal+=1;wordhits+=w in vocab
    return dict(conditional4=lm/mass,emission_per_window=ec/mass,combined=(lm+ec)/mass,windows=mass,sixgram_seen=hits/total if total else None,sixgram_hits=hits,sixgram_windows=total,word_seen=wordhits/wordtotal if wordtotal else None,word_hits=wordhits,complete_words_ge3=wordtotal)

def generate_control(name,plain,N,primary,mask,seed):
    rng=np.random.default_rng(seed);lm=np.array(read('data/homophonic/language.json')['unigram'][:23]);lm/=lm.sum()
    extras=np.floor((N-23)*lm).astype(int)
    for i in np.argsort(-((N-23)*lm-extras),kind='stable')[:N-23-int(extras.sum())]:extras[i]+=1
    lengths=extras+1;labels=list(map(int,rng.permutation(N)));groups=[];truth={};off=0
    for letter,n in enumerate(lengths):
        group=labels[off:off+n];off+=n;groups.append(group)
        for c in group:truth[c]=letter
    raw={}
    for split in ['train','test']:
        rows=[]
        for row in plain['plain_'+split]:
            encoded=[]
            for l in row:
                if l==23:encoded.append(' ');continue
                if rng.random()<mask:encoded.append('?');continue
                group=groups[l]
                if len(group)==1 or primary is not None and rng.random()<primary:c=group[0]
                else:c=int(rng.choice(group if primary is None else group[1:]))
                encoded.append(c)
            rows.append(encoded)
        raw[split]=rows
    counts=Counter(c for r in raw['train'] for c in r if isinstance(c,int));glyphs=sorted(counts,key=lambda c:(-counts[c],c));lookup={c:i for i,c in enumerate(glyphs)};G=len(glyphs)
    out=dict(name=name,N_requested=N,glyph_count=G,cipher_labels=glyphs,true_key=[truth[c] for c in glyphs]+[23],seed=seed,primary_probability=primary,mask_probability=mask,homophones_per_letter=lengths.tolist())
    for split in ['train','test']:
        out[split]=[[G if c==' ' else lookup.get(c,G+1) for c in row] for row in raw[split]]
        out['plain_'+split]=plain['plain_'+split]
        out[split+'_unseen_positions']=sum(isinstance(c,int) and c not in lookup for row in raw[split] for c in row)
    dump(f'data/homophonic/{name}.json',out);return out

def accuracy(data,key):
    G=data['glyph_count'];right=total=wordright=wordtotal=0;letterall=0
    for row,truth in zip(data['test'],data['plain_test']):
        word=[]
        for c,t in zip(row,truth):
            if t<23:letterall+=1
            if c<G:total+=1;right+=key[c]==t
            if c==G:
                if word and all(a<G for a,b in word):wordtotal+=1;wordright+=all(key[a]==b for a,b in word)
                word=[]
            else:word.append((c,t))
    return dict(letter_accuracy=right/total,correct_letters=right,predicted_letters=total,total_letters=letterall,predicted_fraction=total/letterall,word_accuracy=wordright/wordtotal,correct_words=wordright,complete_words=wordtotal,key_accuracy=sum(a==b for a,b in zip(key[:-1],data['true_key'][:-1]))/G)

def controls(version):
    config=read(f'data/homophonic/config-{version}.json');records=[]
    jobs=[('uniform60','control-cicero',60,None,0),('medical127','control-celsus',127,.90,0),('tail127','control-cicero',127,.97,.02)]
    for trial,(label,source,N,primary,mask) in enumerate(jobs):
        name=f'{version}-{label}';data=generate_control(name,read(f'data/latin/{source}.json'),N,primary,mask,config['seed']+trial)
        G=data['glyph_count'];fit=solve(name,data['train'],G,config,config['seed']+100+trial);key=fit['selected']['key']
        trueprob=emission(cipher_counts(data['train'],G),data['true_key'])
        rec=dict(name=name,glyph_count=G,accuracy=accuracy(data,key),selected_restart=fit['selected']['restart'],train_objective=fit['selected']['objective'],test_metrics=metrics(data['test'],key,fit['emission_probabilities'],G),true_key_metrics=metrics(data['test'],data['true_key'],trueprob,G),restart_accuracy=[dict(restart=r['restart'],objective=r['objective'],**accuracy(data,r['key'])) for r in fit['restarts']],data_sha256=digest(f'data/homophonic/{name}.json'),fit_sha256=digest(f'results/homophonic/{name}.frozen.json'))
        records.append(rec);dump(f'results/homophonic/calibration-{version}.json',dict(status='IN_PROGRESS',controls=records));print('CONTROL',name,rec['accuracy'],flush=True)
    passed=all(r['accuracy']['letter_accuracy']>=.95 for r in records)
    dump(f'results/homophonic/calibration-{version}.json',dict(status='PASS' if passed else 'FAIL',controls=records,criterion='>=95% correct predicted test letters for each condition; coverage separately reported'))
    if not passed:raise RuntimeError('Calibration failed. Preserve attempt; tune on controls only.')
    files=['notes/experiment-011-design.md','src/homophonic.py','src/homophonic_search.cpp','src/substitution.py','src/symbol_sequences.py','data/homophonic/manuscript.json','data/homophonic/language.json','data/homophonic/conditional4.bin',f'data/homophonic/config-{version}.json',f'results/homophonic/calibration-{version}.json']
    dump('results/homophonic/solver-freeze.json',dict(status='FROZEN_BEFORE_MANUSCRIPT_SEARCH',version=version,hashes=snap(files)))

def attack():
    freeze=read('results/homophonic/solver-freeze.json')
    for p,h in freeze['hashes'].items():assert digest(p)==h,p
    config=read(f"data/homophonic/config-{freeze['version']}.json");man=read('data/homophonic/manuscript.json');G=len(man['glyphs'])
    splits={s:[r['encoded'] for r in man['rows'] if r['split']==s] for s in ['train','development','reserve']}
    fit=solve('voynich',splits['train'],G,config,config['seed']+1000);key=fit['selected']['key'];prob=fit['emission_probabilities']
    out=dict(status='NOT_A_TRANSLATION',key=key,emission_probabilities=prob,glyphs=man['glyphs'],frozen_key_sha256=digest('results/homophonic/voynich.frozen.json'),metrics={s:metrics(rows,key,prob,G) for s,rows in splits.items()},nulls=[])
    # No later key selection uses these diagnostics.
    out['restart_reserve_metrics']=[dict(restart=r['restart'],objective=r['objective'],key=r['key'],metrics=metrics(splits['reserve'],r['key'],emission(np.array(fit['train_counts']),r['key']),G)) for r in fit['restarts']]
    rng=np.random.default_rng(config['seed']+2000);randomscores=[]
    ct=windows(splits['reserve'],4,G);qs=np.array(list(ct));n=np.array(list(ct.values()));model=np.fromfile(DATA/'conditional4.bin',dtype=np.float64)
    for i in range(999):
        perm=rng.permutation(23);k=np.array([int(perm[v]) for v in key[:-1]]+[23]);d=k[qs];idx=((d[:,0]*24+d[:,1])*24+d[:,2])*24+d[:,3]
        randomscores.append(float(np.dot(model[idx],n)/n.sum()))
    out['random_label_scores']=randomscores;out['random_label_p']=(1+sum(x>=out['metrics']['reserve']['conditional4'] for x in randomscores))/1000
    old=read('results/substitution/voynich.frozen-key.json')['selected']['key'];common=read('data/substitution-manuscript.json')['glyphs'];lookup={g:i for i,g in enumerate(common)}
    shared=[[23 if a==' ' else lookup.get(a,24) for a in r['atoms']] for r in man['rows'] if r['split']=='reserve']
    newcommon=[key[man['glyphs'].index(g)] for g in common]+[23]
    out['common23_candidate_comparison']={label:metrics(shared,k,[1]*23,23) for label,k in [('old_frozen_bijective_key',old),('new_frozen_homophonic_key',newcommon)]}
    decoded=[dict(locus=r['locus'],split=r['split'],source=r['raw'],output=''.join('□' if c==G+1 else ' ' if c==G else ALPH[key[c]] for c in r['encoded'])) for r in man['rows']]
    dump('results/homophonic/voynich-decoded.json',dict(status='NOT_A_TRANSLATION',glyph_to_letter=dict(zip(man['glyphs'],[ALPH[k] for k in key[:-1]])),rows=decoded))
    dump('results/homophonic/attack.json',out)
    for trial in range(9):
        nrng=np.random.default_rng(config['seed']+3000+trial)
        def shuffled(rows):
            result=[]
            for row in rows:
                idx=[i for i,c in enumerate(row) if c<G];rr=list(row);values=nrng.permutation([row[i] for i in idx])
                for i,c in zip(idx,values):rr[i]=int(c)
                result.append(rr)
            return result
        tr=shuffled(splits['train']);te=shuffled(splits['reserve']);f=solve(f'null-{trial}',tr,G,config,config['seed']+4000+trial)
        out['nulls'].append(dict(trial=trial,shuffle_seed=config['seed']+3000+trial,train_objective=f['selected']['objective'],reserve_metrics=metrics(te,f['selected']['key'],f['emission_probabilities'],G),fit_sha256=digest(f'results/homophonic/null-{trial}.frozen.json')))
        dump('results/homophonic/attack.json',out)
    out['null_conditional_p']=(1+sum(r['reserve_metrics']['conditional4']>=out['metrics']['reserve']['conditional4'] for r in out['nulls']))/10
    out['complete']=True;dump('results/homophonic/attack.json',out);print('ATTACK COMPLETE',out['metrics'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','controls','attack']);p.add_argument('--version',default='v1');a=p.parse_args()
    if a.stage=='prepare':prepare()
    elif a.stage=='controls':controls(a.version)
    else:attack()
