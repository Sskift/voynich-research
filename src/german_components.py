"""Experiment014: historical German on fixed EVA components, not a translation."""
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,re,subprocess
import numpy as np
from transcription_units import parse_source,split_words

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/german-components';OUT=ROOT/'results/german-components'
ALPH='abcdefghijklmnopqrstuvwxyzäæëößü';L=len(ALPH);A=L+1
def read(p):return json.loads((ROOT/p).read_text())
def dump(p,o):(ROOT/p).write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n')
def digest(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def snap(paths):return {p:digest(p) for p in paths}
def windows(rows,n,G):
    c=Counter()
    for row in rows:
        for i in range(len(row)-n+1):
            q=tuple(row[i:i+n])
            if max(q)<=G:c[q]+=1
    return c
def words_to_rows(words,group=None):
    groups=[words] if group is None else [words[i:i+group] for i in range(0,len(words),group)]
    return [[L]+[a for w in group for a in [*(ALPH.index(c) if c in ALPH else L+1 for c in w),L]] for group in groups]

def prepare():
    DATA.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    brant=read('data/language-screen/brant-reading.json');rows=[row for r in brant['records'] for row in words_to_rows(r['words'])];counts=[windows(rows,n,L) for n in range(1,5)];probs=[]
    for n,ct in enumerate(counts,1):
        arr=np.zeros((A,)*n)
        for q,v in ct.items():arr[q]=v
        prob=(arr+.1)/(arr.sum()+A*.1) if n==1 else (arr+10*probs[-1][None,...])/(arr.sum(axis=-1,keepdims=True)+10)
        assert np.allclose(prob.sum(axis=-1),1);probs.append(prob)
    np.log(probs[-1]).tofile(DATA/'conditional4.bin')
    vocab=Counter(''.join(w) for r in brant['records'] for w in r['words'] if all(c in ALPH for c in w));six=Counter()
    for row in rows:
        s=''.join(' ' if c==L else '□' if c>L else ALPH[c] for c in row)
        six.update(s[i:i+6] for i in range(len(s)-5) if '□' not in s[i:i+6])
    dump('data/german-components/language.json',dict(alphabet=ALPH,unigram=probs[0].tolist(),word_counts=dict(vocab),sixgram_counts=dict(six),source='data/language-screen/brant-reading.json',source_sha256=digest('data/language-screen/brant-reading.json'),lm_sha256=digest('data/german-components/conditional4.bin'),known_letter_units=sum(c<L for row in rows for c in row),opaque_unit_masks=sum(c>L for row in rows for c in row),context_mass=10,unigram_smoothing=.1))
    tepl=read('data/language-screen/tepl.json');ws=[w for r in tepl['records'] for w in r['words']];cut=int(.6*len(ws));assert all(c in ALPH for w in ws for c in w)
    dump('data/german-components/plain-control.json',dict(source='data/language-screen/tepl.json',source_sha256=digest('data/language-screen/tepl.json'),split_word_index=cut,total_word_tokens=len(ws),plain_train=words_to_rows(ws[:cut],10),plain_test=words_to_rows(ws[cut:],10)))
    selected=read('data/transcription-units/selection.json');native=[]
    for row in selected['rows']:
        z=row['sources']['ZL']
        if z['eligible']:native.append(dict(locus=row['locus'],split=row['split'],raw=z['raw'],words=z['words']['primary']))
    raw=parse_source('ZL');header=next(line for line in (ROOT/'data/raw/ZL3b-n.txt').read_text().splitlines() if line.startswith('<f1r>'))
    assert '$L=A' in header and '$H=1' in header and '$I=T' in header
    for loc,r in raw.items():
        if r['page']=='f1r' and r['eligible']:native.append(dict(locus=loc,split='score_reserve_f1r',raw=r['raw'],words=split_words(r['raw'])))
    ct=Counter(c for r in native if r['split']=='train' for w in r['words'] if w['known'] for c in w['atoms']);glyphs=sorted(ct,key=lambda g:(-ct[g],g));lookup={g:i for i,g in enumerate(glyphs)};G=len(glyphs)
    for row in native:
        row['encoded']=[G]
        for w in row['words']:
            row['encoded'].extend(lookup.get(a,G+1) for a in w['atoms']) if w['known'] else row['encoded'].append(G+1)
            row['encoded'].append(G)
    coverage={}
    for split in ['train','development','reserve','score_reserve_f1r']:
        rs=[r for r in native if r['split']==split]
        coverage[split]=dict(loci=len(rs),known_input_atoms=sum(len(w['atoms']) for r in rs for w in r['words'] if w['known']),unreadable_whole_words=sum(not w['known'] for r in rs for w in r['words']),training_seen_atoms=sum(c<G for r in rs for c in r['encoded']),unseen_known_atoms=sum(sum(a not in lookup for a in w['atoms']) for r in rs for w in r['words'] if w['known']))
    assert coverage['train']['known_input_atoms']==29865 and G==63
    dump('data/german-components/manuscript.json',dict(status='FROZEN_BEFORE_GERMAN_KEY_SEARCH',glyphs=glyphs,rows=native,coverage=coverage,reserve_note='f1r one short T page; global score new, source/illustrations previously seen. Old H/P are development.',input_hashes=snap(['data/transcription-units/selection-freeze.json','data/raw/ZL3b-n.txt'])))
    dump('data/german-components/config-v1.json',dict(version='v1',seed=20261014,restarts=16,cycles=3,steps_per_cycle=25000,temperature=[.02,.00005],emission_weight=1,alphabet=ALPH))
    print('PREPARED',G,coverage,flush=True)

def cipher_counts(rows,G):return np.array([sum(row.count(i) for row in rows) for i in range(G)],dtype=int)
def emission(counts,key):
    totals=np.bincount(key[:len(counts)],weights=counts,minlength=L)
    return [float(counts[i]/totals[key[i]]) for i in range(len(counts))]
def solve(name,rows,G,config,seed):
    counts=cipher_counts(rows,G);assert min(counts)>0
    ct=windows(rows,4,G);base=OUT/name
    base.with_suffix('.grams.txt').write_text(''.join(' '.join(map(str,q))+f' {n}\n' for q,n in sorted(ct.items())))
    base.with_suffix('.counts.txt').write_text(str(G)+'\n'+' '.join(map(str,counts))+'\n')
    lp=np.array(read('data/german-components/language.json')['unigram'][:L]);lp/=lp.sum();remaining=lp*counts.sum();initial=[0]*G
    for i in np.argsort(-counts,kind='stable'):
        k=int(remaining.argmax());initial[int(i)]=k;remaining[k]-=counts[i]
    initial.append(L);base.with_suffix('.initial.txt').write_text(' '.join(map(str,initial))+'\n')
    cmd=[str(ROOT/'bin/component_search'),str(DATA/'conditional4.bin'),str(base.with_suffix('.grams.txt')),str(base.with_suffix('.counts.txt')),str(base.with_suffix('.initial.txt')),str(seed),str(config['restarts']),str(config['steps_per_cycle']),str(config['cycles']),str(base.with_suffix('.keys.txt')),str(A)]
    with base.with_suffix('.log').open('w') as f:subprocess.run(cmd,check=True,stderr=f,stdout=f)
    restarts=[]
    for line in base.with_suffix('.keys.txt').read_text().splitlines():
        cells=line.split();restarts.append(dict(restart=int(cells[0]),objective=float(cells[1]),key=list(map(int,cells[2:]))))
    assert len(restarts)==config['restarts'];selected=max(restarts,key=lambda r:r['objective']);prob=emission(counts,selected['key'])
    out=dict(name=name,seed=seed,restarts=restarts,selected=selected,emission_probabilities=prob,train_counts=counts.tolist(),train_windows=sum(ct.values()),model_sha256=digest('data/german-components/conditional4.bin'))
    dump(f'results/german-components/{name}.frozen.json',out);print('FIT',name,selected['objective'],flush=True);return out

def metrics(rows,key,prob,G):
    lm=np.fromfile(DATA/'conditional4.bin',dtype=np.float64);ct=windows(rows,4,G);mass=sum(ct.values());total=0
    for q,n in ct.items():
        a,b,c,d=[key[x] for x in q];total+=n*lm[((a*A+b)*A+c)*A+d]
    ec=sum(n*np.log(prob[i]) for i,n in enumerate(cipher_counts(rows,G)) if n)
    model=read('data/german-components/language.json');six=model['sixgram_counts'];vocab=model['word_counts'];hits=sixn=wh=wn=0
    for row in rows:
        s=''.join('□' if c>G else ' ' if c==G else ALPH[key[c]] for c in row)
        for i in range(len(s)-5):
            q=s[i:i+6]
            if '□' not in q:sixn+=1;hits+=q in six
        for w in s.split():
            if len(w)>=3 and '□' not in w:wn+=1;wh+=w in vocab
    return dict(conditional4=float(total/mass),emission_per_window=float(ec/mass),combined=float((total+ec)/mass),windows=mass,sixgram_seen=hits/sixn if sixn else None,sixgram_hits=hits,sixgram_windows=sixn,word_seen=wh/wn if wn else None,word_hits=wh,complete_words_ge3=wn)

def generate_control(name,N,primary,mask,seed):
    plain=read('data/german-components/plain-control.json');rng=np.random.default_rng(seed);lm=np.array(read('data/german-components/language.json')['unigram'][:L]);lm/=lm.sum()
    extras=np.floor((N-L)*lm).astype(int)
    for i in np.argsort(-((N-L)*lm-extras),kind='stable')[:N-L-int(extras.sum())]:extras[i]+=1
    lengths=extras+1;labels=list(map(int,rng.permutation(N)));groups=[];truth={};off=0
    for letter,n in enumerate(lengths):
        group=labels[off:off+n];off+=n;groups.append(group)
        for c in group:truth[c]=letter
    raw={}
    for split in ['train','test']:
        rows=[]
        for row in plain['plain_'+split]:
            encoded=[]
            for letter in row:
                if letter==L:encoded.append(' ');continue
                if rng.random()<mask:encoded.append('?');continue
                group=groups[letter]
                if len(group)==1 or primary is not None and rng.random()<primary:c=group[0]
                else:c=int(rng.choice(group if primary is None else group[1:]))
                encoded.append(c)
            rows.append(encoded)
        raw[split]=rows
    counts=Counter(c for r in raw['train'] for c in r if isinstance(c,int));glyphs=sorted(counts,key=lambda c:(-counts[c],c));lookup={c:i for i,c in enumerate(glyphs)};G=len(glyphs)
    out=dict(name=name,N_requested=N,glyph_count=G,cipher_labels=glyphs,true_key=[truth[c] for c in glyphs]+[L],seed=seed,primary_probability=primary,mask_probability=mask,homophones_per_letter=lengths.tolist())
    for split in ['train','test']:
        out[split]=[[G if c==' ' else lookup.get(c,G+1) for c in row] for row in raw[split]];out['plain_'+split]=plain['plain_'+split]
        out[split+'_unseen_positions']=sum(isinstance(c,int) and c not in lookup for row in raw[split] for c in row)
    dump(f'data/german-components/{name}.json',out);return out

def accuracy(data,key):
    G=data['glyph_count'];right=total=wr=wn=all_letters=0
    for row,truth in zip(data['test'],data['plain_test']):
        word=[]
        for c,t in zip(row,truth):
            if t<L:all_letters+=1
            if c<G:total+=1;right+=key[c]==t
            if c==G:
                if word and all(a<G for a,b in word):wn+=1;wr+=all(key[a]==b for a,b in word)
                word=[]
            else:word.append((c,t))
    return dict(letter_accuracy=right/total,correct_letters=right,predicted_letters=total,total_letters=all_letters,predicted_fraction=total/all_letters,word_accuracy=wr/wn,correct_words=wr,complete_words=wn,key_accuracy=sum(a==b for a,b in zip(key[:-1],data['true_key'][:-1]))/G)

def controls():
    config=read('data/german-components/config-v1.json');records=[]
    for trial,(name,N,primary,mask) in enumerate([('v1-bijection32',32,None,0),('v1-uniform64',64,None,0),('v1-tail63',63,.97,.02)]):
        data=generate_control(name,N,primary,mask,config['seed']+trial);G=data['glyph_count'];fit=solve(name,data['train'],G,config,config['seed']+100+trial);key=fit['selected']['key']
        trueprob=emission(cipher_counts(data['train'],G),data['true_key'])
        rec=dict(name=name,glyph_count=G,accuracy=accuracy(data,key),selected_restart=fit['selected']['restart'],train_objective=fit['selected']['objective'],test_metrics=metrics(data['test'],key,fit['emission_probabilities'],G),true_key_metrics=metrics(data['test'],data['true_key'],trueprob,G),restart_accuracy=[dict(restart=r['restart'],objective=r['objective'],**accuracy(data,r['key'])) for r in fit['restarts']],data_sha256=digest(f'data/german-components/{name}.json'),fit_sha256=digest(f'results/german-components/{name}.frozen.json'))
        records.append(rec);dump('results/german-components/calibration-v1.json',dict(status='IN_PROGRESS',controls=records));print('CONTROL',name,rec['accuracy'],flush=True)
    passed=all(r['accuracy']['letter_accuracy']>=.95 for r in records)
    dump('results/german-components/calibration-v1.json',dict(status='PASS' if passed else 'FAIL',controls=records,criterion='>=95% correct predictable test letters in each selected condition; coverage separate; all controls are one Tepl work.'))
    if not passed:raise RuntimeError('Calibration failed: preserve it; do not search manuscript keys.')
    files=['notes/experiment-014-design.md','src/german_components.py','src/component_search.cpp','bin/component_search','src/transcription_units.py','data/german-components/manuscript.json','data/german-components/language.json','data/german-components/conditional4.bin','data/german-components/config-v1.json','results/german-components/calibration-v1.json','data/german-components/plain-control.json','data/language-screen/source-freeze.json']
    dump('results/german-components/solver-freeze.json',dict(status='FROZEN_BEFORE_MANUSCRIPT_KEY_SEARCH',hashes=snap(files)))

def attack():
    for p,h in read('results/german-components/solver-freeze.json')['hashes'].items():assert digest(p)==h,p
    config=read('data/german-components/config-v1.json');man=read('data/german-components/manuscript.json');G=len(man['glyphs']);splits={s:[r['encoded'] for r in man['rows'] if r['split']==s] for s in man['coverage']};target='score_reserve_f1r'
    fit=solve('voynich',splits['train'],G,config,config['seed']+1000);key=fit['selected']['key'];prob=fit['emission_probabilities']
    out=dict(status='NOT_A_TRANSLATION',complete=False,key=key,emission_probabilities=prob,glyphs=man['glyphs'],frozen_key_sha256=digest('results/german-components/voynich.frozen.json'),metrics={s:metrics(rows,key,prob,G) for s,rows in splits.items()},nulls=[])
    out['restart_reserve_metrics']=[dict(restart=r['restart'],objective=r['objective'],key=r['key'],metrics=metrics(splits[target],r['key'],emission(np.array(fit['train_counts']),r['key']),G)) for r in fit['restarts']]
    rng=np.random.default_rng(config['seed']+2000);scores=[];ct=windows(splits[target],4,G);qs=np.array(list(ct));n=np.array(list(ct.values()));model=np.fromfile(DATA/'conditional4.bin',dtype=np.float64)
    for _ in range(999):
        perm=rng.permutation(L);k=np.array([int(perm[v]) for v in key[:-1]]+[L]);d=k[qs];idx=((d[:,0]*A+d[:,1])*A+d[:,2])*A+d[:,3];scores.append(float(np.dot(model[idx],n)/n.sum()))
    out['random_label_scores']=scores;out['random_label_p']=(1+sum(s>=out['metrics'][target]['conditional4'] for s in scores))/1000
    decoded=[dict(locus=r['locus'],split=r['split'],source=r['raw'],output=''.join('□' if c>G else ' ' if c==G else ALPH[key[c]] for c in r['encoded'])) for r in man['rows']]
    dump('results/german-components/voynich-decoded.json',dict(status='NOT_A_TRANSLATION',glyph_to_letter=dict(zip(man['glyphs'],[ALPH[k] for k in key[:-1]])),rows=decoded));dump('results/german-components/attack.json',out)
    print('ACTUAL SCORED',out['metrics'],flush=True)
    for trial in range(9):
        rng=np.random.default_rng(config['seed']+3000+trial)
        def shuffled(rows):
            rr=[]
            for row in rows:
                idx=[i for i,c in enumerate(row) if c<G];new=list(row);vals=rng.permutation([row[i] for i in idx])
                for i,c in zip(idx,vals):new[i]=int(c)
                rr.append(new)
            return rr
        tr=shuffled(splits['train']);te=shuffled(splits[target]);nf=solve(f'null-{trial}',tr,G,config,config['seed']+4000+trial)
        out['nulls'].append(dict(trial=trial,shuffle_seed=config['seed']+3000+trial,train_objective=nf['selected']['objective'],reserve_metrics=metrics(te,nf['selected']['key'],nf['emission_probabilities'],G),fit_sha256=digest(f'results/german-components/null-{trial}.frozen.json')));dump('results/german-components/attack.json',out)
    out['null_conditional_p']=(1+sum(r['reserve_metrics']['conditional4']>=out['metrics'][target]['conditional4'] for r in out['nulls']))/10;out['complete']=True;dump('results/german-components/attack.json',out);print('ATTACK COMPLETE',out['metrics'][target],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','controls','attack']);a=p.parse_args()
    {'prepare':prepare,'controls':controls,'attack':attack}[a.stage]()
