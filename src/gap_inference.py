"""Experiment018: calibrated glyph-to-character/SPACE inference."""
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,math,subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/gap-inference';OUT=ROOT/'results/gap-inference'
LANGS=['la','de','it'];NAMES=['bijective','uniform127','tail127-position-mask']
def read(p):return json.loads((ROOT/p).read_text())
def dump(p,x):
    p=ROOT/p;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def snap(paths):return {p:sha(p) for p in paths}
def verify(p):
    for f,h in read(p)['hashes'].items():assert sha(f)==h,f
def verify_prepared():
    verify('data/gap-inference/input-freeze.json')
    manifest=read('data/gap-inference/prepared-manifest.json')
    assert manifest['status']=='PREPARED_NOT_FITTED'
    assert manifest['input_hashes']==read('data/gap-inference/input-freeze.json')['hashes']
    required={'data/gap-inference/config.json','data/gap-inference/manuscript.json','results/gap-inference/gate.json','results/gap-inference/solver-selfcheck.json'}
    for lang in LANGS:
        required.update({f'data/gap-inference/{lang}/language.json',f'data/gap-inference/{lang}/conditional.bin'})
        required.update(f'data/gap-inference/{lang}/controls/{name}.json' for name in NAMES)
    assert required <= set(manifest['hashes']),sorted(required-set(manifest['hashes']))
    verify('data/gap-inference/prepared-manifest.json')
def seedbase(lang):return 20261102+LANGS.index(lang)*10000
def segments(rows):
    result=[]
    for r in rows:
        seq=[]
        for c in r['cipher']:
            if c<0:
                if seq:result.append(seq);seq=[]
            else:seq.append(c)
        if seq:result.append(seq)
    return result
def grams(rows):
    cs=[Counter() for _ in range(4)]
    for row in segments(rows):
        for i in range(len(row)):
            n=min(4,i+1);cs[n-1][tuple(row[i-n+1:i+1])]+=1
    return cs
def counts(rows,G):return np.bincount([c for r in rows for c in r['cipher'] if c>=0],minlength=G)
def emission(cc,key,A):
    total=np.bincount(key,weights=cc,minlength=A)
    return [float(cc[c]/total[v]) for c,v in enumerate(key)]
def models(lang):
    info=read(f'data/gap-inference/{lang}/language.json');A=len(info['alphabet']);arr=np.fromfile(DATA/lang/'conditional.bin',dtype='<f8');off=0;ls=[]
    for n in range(1,5):ls.append(arr[off:off+A**n].reshape((A,)*n));off+=A**n
    assert off==len(arr)
    return info,ls

def solve(lang,name,rows,G,seed):
    config=read('data/gap-inference/config.json');info,lm=models(lang);A=len(info['alphabet']);cc=counts(rows,G);assert len(cc)==G and min(cc)>0
    gc=grams(rows);assert sum(sum(c.values()) for c in gc)==int(sum(cc));base=OUT/lang/name;base.parent.mkdir(parents=True,exist_ok=True)
    base.with_suffix('.grams.txt').write_text(''.join(f'{n} {v} '+' '.join(map(str,q))+'\n' for n,ct in enumerate(gc,1) for q,v in sorted(ct.items())))
    base.with_suffix('.counts.txt').write_text(str(G)+'\n'+'\n'.join(map(str,cc))+'\n')
    remain=np.array(info['unigram'])*cc.sum();initial=[0]*G
    for c in np.argsort(-cc,kind='stable'):
        v=int(remain.argmax());initial[int(c)]=v;remain[v]-=cc[c]
    rng=np.random.default_rng(seed+100000);initials=[]
    for r in range(config['restarts']):
        key=list(initial)
        if r:
            perm=rng.permutation(A);key=[int(perm[v]) for v in key]
            if r%2==0:key=[int(rng.integers(A)) if rng.random()<.3 else v for v in key]
        initials.append(key)
    base.with_suffix('.initial.txt').write_text(''.join(' '.join(map(str,key))+'\n' for key in initials))
    cmd=[str(ROOT/'bin/gap_search'),str(DATA/lang/'conditional.bin'),str(base.with_suffix('.grams.txt')),str(base.with_suffix('.counts.txt')),str(base.with_suffix('.initial.txt')),str(seed),str(config['restarts']),str(config['steps_per_cycle']),str(config['cycles']),str(config['polish']),str(base.with_suffix('.keys.txt')),str(A)]
    with base.with_suffix('.log').open('w') as f:subprocess.run(cmd,check=True,stdout=f,stderr=f)
    restarts=[]
    for line in base.with_suffix('.keys.txt').read_text().splitlines():
        cells=line.split();restarts.append(dict(restart=int(cells[0]),objective=float(cells[1]),key=list(map(int,cells[2:]))))
    assert len(restarts)==config['restarts'];selected=max(restarts,key=lambda r:r['objective']);phi=emission(cc,selected['key'],A)
    result=dict(language=lang,name=name,seed=seed,G=G,A=A,train_counts=cc.tolist(),predictions=int(cc.sum()),segments=len(segments(rows)),restarts=restarts,selected=selected,emission_probabilities=phi,model_sha256=sha(f'data/gap-inference/{lang}/conditional.bin'))
    dump(f'results/gap-inference/{lang}/{name}.fit.json',result);print('FIT',lang,name,selected['restart'],selected['objective'],flush=True);return result

def metrics(lang,rows,key,phi):
    info,lm=models(lang);alphabet=info['alphabet'];space=alphabet.index(' ');cs=grams(rows);language=0.;N=sum(sum(ct.values()) for ct in cs)
    for n,ct in enumerate(cs,1):
        for q,v in ct.items():language+=v*float(lm[n-1][tuple(key[c] for c in q)])
    cc=Counter(c for r in rows for c in r['cipher'] if c>=0);emiss=sum(n*math.log(phi[c]) for c,n in cc.items());decoded_space=sum(n for c,n in cc.items() if key[c]==space)
    words=[];adjacent_spaces=0
    for row in segments(rows):
        text=''.join(alphabet[key[c]] for c in row);ss=[i for i,c in enumerate(text) if c==' ']
        for a,b in zip(ss,ss[1:]):
            if b==a+1:adjacent_spaces+=1
            else:words.append(text[a+1:b])
    longer=[w for w in words if len(w)>=3];vocab=info['word_counts'];hits=sum(w in vocab for w in words);hits3=sum(w in vocab for w in longer)
    return dict(objective=(language+emiss)/N,language_per_character=language/N,emission_per_character=emiss/N,predicted_characters=N,segments=len(segments(rows)),predicted_spaces=decoded_space,space_fraction=decoded_space/N,adjacent_space_pairs=adjacent_spaces,complete_predicted_words=len(words),word_hits=hits,word_hit_fraction=hits/len(words) if words else None,complete_predicted_words_ge3=len(longer),word_hits_ge3=hits3,word_hit_fraction_ge3=hits3/len(longer) if longer else None)

def accuracy(data,key):
    space=data['alphabet'].index(' ');correct=predicted=allchars=letterright=letters=tp=fp=fn=wordright=words=allwords=0
    for row in data['test']:
        cipher=row['cipher'];plain=row['plain'];assert len(cipher)==len(plain);allchars+=len(plain)
        for c,t in zip(cipher,plain):
            if c<0:continue
            p=key[c];predicted+=1;correct+=p==t
            if t!=space:letters+=1;letterright+=p==t
            tp+=p==space and t==space;fp+=p==space and t!=space;fn+=p!=space and t==space
        boundaries=[i for i,c in enumerate(plain) if c==space]
        for a,b in zip(boundaries,boundaries[1:]):
            if b<=a+1:continue
            allwords+=1
            if all(c>=0 for c in cipher[a:b+1]):
                words+=1;wordright+=all(key[c]==t for c,t in zip(cipher[a:b+1],plain[a:b+1]))
    return dict(character_accuracy=correct/predicted,correct_characters=correct,predictable_characters=predicted,total_characters=allchars,character_coverage=predicted/allchars,letter_accuracy=letterright/letters,correct_letters=letterright,predictable_letters=letters,space_precision=tp/(tp+fp) if tp+fp else 0.,space_recall=tp/(tp+fn) if tp+fn else 0.,space_true_positives=tp,space_false_positives=fp,space_false_negatives=fn,word_accuracy=wordright/words if words else None,correct_words=wordright,complete_predictable_words=words,complete_true_words=allwords,complete_word_coverage=words/allwords if allwords else None,key_accuracy=sum(a==b for a,b in zip(key,data['true_key']))/len(key))

def calibrate():
    verify_prepared();assert read('results/gap-inference/solver-selfcheck.json')['status']=='PASS';gate=read('results/gap-inference/gate.json');summary={}
    for lang in LANGS:
        if not gate['primary_decisions'][lang]['admitted']:
            summary[lang]=dict(status='NOT_FIT_GATE_REJECTED');dump('results/gap-inference/calibration-index.json',summary);continue
        records=[]
        for i,name in enumerate(NAMES):
            data=read(f'data/gap-inference/{lang}/controls/{name}.json');fit=solve(lang,name,data['train'],data['G'],seedbase(lang)+100+i);key=fit['selected']['key'];phi=fit['emission_probabilities'];ac=accuracy(data,key)
            tr=metrics(lang,data['train'],key,phi);assert abs(tr['objective']-fit['selected']['objective'])<1e-8
            truthphi=emission(np.array(fit['train_counts']),data['true_key'],fit['A']);truthrows=data['test']
            record=dict(name=name,selected_restart=fit['selected']['restart'],accuracy=ac,passed=all(ac[k]>=.95 for k in ['character_accuracy','letter_accuracy','space_precision','space_recall']),train_metrics=tr,test_metrics=metrics(lang,data['test'],key,phi),true_test_metrics=metrics(lang,truthrows,data['true_key'],truthphi),restart_accuracy=[dict(restart=r['restart'],objective=r['objective'],**accuracy(data,r['key'])) for r in fit['restarts']],data_sha256=sha(f'data/gap-inference/{lang}/controls/{name}.json'))
            records.append(record);dump(f'results/gap-inference/{lang}/calibration.json',dict(status='IN_PROGRESS',controls=records));print('CONTROL',lang,name,ac,'PASS',record['passed'],flush=True)
        passed=all(r['passed'] for r in records);result=dict(status='PASS' if passed else 'FAIL',controls=records,criterion='Each selected condition: predictable character accuracy, letter accuracy, SPACE precision and recall all >=.95; coverage and complete words separate.')
        dump(f'results/gap-inference/{lang}/calibration.json',result);summary[lang]=dict(status=result['status']);dump('results/gap-inference/calibration-index.json',summary)
        if passed:
            paths=['notes/experiment-018-design.md','src/gap_inputs.py','src/gap_inference.py','src/gap_search.cpp','bin/gap_search','src/check_gap_solver.py','results/gap-inference/solver-selfcheck.json','data/gap-inference/input-freeze.json','data/gap-inference/config.json','data/gap-inference/manuscript.json','results/gap-inference/gate.json',f'data/gap-inference/{lang}/language.json',f'data/gap-inference/{lang}/conditional.bin',f'results/gap-inference/{lang}/calibration.json']
            paths += [f'data/gap-inference/{lang}/controls/{n}.json' for n in NAMES]
            paths += [f'results/gap-inference/{lang}/{n}.fit.json' for n in NAMES]
            paths += ['data/gap-inference/prepared-manifest.json']
            dump(f'results/gap-inference/{lang}/solver-freeze.json',dict(status='FROZEN_BEFORE_ACTUAL_GC_SEARCH',hashes=snap(paths)))
    print('CALIBRATION COMPLETE',summary,flush=True)

def random_scores(lang,rows,key,phi):
    info,lm=models(lang);A=len(info['alphabet']);parts=[];N=0
    for n,ct in enumerate(grams(rows),1):
        if ct:
            q=np.array(list(ct));weights=np.array(list(ct.values()));parts.append((n,np.array(key)[q],weights));N+=int(weights.sum())
    cc=Counter(c for r in rows for c in r['cipher'] if c>=0);constant=sum(n*math.log(phi[c]) for c,n in cc.items())/N
    rng=np.random.default_rng(seedbase(lang)+4000);scores=[]
    for _ in range(999):
        perm=rng.permutation(A);total=0.
        for n,q,weights in parts:
            d=perm[q];total+=float(np.dot(lm[n-1][tuple(d[:,i] for i in range(n))],weights))
        scores.append(total/N+constant)
    return scores

def attack():
    verify_prepared();index=read('results/gap-inference/calibration-index.json');man=read('data/gap-inference/manuscript.json');G=man['G'];rows={s:[r for r in man['rows'] if r['split']==s] for s in ['train','development','reserve']}
    for lang in LANGS:
        if index[lang]['status']!='PASS':continue
        verify(f'results/gap-inference/{lang}/solver-freeze.json');fit=solve(lang,'voynich',rows['train'],G,seedbase(lang)+1000);key=fit['selected']['key'];phi=fit['emission_probabilities'];info,_=models(lang);alphabet=info['alphabet']
        out=dict(status='NOT_A_TRANSLATION',complete=False,language=lang,key=key,emission_probabilities=phi,fit_sha256=sha(f'results/gap-inference/{lang}/voynich.fit.json'),metrics={s:metrics(lang,rr,key,phi) for s,rr in rows.items()},nulls=[],coverage=man['coverage'],scope='All79H/16H/16P are development. Visible word gaps omitted; unknown gaps/physical rows cut context. Predicted spaces are preserved.')
        assert abs(out['metrics']['train']['objective']-fit['selected']['objective'])<1e-8
        ranked=sorted(fit['restarts'],key=lambda r:r['objective'],reverse=True);first,second=ranked[:2]
        differing=[c for c,(a,b) in enumerate(zip(first['key'],second['key'])) if a!=b]
        out['key_diagnostics']=dict(unique_keys=len({tuple(r['key']) for r in ranked}),top_restarts=[first['restart'],second['restart']],top_objective_gap=first['objective']-second['objective'],different_glyph_types=len(differing),different_training_positions=sum(fit['train_counts'][c] for c in differing),training_positions=fit['predictions'],scope='Post-fit description; no restart reselection or global-optimum claim')
        out['restart_development_metrics']=[dict(restart=r['restart'],objective=r['objective'],metrics=metrics(lang,rows['reserve'],r['key'],emission(np.array(fit['train_counts']),r['key'],fit['A']))) for r in fit['restarts']]
        out['random_label_scores']=random_scores(lang,rows['reserve'],key,phi);out['random_label_p']=(1+sum(s>=out['metrics']['reserve']['objective'] for s in out['random_label_scores']))/1000
        decoded=[dict(locus=r['locus'],split=r['split'],source=r['raw'],output=''.join('□' if c<0 else alphabet[key[c]] for c in r['cipher']),output_characters=['□' if c<0 else alphabet[key[c]] for c in r['cipher']]) for r in man['rows']]
        dump(f'results/gap-inference/{lang}/voynich-decoded.json',dict(status='NOT_A_TRANSLATION',glyph_to_character=dict(zip(man['glyphs'],[alphabet[x] for x in key])),rows=decoded));dump(f'results/gap-inference/{lang}/attack.json',out);print('ACTUAL',lang,out['metrics'],flush=True)
        for trial in range(9):
            rng=np.random.default_rng(seedbase(lang)+2000+trial)
            def shuffled(rr):
                result=[]
                for r in rr:
                    cipher=r['cipher'];values=iter(map(int,rng.permutation([c for c in cipher if c>=0])));result.append(dict(cipher=[next(values) if c>=0 else -1 for c in cipher]))
                return result
            train=shuffled(rows['train']);test=shuffled(rows['reserve']);nf=solve(lang,f'null-{trial}',train,G,seedbase(lang)+3000+trial)
            rec=dict(trial=trial,shuffle_seed=seedbase(lang)+2000+trial,train_objective=nf['selected']['objective'],development_metrics=metrics(lang,test,nf['selected']['key'],nf['emission_probabilities']),fit_sha256=sha(f'results/gap-inference/{lang}/null-{trial}.fit.json'));out['nulls'].append(rec);dump(f'results/gap-inference/{lang}/attack.json',out);print('NULL',lang,trial,rec['development_metrics']['objective'],flush=True)
        out['shuffle_p']=(1+sum(r['development_metrics']['objective']>=out['metrics']['reserve']['objective'] for r in out['nulls']))/10;out['complete']=True;dump(f'results/gap-inference/{lang}/attack.json',out);print('ATTACK COMPLETE',lang,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['calibrate','attack']);a=p.parse_args();globals()[a.command]()
