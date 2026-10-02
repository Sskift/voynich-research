"""Apply frozen015 solver to old development pages; no new reserve consumed."""
from collections import Counter
import argparse,math
import numpy as np
import variable_units as v

def prepare():
    v.verify('data/variable-units-v1/calibration-freeze.json')
    assert v.read('results/variable-units-v1/calibration.json')['status']=='PASS'
    assert v.read('results/audit-011-variable-units.json')['status']=='PASS'
    old=v.read('data/homophonic/manuscript.json');G=len(old['glyphs']);rows=[]
    for r in old['rows']:
        words=[];word=[]
        for c in r['encoded']:
            if c==G:
                if word:words.append(word);word=[]
            else:word.append(c if c<G else -1)
        assert not word and all(len(w)<=v.MAXG for w in words)
        rows.append(dict(locus=r['locus'],page=r['page'],split=r['split'],raw=r['raw'],atoms=r['atoms'],words=words))
    coverage={}
    for split in ['train','development','reserve']:
        rs=[r for r in rows if r['split']==split];ws=[w for r in rs for w in r['words']];complete=[w for w in ws if min(w)>=0]
        coverage[split]=dict(loci=len(rs),word_slots=len(ws),complete_words=len(complete),incomplete_words=len(ws)-len(complete),complete_input_units=sum(map(len,complete)),all_encoded_positions=sum(map(len,ws)),known_units_outside_complete_words=sum(c>=0 for w in ws if min(w)<0 for c in w),unknown_positions=sum(c<0 for w in ws for c in w))
    assert coverage['train']['complete_input_units']==23721
    v.dump('data/variable-units-v1/manuscript.json',dict(status='FROZEN_BEFORE_015_MANUSCRIPT_SCORES',glyphs=old['glyphs'],rows=rows,coverage=coverage,source_sha256=v.digest('data/homophonic/manuscript.json'),scope_note='79H training,16H and16P development. reserve is a legacy split label, NOT fresh confirmation. f49v/f58r/f58v remain unscored.'))
    paths=['src/variable_units_attack.py','data/variable-units-v1/manuscript.json','data/variable-units-v1/calibration-freeze.json','results/variable-units-v1/calibration.json','src/audit_variable_units.py','results/audit-011-variable-units.json','data/homophonic/manuscript.json']
    paths += [f'results/variable-units-v1/{n}.fit.json' for n in ['letter23','pair127','homophone127']]
    v.dump('results/variable-units-v1/solver-freeze.json',dict(status='FROZEN_BEFORE_MANUSCRIPT_KEY_SEARCH',hashes=v.snap(paths)))
    print('ACTUAL INPUT FROZEN',coverage,flush=True)

def randomlabels(words,key,phi):
    lm=np.fromfile(v.DATA/'word-conditional4.bin').reshape((v.A,)*4)
    ct=Counter()
    for w in words:
        if min(w)<0:continue
        plain=''.join(v.UNITS[key[c]] for c in w);seq=[v.L]*3+[v.ALPH.index(c) for c in plain]+[v.L]
        ct.update(tuple(seq[i-3:i+1]) for i in range(3,len(seq)))
    qs=np.array(list(ct));weights=np.array(list(ct.values()));base=v.metrics(words,key,phi)
    constant=base['length_normalization_per_input_unit']+base['emission_per_input_unit']
    rng=np.random.default_rng(20263025);scores=[]
    for _ in range(999):
        perm=np.r_[rng.permutation(v.L),v.L];d=perm[qs]
        scores.append(float(np.dot(lm[d[:,0],d[:,1],d[:,2],d[:,3]],weights)/base['input_units']+constant))
    return scores

def attack():
    v.verify('data/variable-units-v1/calibration-freeze.json');v.verify('results/variable-units-v1/solver-freeze.json')
    man=v.read('data/variable-units-v1/manuscript.json');G=len(man['glyphs'])
    rs={s:[r for r in man['rows'] if r['split']==s] for s in man['coverage']}
    ws={s:[w for r in rows for w in r['words']] for s,rows in rs.items()}
    fit=v.solve('voynich',ws['train'],G,20262025);key=fit['selected']['key'];phi=fit['emission_probabilities']
    decoded=[dict(locus=r['locus'],split=r['split'],source=r['raw'],words=[''.join(v.UNITS[key[c]] if c>=0 else '□' for c in w) for w in r['words']],whole_word_score_mask=[min(w)<0 for w in r['words']]) for r in man['rows']]
    v.dump('results/variable-units-v1/voynich-decoded.json',dict(status='NOT_A_TRANSLATION',glyph_to_output=dict(zip(man['glyphs'],[v.UNITS[u] for u in key])),rows=decoded))
    out=dict(status='NOT_A_TRANSLATION',complete=False,key=key,emission_probabilities=phi,fit_sha256=v.digest('results/variable-units-v1/voynich.fit.json'),metrics={s:v.metrics(words,key,phi) for s,words in ws.items()},coverage=man['coverage'],nulls=[],scope_note=man['scope_note'])
    assert abs(out['metrics']['train']['objective']-fit['selected']['objective'])<1e-8
    out['restart_development_metrics']=[dict(restart=r['restart'],objective=r['objective'],metrics=v.metrics(ws['reserve'],r['key'],v.emission(np.array(fit['train_counts']),r['key']))) for r in fit['restarts']]
    out['random_label_scores']=randomlabels(ws['reserve'],key,phi)
    out['random_label_p']=(1+sum(s>=out['metrics']['reserve']['objective'] for s in out['random_label_scores']))/1000
    v.dump('results/variable-units-v1/attack.json',out);print('ACTUAL SCORED',out['metrics'],flush=True)
    for trial in range(9):
        rng=np.random.default_rng(20264025+trial)
        def shuffle(rows):
            result=[]
            for r in rows:
                flat=[c for w in r['words'] for c in w if c>=0];values=iter(map(int,rng.permutation(flat)))
                result.extend([[next(values) if c>=0 else -1 for c in w] for w in r['words']])
            return result
        train=shuffle(rs['train']);test=shuffle(rs['reserve']);nf=v.solve(f'null-{trial}',train,G,20265025+trial)
        out['nulls'].append(dict(trial=trial,shuffle_seed=20264025+trial,train_objective=nf['selected']['objective'],development_metrics=v.metrics(test,nf['selected']['key'],nf['emission_probabilities']),fit_sha256=v.digest(f'results/variable-units-v1/null-{trial}.fit.json')))
        v.dump('results/variable-units-v1/attack.json',out);print('NULL',trial,out['nulls'][-1]['development_metrics'],flush=True)
    out['shuffle_p']=(1+sum(r['development_metrics']['objective']>=out['metrics']['reserve']['objective'] for r in out['nulls']))/10
    out['complete']=True;v.dump('results/variable-units-v1/attack.json',out);print('ATTACK COMPLETE',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','attack']);a=p.parse_args();globals()[a.command]()
