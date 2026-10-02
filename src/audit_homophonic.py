"""Independent input, probability, emission, recovery and null recount for011."""
from pathlib import Path
from collections import Counter
import hashlib,json,math,re
import numpy as np

ROOT=Path(__file__).resolve().parents[1];ALPH='abcdefghiklmnopqrstuxyz'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def close(x,y):assert math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-9),(x,y)

def raw_input():
    logical=[];pending=''
    for s in (ROOT/'data/raw/GC2a-n.txt').read_text().splitlines():
        if not s or s.startswith('#'):continue
        if pending:assert s[0]=='/';s=pending+s[1:].lstrip();pending=''
        if s.endswith('/'):pending=s[:-1]
        else:logical.append(s)
    assert not pending
    records=[];meta={};page=''
    for s in logical:
        m=re.fullmatch(r'<(f[^.,>]+)>\s*(.*)',s)
        if m:page=m[1];meta=dict(re.findall(r'\$([A-Z])=(\S)',m[2]));continue
        m=re.fullmatch(r'<(f[^.>]+)\.(\d+),([^>]+)>\s*(.*)',s);assert m and page==m[1]
        meta.update(dict(re.findall(r'<@([A-Z])=(.)>',m[4])))
        records.append(dict(page=page,locus=page+'.'+m[2],code=m[3],raw=m[4],meta=meta.copy()))
    result={}
    for i,r in enumerate(records):
        if r['meta'].get('L')!='A' or r['meta'].get('H')!='1' or r['meta'].get('I') not in {'H','P'}:continue
        if r['code'] not in {'@P0','+P0','*P0'}:continue
        if i+1<len(records) and records[i+1]['page']==r['page'] and records[i+1]['code'][0] in '=~' and records[i+1]['code'][1]=='P':continue
        s=r['raw'].replace('<->','.').replace('<~>','.');s=re.sub('<[^>]*>','',s);s=''.join(s.split()).replace('{','').replace('}','')
        atoms=[' '];first=True
        for w in re.split('[.,]',s):
            if not w:continue
            if not first:atoms.append(' ')
            first=False
            if '?' in w or '[' in w:atoms.append('?');continue
            j=0
            while j<len(w):
                if w[j]=='@':
                    token=w[j:j+5];assert re.fullmatch(r'@\d{3};',token) and 128<=int(token[1:4])<=255;atoms.append(token);j+=5
                else:atoms.append(w[j]);j+=1
        atoms.append(' ');r['atoms']=atoms;result[r['locus']]=r
    return result

def count4(rows,G):
    ct=Counter()
    for row in rows:
        for q in zip(row,row[1:],row[2:],row[3:]):
            if G+1 not in q:ct[q]+=1
    return ct

def emission(counts,key):
    totals=Counter()
    for i,n in enumerate(counts):totals[key[i]]+=n
    result=[counts[i]/totals[key[i]] for i in range(len(counts))]
    for letter in totals:close(sum(p for i,p in enumerate(result) if key[i]==letter),1)
    return result

def evaluate(rows,key,prob,G,model,lang):
    ct=count4(rows,G);q=np.array(list(ct));n=np.array(list(ct.values()));plain=np.array(key)[q]
    index=((plain[:,0]*24+plain[:,1])*24+plain[:,2])*24+plain[:,3];lm=float(np.dot(model[index],n));mass=int(n.sum())
    c=Counter(a for row in rows for a in row if a<G);channel=sum(n*math.log(prob[a]) for a,n in c.items())
    hit6=total6=hitw=totalw=0
    for row in rows:
        s=''.join('□' if a==G+1 else ' ' if a==G else ALPH[key[a]] for a in row)
        for gram in zip(s,s[1:],s[2:],s[3:],s[4:],s[5:]):
            if '□' not in gram:total6+=1;hit6+=''.join(gram) in lang['sixgram_counts']
        for word in s.split():
            if len(word)>=3 and '□' not in word:totalw+=1;hitw+=word in lang['word_counts']
    return dict(conditional4=lm/mass,emission_per_window=channel/mass,combined=(lm+channel)/mass,windows=mass,sixgram_seen=hit6/total6 if total6 else None,sixgram_hits=hit6,sixgram_windows=total6,word_seen=hitw/totalw if totalw else None,word_hits=hitw,complete_words_ge3=totalw)

def metrics_equal(a,b):
    assert set(a)==set(b)
    for k in a:
        if a[k] is None:assert b[k] is None
        else:close(a[k],b[k])

def main():
    checks=[];freeze=read('results/homophonic/solver-freeze.json')
    for p,h in freeze['hashes'].items():assert sha(p)==h,p
    previous=read('results/substitution/solver-freeze.json')
    for p,h in previous['hashes'].items():assert sha(p)==h,p
    checks.append('previous and current frozen model/code/design/input hashes unchanged')
    language=read('data/homophonic/language.json');counts=[Counter() for _ in range(4)];vocab=Counter();six=Counter()
    for p,h in language['training_files'].items():
        assert sha(p)==h
        for line in (ROOT/p).read_text().splitlines():
            vocab.update(line.split());s=' '+line+' '
            for n in range(1,5):counts[n-1].update(s[i:i+n] for i in range(len(s)-n+1))
            six.update(s[i:i+6] for i in range(len(s)-5))
    assert dict(vocab)==language['word_counts'] and dict(six)==language['sixgram_counts']
    p=None
    for n,ct in enumerate(counts,1):
        arr=np.zeros((24,)*n)
        for s,freq in ct.items():arr[tuple(23 if c==' ' else ALPH.index(c) for c in s)]=freq
        p=(arr+.1)/(arr.sum()+2.4) if n==1 else (arr+10*p[np.newaxis,...])/(arr.sum(axis=-1)[...,np.newaxis]+10)
        assert np.allclose(p.sum(axis=-1),1,rtol=0,atol=1e-14)
    model=np.fromfile(ROOT/'data/homophonic/conditional4.bin',dtype=np.float64)
    assert np.allclose(np.log(p).ravel(),model,rtol=0,atol=1e-13)
    checks.append('all331776 conditional scores, context normalization, word and sixgram reference counts')
    man=read('data/homophonic/manuscript.json');G=len(man['glyphs']);raw=raw_input();assert len(raw)==len(man['rows'])
    counts=Counter(a for r in man['rows'] if r['split']=='train' for a in raw[r['locus']]['atoms'] if a not in {' ','?'})
    assert sorted(counts,key=lambda x:(-counts[x],x))==man['glyphs']
    lookup={a:i for i,a in enumerate(man['glyphs'])};splits={s:[] for s in ['train','development','reserve']}
    for r in man['rows']:
        x=raw[r['locus']];assert x['raw']==r['raw'] and x['atoms']==r['atoms']
        if x['meta']['I']=='P':expected='reserve'
        else:expected='development' if int(hashlib.sha256(x['page'].encode()).hexdigest(),16)%5==0 else 'train'
        assert r['split']==expected
        coded=[G if a==' ' else lookup.get(a,G+1) for a in x['atoms']];assert coded==r['encoded'];splits[expected].append(coded)
    assert set(man['coverage']['reserve']['pages'])==set('f88r f88v f89r1 f89r2 f89v2 f89v1 f99r f99v f100r f100v f101r f101v f102r1 f102r2 f102v2 f102v1'.split())
    for split,rows in splits.items():
        cov=man['coverage'][split];assert len(rows)==cov['loci'];assert sum(count4(rows,G).values())==cov['windows']
        assert sum(c<G for row in rows for c in row)==cov['known'];assert sum(c==G+1 for row in rows for c in row)==cov['masked']
    checks.append('all1388 raw native loci,127 train glyphs, masks,16 reserved P pages, prior H split')
    recounted=0
    def fitcheck(name,rows,size):
        nonlocal recounted
        f=read(f'results/homophonic/{name}.frozen.json');count=Counter(c for row in rows for c in row if c<size);frequencies=[count[i] for i in range(size)]
        assert frequencies==f['train_counts'] and min(frequencies)>0
        q=count4(rows,size);saved={}
        for line in (ROOT/f'results/homophonic/{name}.grams.txt').read_text().splitlines():
            *gram,n=map(int,line.split());saved[tuple(gram)]=n
        assert q==saved
        for r in f['restarts']:
            assert len(r['key'])==size+1 and r['key'][size]==23 and all(0<=k<23 for k in r['key'][:size])
            prob=emission(frequencies,r['key']);value=evaluate(rows,r['key'],prob,size,model,language)
            close(r['objective'],value['combined']);recounted+=1
        assert max(f['restarts'],key=lambda x:x['objective'])==f['selected']
        for a,b in zip(emission(frequencies,f['selected']['key']),f['emission_probabilities']):close(a,b)
        return f
    cal=read(f"results/homophonic/calibration-{freeze['version']}.json");assert cal['status']=='PASS';caldetails=[]
    for rec in cal['controls']:
        name=rec['name'];data=read(f'data/homophonic/{name}.json');size=data['glyph_count'];assert sha(f'data/homophonic/{name}.json')==rec['data_sha256']
        f=fitcheck(name,data['train'],size);assert sha(f'results/homophonic/{name}.frozen.json')==rec['fit_sha256']
        for split in ['train','test']:
            for row,truth in zip(data[split],data['plain_'+split]):
                assert len(row)==len(truth)
                for a,b in zip(row,truth):
                    if a<=size:assert data['true_key'][a]==b
        for restart,r in zip(f['restarts'],rec['restart_accuracy']):
            pairs=[(a,b) for row,truth in zip(data['test'],data['plain_test']) for a,b in zip(row,truth) if a<size]
            right=sum(restart['key'][a]==b for a,b in pairs);close(right/len(pairs),r['letter_accuracy'])
        key=f['selected']['key'];a=rec['accuracy'];selected=rec['restart_accuracy'][f['selected']['restart']]
        for k in a:close(a[k],selected[k])
        pairs=[(c,t) for row,truth in zip(data['test'],data['plain_test']) for c,t in zip(row,truth) if c<size]
        predicted=len(pairs);correct=sum(key[c]==t for c,t in pairs)
        all_letters=sum(t<23 for row in data['plain_test'] for t in row)
        assert a['predicted_letters']==predicted and a['correct_letters']==correct and a['total_letters']==all_letters
        close(a['predicted_fraction'],predicted/all_letters)
        close(a['key_accuracy'],sum(key[c]==data['true_key'][c] for c in range(size))/size)
        total=right=0
        for row,truth in zip(data['test'],data['plain_test']):
            word=[]
            for c,t in zip(row,truth):
                if c==size:
                    if word and all(x<size for x,y in word):total+=1;right+=all(key[x]==y for x,y in word)
                    word=[]
                else:word.append((c,t))
        close(right/total,a['word_accuracy']);assert a['correct_words']==right and a['complete_words']==total
        assert a['letter_accuracy']>=.95
        metrics_equal(evaluate(data['test'],key,f['emission_probabilities'],size,model,language),rec['test_metrics'])
        trueprob=emission(f['train_counts'],data['true_key']);metrics_equal(evaluate(data['test'],data['true_key'],trueprob,size,model,language),rec['true_key_metrics'])
        caldetails.append(dict(name=name,observed_training_glyphs=size,correct=a['correct_letters'],predicted=a['predicted_letters'],all_test_letters=a['total_letters']))
    checks.append('three many-to-one known ciphers; all48 restart letter accuracies and selected fullword metrics')
    attack=read('results/homophonic/attack.json');assert attack['complete'] and len(attack['nulls'])==9
    fit=fitcheck('voynich',splits['train'],G);assert sha('results/homophonic/voynich.frozen.json')==attack['frozen_key_sha256']
    assert fit['selected']['key']==attack['key'];key=attack['key'];prob=attack['emission_probabilities']
    assert prob==fit['emission_probabilities']
    for split,rows in splits.items():metrics_equal(evaluate(rows,key,prob,G,model,language),attack['metrics'][split])
    for f,r in zip(fit['restarts'],attack['restart_reserve_metrics']):
        assert f['key']==r['key'];metrics_equal(evaluate(splits['reserve'],f['key'],emission(fit['train_counts'],f['key']),G,model,language),r['metrics'])
    config=read(f"data/homophonic/config-{freeze['version']}.json");rng=np.random.default_rng(config['seed']+2000)
    ct=count4(splits['reserve'],G);q=np.array(list(ct));n=np.array(list(ct.values()))
    for expected in attack['random_label_scores']:
        permutation=rng.permutation(23);k=np.array([int(permutation[v]) for v in key[:G]]+[23]);x=k[q];idx=((x[:,0]*24+x[:,1])*24+x[:,2])*24+x[:,3]
        close(expected,float(np.dot(model[idx],n)/n.sum()))
    close(attack['random_label_p'],(1+sum(s>=attack['metrics']['reserve']['conditional4'] for s in attack['random_label_scores']))/1000)
    for trial in attack['nulls']:
        rng=np.random.default_rng(trial['shuffle_seed']);both=[]
        for rows in [splits['train'],splits['reserve']]:
            new=[]
            for row in rows:
                positions=np.flatnonzero(np.array(row)<G);values=rng.permutation(np.array(row)[positions]);rr=list(row)
                for i,c in zip(positions,values):rr[int(i)]=int(c)
                assert Counter(rr)==Counter(row);new.append(rr)
            both.append(new)
        name=f"null-{trial['trial']}";f=fitcheck(name,both[0],G);assert sha(f'results/homophonic/{name}.frozen.json')==trial['fit_sha256']
        close(f['selected']['objective'],trial['train_objective']);metrics_equal(evaluate(both[1],f['selected']['key'],f['emission_probabilities'],G,model,language),trial['reserve_metrics'])
    close(attack['null_conditional_p'],(1+sum(r['reserve_metrics']['conditional4']>=attack['metrics']['reserve']['conditional4'] for r in attack['nulls']))/10)
    checks.append('all208 fitted key/emission objectives; all999 random label permutations; all9 regenerated shuffled inputs and reserve diagnostics')
    old=read('data/substitution-manuscript.json');oldglyphs=old['glyphs'];common=[]
    for r in man['rows']:
        if r['split']=='reserve':common.append([23 if a==' ' else oldglyphs.index(a) if a in oldglyphs else 24 for a in r['atoms']])
    keys={'old_frozen_bijective_key':read('results/substitution/voynich.frozen-key.json')['selected']['key'],'new_frozen_homophonic_key':[key[man['glyphs'].index(a)] for a in oldglyphs]+[23]}
    for label,k in keys.items():metrics_equal(evaluate(common,k,[1]*23,23,model,language),attack['common23_candidate_comparison'][label])
    decoded=read('results/homophonic/voynich-decoded.json');assert len(decoded['rows'])==len(man['rows'])
    for r,d in zip(man['rows'],decoded['rows']):
        assert r['locus']==d['locus'] and d['output']==''.join('□' if c==G+1 else ' ' if c==G else ALPH[key[c]] for c in r['encoded'])
    checks.append('common23 mask matched comparison and all1388 mechanical outputs; no translation supplied')
    result=dict(status='PASS',checks=checks,training_keys_independently_rescored=recounted,controls=caldetails,completion='UNPROVEN',scope_limit='Does not certify a global optimum, allographic equivalence, language identity or manuscript meaning.')
    (ROOT/'results/audit-008.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
