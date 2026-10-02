"""Independent actual015 audit, using the already independently checked LM."""
from pathlib import Path
from collections import Counter
import hashlib,json,math,re
import numpy as np
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data/variable-units-v1';O=ROOT/'results/variable-units-v1'
ALPH='abcdefghiklmnopqrstuxyz';U=list(ALPH)+[a+b for a in ALPH for b in ALPH]
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-10),(a,b)
def equal(a,b):
    assert a.keys()==b.keys()
    for k in a:close(a[k],b[k])

def main():
    man=read('data/variable-units-v1/manuscript.json');attack=read('results/variable-units-v1/attack.json');assert attack['complete'] and len(attack['nulls'])==9
    frozen=0
    for f in ['data/variable-units-v1/calibration-freeze.json','results/variable-units-v1/solver-freeze.json']:
        for p,h in read(f)['hashes'].items():assert sha(p)==h,p;frozen+=1
    G=len(man['glyphs']);lookup={c:i for i,c in enumerate(man['glyphs'])};source={}
    lines=(ROOT/'data/raw/GC2a-n.txt').read_text().splitlines();i=0
    while i<len(lines):
        m=re.fullmatch(r'<(f[^.,>]+\.\d+),[^>]+>\s*(.*)',lines[i]);i+=1
        if not m:continue
        loc,raw=m.groups()
        while raw.rstrip().endswith('/'):
            assert lines[i].startswith('/');raw=raw.rstrip()[:-1]+lines[i][1:].lstrip();i+=1
        source[loc]=raw
    splits={s:[] for s in ['train','development','reserve']};byrows={s:[] for s in splits}
    for r in man['rows']:
        raw=source[r['locus']];assert raw==r['raw']
        text=re.sub(r'\s+','',re.sub(r'<[^>]*>','',raw.replace('<->','.').replace('<~>','.'))).replace('{','').replace('}','')
        ww=[];atoms=[' ']
        for w in re.split('[.,]',text):
            if not w:continue
            aa=['?'] if '?' in w or '[' in w else re.findall(r'@\d{3};|.',w)
            ww.append([lookup.get(c,-1) for c in aa]);atoms.extend(aa);atoms.append(' ')
        assert ww==r['words'] and atoms==r['atoms'],r['locus']
        splits[r['split']].extend(ww);byrows[r['split']].append(ww)
    for s,ws in splits.items():
        complete=[w for w in ws if min(w)>=0]
        cov=dict(loci=len(byrows[s]),word_slots=len(ws),complete_words=len(complete),incomplete_words=len(ws)-len(complete),complete_input_units=sum(map(len,complete)),all_encoded_positions=sum(map(len,ws)),known_units_outside_complete_words=sum(c>=0 for w in ws if min(w)<0 for c in w),unknown_positions=sum(c<0 for w in ws for c in w))
        assert cov==man['coverage'][s]==attack['coverage'][s]
    lm=np.fromfile(D/'word-conditional4.bin').reshape((24,)*4);terms=np.fromfile(D/'length-terms.bin').reshape(65,129);lang=read('data/variable-units-v1/language.json')
    for p,h in lang['source_hashes'].items():assert sha(p)==h,p;frozen+=1
    for name in ['letter23','pair127','homophone127']:
        control=read(f'data/variable-units-v1/{name}.json');assert sha(control['source'])==control['source_sha256'];frozen+=1
    def phi(cc,key):
        total=Counter()
        for c,u in enumerate(key):total[u]+=cc[c]
        return [cc[c]/total[u] for c,u in enumerate(key)]
    def evaluate(ws,key,prob):
        lp=tp=ep=mass=out=wn=hits=0
        for w,num in Counter(tuple(w) for w in ws if min(w)>=0).items():
            text=''.join(U[key[c]] for c in w);seq=[23,23,23]+[ALPH.index(c) for c in text]+[23]
            lp+=num*sum(float(lm[tuple(seq[i-3:i+1])]) for i in range(3,len(seq)))
            tp+=num*terms[len(w),len(text)];ep+=num*sum(math.log(prob[c]) for c in w)
            mass+=num*len(w);out+=num*len(text);wn+=num;hits+=num*(text in lang['word_counts'])
        return dict(objective=(lp+tp+ep)/mass,word_language_per_input_unit=lp/mass,length_normalization_per_input_unit=tp/mass,emission_per_input_unit=ep/mass,input_units=mass,output_letters=out,expansion_ratio=out/mass,complete_words=wn,word_hits=hits,source_word_fraction=hits/wn)
    checked=0
    def fitcheck(name,words):
        nonlocal checked
        fit=read(f'results/variable-units-v1/{name}.fit.json');cc=Counter(c for w in words for c in w);freq=[cc[i] for i in range(G)];assert freq==fit['train_counts']
        grams=Counter(tuple(w) for w in words);archived={}
        for line in (O/(name+'.words.txt')).read_text().splitlines():
            xs=list(map(int,line.split()));assert xs[0]==len(xs)-2;archived[tuple(xs[2:])]=xs[1]
        assert grams==archived and fit['input_units']==sum(freq) and fit['word_tokens']==len(words) and fit['word_types']==len(grams)
        for r in fit['restarts']:
            key=r['key'];assert len(key)==G and all(0<=u<len(U) for u in key)
            close(evaluate(words,key,phi(freq,key))['objective'],r['objective']);checked+=1
        assert max(fit['restarts'],key=lambda r:r['objective'])==fit['selected']
        for p,q in zip(phi(freq,fit['selected']['key']),fit['emission_probabilities']):close(p,q)
        return fit
    fit=fitcheck('voynich',splits['train']);key=fit['selected']['key'];prob=fit['emission_probabilities']
    assert key==attack['key'] and prob==attack['emission_probabilities'] and sha('results/variable-units-v1/voynich.fit.json')==attack['fit_sha256']
    for s,words in splits.items():equal(evaluate(words,key,prob),attack['metrics'][s])
    for r,stored in zip(fit['restarts'],attack['restart_development_metrics']):
        assert r['restart']==stored['restart'];close(r['objective'],stored['objective']);equal(evaluate(splits['reserve'],r['key'],phi(fit['train_counts'],r['key'])),stored['metrics'])
    # Sum all repeated ngrams without word-type aggregation for the random renames.
    ct=Counter()
    for word in splits['reserve']:
        if min(word)<0:continue
        s=[23]*3+[ALPH.index(c) for x in word for c in U[key[x]]]+[23]
        ct.update(tuple(s[i-3:i+1]) for i in range(3,len(s)))
    qs=np.array(list(ct));multipliers=np.array(list(ct.values()));rng=np.random.default_rng(20263025)
    constant=attack['metrics']['reserve']['length_normalization_per_input_unit']+attack['metrics']['reserve']['emission_per_input_unit']
    for j,expected in enumerate(attack['random_label_scores']):
        perm=list(map(int,rng.permutation(23)))+[23];q=np.array(perm)[qs]
        actual=float(np.dot(lm[q[:,0],q[:,1],q[:,2],q[:,3]],multipliers)/attack['metrics']['reserve']['input_units']+constant);close(actual,expected)
        if j<3:
            kk=[U.index(''.join(ALPH[perm[ALPH.index(c)]] for c in U[u])) for u in key]
            close(evaluate(splits['reserve'],kk,prob)['objective'],expected)
    close(attack['random_label_p'],(1+sum(s>=attack['metrics']['reserve']['objective'] for s in attack['random_label_scores']))/1000)
    for rec in attack['nulls']:
        rng=np.random.default_rng(rec['shuffle_seed']);both=[]
        for split in ['train','reserve']:
            changed=[]
            for row in byrows[split]:
                flat=[c for w in row for c in w];where=[i for i,c in enumerate(flat) if c>=0];new=list(flat)
                for i,c in zip(where,rng.permutation([flat[i] for i in where])):new[i]=int(c)
                assert Counter(new)==Counter(flat);offset=0
                for w in row:changed.append(new[offset:offset+len(w)]);offset+=len(w)
            both.append(changed)
        name='null-'+str(rec['trial']);nf=fitcheck(name,both[0]);assert sha(f'results/variable-units-v1/{name}.fit.json')==rec['fit_sha256']
        close(nf['selected']['objective'],rec['train_objective']);equal(evaluate(both[1],nf['selected']['key'],nf['emission_probabilities']),rec['development_metrics'])
    close(attack['shuffle_p'],(1+sum(r['development_metrics']['objective']>=attack['metrics']['reserve']['objective'] for r in attack['nulls']))/10)
    decoded=read('results/variable-units-v1/voynich-decoded.json');assert decoded['status']=='NOT_A_TRANSLATION';assert decoded['glyph_to_output']==dict(zip(man['glyphs'],[U[u] for u in key]))
    assert len(decoded['rows'])==len(man['rows'])
    for r,out in zip(man['rows'],decoded['rows']):
        assert out==dict(locus=r['locus'],split=r['split'],source=r['raw'],words=[''.join(U[key[c]] if c>=0 else '□' for c in w) for w in r['words']],whole_word_score_mask=[min(w)<0 for w in r['words']])
    result=dict(status='PASS',meaning='Reproduction audit, not a translation or global-optimum certificate.',frozen_hashes=frozen,native_loci=len(man['rows']),recomputed_fit_objectives=checked,random_label_scores=len(attack['random_label_scores']),regenerated_shuffled_pairs=len(attack['nulls']),mechanical_outputs=len(decoded['rows']))
    (ROOT/'results/audit-011-attack.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
