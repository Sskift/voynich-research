"""Independent experiment010 recount; does not import search/preparation code."""
from pathlib import Path
from collections import Counter
import hashlib,json,math,re
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
ALPH='abcdefghiklmnopqrstuxyz'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def close(a,b):assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-9),(a,b)

def native_input():
    lines=(ROOT/'data/raw/GC2a-n.txt').read_text().splitlines();logical=[];pending=''
    for line in lines:
        if not line or line.startswith('#'):continue
        if pending:
            assert line.startswith('/');line=pending+line[1:].lstrip();pending=''
        if line.endswith('/'):pending=line[:-1]
        else:logical.append(line)
    assert not pending
    records=[];meta={};page=''
    for line in logical:
        m=re.fullmatch(r'<(f[^.,>]+)>\s*(.*)',line)
        if m:page=m[1];meta=dict(re.findall(r'\$([A-Z])=(\S)',m[2]));continue
        m=re.fullmatch(r'<(f[^.>]+)\.(\d+),([^>]+)>\s*(.*)',line)
        assert m and m[1]==page
        meta.update(dict(re.findall(r'<@([A-Z])=(.)>',m[4])))
        records.append(dict(page=page,locus=page+'.'+m[2],code=m[3],raw=m[4],meta=meta.copy()))
    accepted=[]
    for i,r in enumerate(records):
        if not all(r['meta'].get(k)==v for k,v in {'L':'A','H':'1','I':'H'}.items()):continue
        if r['code'] not in {'@P0','+P0','*P0'}:continue
        if i+1<len(records) and records[i+1]['page']==r['page'] and records[i+1]['code'][0] in '=~' and records[i+1]['code'][1]=='P':continue
        s=r['raw'].replace('<->','.').replace('<~>','.')
        s=re.sub('<[^>]*>','',s);s=''.join(s.split()).replace('{','').replace('}','')
        atoms=[' '];words=[]
        for word in re.split('[.,]',s):
            if not word:continue
            if '?' in word or '[' in word:parts=['?']
            else:
                parts=[];pos=0
                while pos<len(word):
                    if word[pos]=='@':
                        a=word[pos:pos+5];assert re.fullmatch(r'@\d{3};',a) and 128<=int(a[1:4])<=255
                        parts.append(a);pos+=5
                    else:parts.append(word[pos]);pos+=1
            if words:atoms.append(' ')
            atoms.extend(parts);words.append(parts)
        atoms.append(' ');r['atoms']=atoms;r['words']=words
        r['split']='test' if int(hashlib.sha256(r['page'].encode()).hexdigest(),16)%5==0 else 'train'
        accepted.append(r)
    return accepted

def gramcounts(rows):
    counts=Counter()
    for row in rows:
        for a,b,c,d in zip(row,row[1:],row[2:],row[3:]):
            if 24 not in (a,b,c,d):counts[(a,b,c,d)]+=1
    return counts

def scores(g,keys,model):
    q=np.array(list(g),dtype=int);n=np.array(list(g.values()),dtype=float)
    result=[]
    for key in keys:
        assert sorted(key[:23])==list(range(23)) and key[23]==23
        d=np.array(key)[q];idx=((d[:,0]*24+d[:,1])*24+d[:,2])*24+d[:,3]
        result.append(float(np.dot(model[idx],n)/n.sum()))
    return result

def main():
    check=[];freeze=read('results/substitution/solver-freeze.json')
    for p,h in freeze['hashes'].items():assert sha(p)==h,p
    for src in read('data/latin-sources.json')['sources']:assert sha(src['path'])==src['sha256']
    for p,v in read('data/latin-extraction.json')['files'].items():assert sha(p)==v['sha256']
    check.append('frozen design/code/model/data, raw Latin sources and extracted plaintext hashes')
    lmmeta=read('data/latin-model.json');joint=Counter()
    for p,h in lmmeta['training_files'].items():
        assert sha(p)==h
        for line in (ROOT/p).read_text().splitlines():
            text=' '+line+' '
            for i in range(len(text)-3):joint[text[i:i+4]]+=1
    assert sum(joint.values())==lmmeta['gram_count']
    dense=np.full(24**4,.1)
    for word,n in joint.items():
        index=0
        for c in word:index=index*24+(' '+ALPH).index(c)-1 if c!=' ' else index*24+23
        dense[index]+=n
    independently=np.log(dense/dense.sum());model=np.fromfile(ROOT/'data/latin/latin-4gram.bin',dtype=np.float64)
    assert np.allclose(independently,model,rtol=0,atol=1e-13)
    check.append('all331776 Latin fourgram scores reconstructed from text, independent encoding')
    man=read('data/substitution-manuscript.json');raw=native_input();assert len(raw)==len(man['rows'])
    counts=Counter(a for r in raw if r['split']=='train' for a in r['atoms'] if a not in {'?',' '})
    glyphs=sorted(counts,key=lambda c:(-counts[c],c))[:23];assert glyphs==man['glyphs']
    for x,y in zip(raw,man['rows']):
        for field in ['page','locus','raw','atoms','words','split']:assert x[field]==y[field],(field,x['locus'])
        coded=[23 if a==' ' else glyphs.index(a) if a in glyphs else 24 for a in x['atoms']]
        assert coded==y['encoded']
    train=[r['encoded'] for r in man['rows'] if r['split']=='train'];test=[r['encoded'] for r in man['rows'] if r['split']=='test']
    assert set(man['coverage']['train']['pages']).isdisjoint(man['coverage']['test']['pages'])
    for split,rows in [('train',train),('test',test)]:
        assert sum(gramcounts(rows).values())==man['coverage'][split]['score_windows']
        assert sum(c<23 for r in rows for c in r)==man['coverage'][split]['known']
    check.append('all1165 eligible native loci, case/escape units, masks, train-only top23 and page split')
    recounted=0;cal=read('results/substitution/calibration.json')
    def fitted(name,tr):
        nonlocal recounted
        f=read(f'results/substitution/{name}.frozen-key.json');g=gramcounts(tr)
        assert sum(g.values())==f['train_windows'];keys=[r['key'] for r in f['restarts']]
        for r,score in zip(f['restarts'],scores(g,keys,model)):close(r['train_score'],score);recounted+=1
        assert max(f['restarts'],key=lambda r:r['train_score'])==f['selected']
        saved=Counter()
        for line in (ROOT/f'results/substitution/{name}.grams.txt').read_text().splitlines():
            a,b,c,d,n=map(int,line.split());saved[(a,b,c,d)]=n
        assert saved==g
        return f
    errors=[]
    for c in cal['controls']:
        name=c['name'];data=read(f'data/latin/{name}.json');assert sha(f'data/latin/{name}.json')==c['data_sha256']
        f=fitted(name,data['train']);assert sha(f'results/substitution/{name}.frozen-key.json')==c['frozen_key_sha256']
        for split in ['train','test']:
            for row,truth in zip(data[split],data['plain_'+split]):
                assert len(row)==len(truth)
                for a,b in zip(row,truth):
                    if a!=24:assert data['true_key'][a]==b
        key=f['selected']['key'];right=total=goodwords=words=0;wrong=[]
        for row,truth in zip(data['test'],data['plain_test']):
            for a,b in zip(row,truth):
                if a<23:
                    total+=1;right+=key[a]==b
                    if key[a]!=b:wrong.append([ALPH[b],ALPH[key[a]]])
            chunks=[];now=[]
            for a,b in zip(row,truth):
                if a==23:
                    if now:chunks.append(now);now=[]
                else:now.append((a,b))
            for chunk in chunks:
                if any(a==24 for a,b in chunk):continue
                words+=1;goodwords+=all(key[a]==b for a,b in chunk)
        close(c['test_token_accuracy'],right/total);close(c['test_fullword_recovery']['accuracy'],goodwords/words)
        g=gramcounts(data['test']);ss=scores(g,[key,data['true_key']],model)
        close(ss[0],c['test_score']);close(ss[1],c['correct_key_test_score'])
        for r,sc in zip(c['all_restarts'],f['restarts']):
            pairs=[(a,b) for row,truth in zip(data['test'],data['plain_test']) for a,b in zip(row,truth) if a<23]
            close(r['test_token_accuracy'],sum(sc['key'][a]==b for a,b in pairs)/len(pairs))
        errors.append(dict(name=name,test_letters=total,wrong_letters=total-right,substitutions=wrong,unmasked_words=words,wrong_words=words-goodwords))
    assert cal['status']=='PASS' and all(c['test_token_accuracy']>=.95 for c in cal['controls'])
    check.append('three known ciphers, all48 restart test accuracies and complete-word recovery')
    attack=read('results/substitution/attack.json');assert attack['complete'] and len(attack['nulls'])==9
    f=fitted('voynich',train);assert sha('results/substitution/voynich.frozen-key.json')==attack['frozen_key_sha256']
    assert attack['selected_key']==f['selected']['key'];gtest=gramcounts(test)
    close(scores(gtest,[attack['selected_key']],model)[0],attack['test_score'])
    for r,sc in zip(attack['restart_test_scores'],scores(gtest,[r['key'] for r in f['restarts']],model)):close(r['test_score'],sc)
    config=read('data/substitution-config.json');rng=np.random.default_rng(config['seed']+2000)
    randomkeys=[list(map(int,rng.permutation(23)))+[23] for _ in range(999)]
    for a,b in zip(attack['random_key_scores'],scores(gtest,randomkeys,model)):close(a,b)
    close(attack['random_key_p'],(1+sum(s>=attack['test_score'] for s in attack['random_key_scores']))/1000)
    for n in attack['nulls']:
        rng=np.random.default_rng(n['shuffle_seed']);both=[]
        for rows in [train,test]:
            shuffled=[]
            for row in rows:
                idx=np.flatnonzero(np.array(row)<23);values=rng.permutation(np.array(row)[idx]);r=list(row)
                for i,v in zip(idx,values):r[int(i)]=int(v)
                assert Counter(r)==Counter(row);shuffled.append(r)
            both.append(shuffled)
        name=f"null-{n['trial']}";fit=fitted(name,both[0]);assert sha(f'results/substitution/{name}.frozen-key.json')==n['frozen_key_sha256']
        close(fit['selected']['train_score'],n['train_score']);close(scores(gramcounts(both[1]),[fit['selected']['key']],model)[0],n['test_score'])
    close(attack['null_test_p'],(1+sum(n['test_score']>=attack['test_score'] for n in attack['nulls']))/10)
    check.append('all observed/null training keys and scores, all999 random keys, all9 independently regenerated shuffled inputs and test scores')
    decoded=read('results/substitution/voynich-decoded.json');key=attack['selected_key']
    for r,d in zip(man['rows'],decoded['rows']):
        s=''.join('□' if x==24 else ' ' if x==23 else ALPH[key[x]] for x in r['encoded']);assert d['output']==s and d['locus']==r['locus']
    check.append('all mechanical decoded lines retain masks; no semantic translation claim')
    result=dict(status='PASS',checks=check,training_keys_independently_rescored=recounted,control_error_details=errors,completion='UNPROVEN',scope_limit='Does not independently prove a global optimum, authenticate source transcription, or validate Latin/glyph assumptions.')
    (ROOT/'results/audit-007.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
