"""Experiment010: explicit native-glyph subset, known-text calibration, frozen keys."""
from pathlib import Path
from collections import Counter
from html.parser import HTMLParser
import argparse, hashlib, json, re, subprocess, unicodedata
import xml.etree.ElementTree as ET
import numpy as np
from symbol_sequences import raw_loci, clean

ROOT=Path(__file__).resolve().parents[1]
ALPH='abcdefghiklmnopqrstuxyz'
SPACE=23; UNKNOWN=24
OUT=ROOT/'results/substitution'; DATA=ROOT/'data/latin'
NS={'t':'http://www.tei-c.org/ns/1.0'}

def digest(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dump(p,x):
    p=ROOT/p;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def read(p):return json.loads((ROOT/p).read_text())
def snap(paths):return {p:digest(p) for p in paths}
def normalize(s):
    s=''.join(c for c in unicodedata.normalize('NFKD',s.lower()) if not unicodedata.combining(c))
    s=s.replace('j','i').replace('v','u').replace('æ','ae').replace('œ','oe')
    return ' '.join(re.findall('['+ALPH+']+',s))

def tei_text(e):
    name=e.tag.split('}')[-1]
    if name in {'note','head','foreign','rdg','del','num'}:return ' '
    if name=='choice':
        selected=e.find('t:expan',NS)
        if selected is None:selected=e.find('t:reg',NS)
        return tei_text(selected if selected is not None else e[0])
    return (e.text or '')+''.join(tei_text(c)+(c.tail or '') for c in e)

class PlainHTML(HTMLParser):
    def __init__(self):super().__init__();self.stack=[];self.out=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs);skip=d.get('class','') in {'pagenum','ref','sec','chapter','note'} or tag in {'script','style'}
        if tag not in {'img','br','hr','meta','link','input'}:self.stack.append((tag,skip))
        if tag in {'p','br','div','h1','h2','h3'}:self.out.append(' ')
    def handle_endtag(self,tag):
        for i in range(len(self.stack)-1,-1,-1):
            if self.stack[i][0]==tag:self.stack=self.stack[:i];break
    def handle_data(self,s):
        if not any(x[1] for x in self.stack):self.out.append(s)

def plaintext():
    provenance=read('data/latin-sources.json');counts={}
    for fn,kind in [('caesar-gallic','book'),('cicero-catilinam','speech')]:
        root=ET.parse(ROOT/'data/raw/latin'/f'{fn}.xml').getroot()
        divs=root.findall(f'.//t:body/t:div/t:div[@subtype="{kind}"]',NS)
        for div in divs:
            n=int(div.get('n'))
            if fn=='caesar-gallic' and n>7:continue
            paragraphs=[normalize(tei_text(p)) for p in div.findall('.//t:p',NS)]
            p=DATA/f'{fn}-{n}.txt';p.write_text('\n'.join(paragraphs)+'\n')
            counts[str(p.relative_to(ROOT))]=dict(chars=sum(map(len,paragraphs)),words=sum(len(x.split()) for x in paragraphs),sha256=digest(str(p.relative_to(ROOT))))
    s=(ROOT/'data/raw/latin/celsus-1.html').read_text()
    start=s.index('<A CLASS="chapter" NAME="1">');end=s.index('Adnotationes Editoris:')
    h=PlainHTML();h.feed(s[start:end]);text=normalize(''.join(h.out))
    assert text.startswith('sanus homo qui et bene ualet') and 'marx suggests' not in text
    p=DATA/'celsus-1.txt';p.write_text(text+'\n')
    counts[str(p.relative_to(ROOT))]=dict(chars=len(text),words=len(text.split()),sha256=digest(str(p.relative_to(ROOT))))
    dump('data/latin-extraction.json',dict(files=counts,sources=provenance,tei_rule='body p, choice expan/reg, apparatus lem; omit note/head/foreign/rdg/del/num',celsus_rule='chapter1 anchor through before editorial notes heading; omit page/section/reference numbers',license_note='Perseus XML CC BY-SA4.0; public-domain ancient works/editions. Retain source attribution.'))
    print('PLAINTEXT',counts,flush=True)

def gc_rows():
    pages=[]
    for line in (ROOT/'data/raw/GC2a-n.txt').read_text().splitlines():
        m=re.fullmatch(r'<(f[^.>]+)>\s*(.*)',line)
        if m:
            meta=dict(re.findall(r'\$([A-Z])=(\S)',m[2]))
            if all(meta.get(k)==v for k,v in {'L':'A','H':'1','I':'H'}.items()):pages.append(m[1])
    rows=[]
    for page in pages:
        items=list(raw_loci('GC',page).items())
        for j,(num,v) in enumerate(items):
            if v['code'] not in {'@P0','+P0','*P0'}:continue
            if j+1<len(items) and items[j+1][1]['code'][0] in '=~' and items[j+1][1]['code'][1]=='P':continue
            s=clean(v['raw']).replace('{','').replace('}','')
            words=re.split(r'[.,]',s);atoms=[];word_atoms=[]
            for w in words:
                if not w:continue
                if '?' in w or '[' in w:parts=['?']
                else:
                    parts=re.findall(r'@\d{3};|.',w)
                    assert ''.join(parts)==w
                    assert not any(c in ';@[]:<>/' for p in parts if len(p)==1 for c in p)
                word_atoms.append(parts)
                if atoms:atoms.append(' ')
                atoms.extend(parts)
            rows.append(dict(page=page,locus=f'{page}.{num}',raw=v['raw'],words=word_atoms,atoms=[' ']+atoms+[' '],split='test' if int(hashlib.sha256(page.encode()).hexdigest(),16)%5==0 else 'train'))
    return rows

def grams(lines):
    counts=Counter()
    for line in lines:
        for i in range(len(line)-3):
            q=tuple(line[i:i+4])
            if max(q)<24:counts[q]+=1
    return counts

def lm_build():
    lines=[]
    for p in sorted(DATA.glob('caesar-gallic-*.txt')):
        for line in p.read_text().splitlines():
            if line:lines.append([SPACE]+[ALPH.index(c) if c!=' ' else SPACE for c in line]+[SPACE])
    counts=grams(lines);arr=np.full((24,)*4,.1,dtype=np.float64)
    for q,c in counts.items():arr[q]+=c
    arr=np.log(arr/arr.sum());arr.tofile(DATA/'latin-4gram.bin')
    letters=Counter(c for line in lines for c in line if c!=SPACE)
    rank=sorted(range(23),key=lambda x:(-letters[x],x))
    dump('data/latin-model.json',dict(alphabet=ALPH,space=SPACE,order=4,smoothing=.1,training_files={str(p.relative_to(ROOT)):digest(str(p.relative_to(ROOT))) for p in sorted(DATA.glob('caesar-gallic-*.txt'))},gram_count=sum(counts.values()),letter_counts=dict(letters),letter_rank=rank,model_sha256=digest('data/latin/latin-4gram.bin')))

def prepare():
    OUT.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    plaintext();lm_build();rows=gc_rows()
    counts=Counter(a for r in rows if r['split']=='train' for a in r['atoms'] if a not in {' ','?'})
    glyphs=sorted(counts,key=lambda x:(-counts[x],x))[:23];lookup={v:i for i,v in enumerate(glyphs)}
    for r in rows:r['encoded']=[SPACE if a==' ' else lookup.get(a,UNKNOWN) for a in r['atoms']]
    coverage={}
    for split in ['train','test']:
        rs=[r for r in rows if r['split']==split];letters=[c for r in rs for c in r['encoded'] if c!=SPACE]
        coverage[split]=dict(pages=sorted({r['page'] for r in rs}),loci=len(rs),glyph_slots=len(letters),known=sum(c<23 for c in letters),unknown=sum(c==UNKNOWN for c in letters),score_windows=sum(grams([r['encoded'] for r in rs]).values()))
    dump('data/substitution-manuscript.json',dict(status='FROZEN_BEFORE_DECODING_SCORES',input_hashes=snap(['notes/experiment-010-design.md','data/raw/GC2a-n.txt']),glyphs=glyphs,training_glyph_counts=dict(counts),coverage=coverage,rows=rows))
    config=dict(restarts=16,steps_per_cycle=15000,cycles=3,temperatures=[.03,.0001],seed=20261010)
    dump('data/substitution-config.json',config)
    print('MANUSCRIPT COVERAGE',coverage,'glyphs',glyphs,flush=True)

def solve(name,lines,seed,config):
    g=grams(lines);prefix=OUT/name
    prefix.with_suffix('.grams.txt').write_text(''.join(' '.join(map(str,q))+f' {n}\n' for q,n in sorted(g.items())))
    counts=Counter(c for line in lines for c in line if c<23)
    ranked=[0]*24
    for a,b in zip(sorted(range(23),key=lambda x:(-counts[x],x)),read('data/latin-model.json')['letter_rank']):ranked[a]=b
    ranked[23]=23
    prefix.with_suffix('.rank.txt').write_text(' '.join(map(str,ranked))+'\n')
    args=[str(ROOT/'bin/substitution_search'),str(DATA/'latin-4gram.bin'),str(prefix.with_suffix('.grams.txt')),str(prefix.with_suffix('.rank.txt')),str(seed),str(config['restarts']),str(config['steps_per_cycle']),str(prefix.with_suffix('.keys.txt'))]
    with prefix.with_suffix('.log').open('w') as f:subprocess.run(args,stderr=f,check=True)
    results=[]
    for line in prefix.with_suffix('.keys.txt').read_text().splitlines():
        cells=line.split();results.append(dict(restart=int(cells[0]),train_score=float(cells[1]),key=list(map(int,cells[2:]))))
    best=max(results,key=lambda r:r['train_score'])
    record=dict(name=name,seed=seed,config=config,train_windows=sum(g.values()),restarts=results,selected=best,status='KEY_FROZEN_BEFORE_HELDOUT_SCORING')
    dump(f'results/substitution/{name}.frozen-key.json',record)
    print('SOLVED',name,best['train_score'],flush=True)
    return record

def score(lines,key):
    lm=np.fromfile(DATA/'latin-4gram.bin',dtype=np.float64).reshape((24,)*4);g=grams(lines)
    return sum(n*lm[tuple(key[c] for c in q)] for q,n in g.items())/sum(g.values())

def chunk_text(text):
    words=text.split();return [[SPACE]+[SPACE if c==' ' else ALPH.index(c) for c in ' '.join(words[i:i+8])]+[SPACE] for i in range(0,len(words),8)]

def word_accuracy(cipher,plain,key):
    good=total=0
    for row,truth in zip(cipher,plain):
        start=1
        for end in range(1,len(row)):
            if row[end]!=SPACE:continue
            if end>start and UNKNOWN not in row[start:end]:
                total+=1;good+=all(key[a]==b for a,b in zip(row[start:end],truth[start:end]))
            start=end+1
    return dict(unmasked_words=total,exact_words=good,accuracy=good/total)

def controls():
    config=read('data/substitution-config.json');man=read('data/substitution-manuscript.json')
    cov=man['coverage']['train'];fraction=cov['unknown']/cov['glyph_slots'];out=[]
    cicero=' '.join(p.read_text() for p in sorted(DATA.glob('cicero-catilinam-*.txt')))
    celsus=(DATA/'celsus-1.txt').read_text()
    for trial,(name,text,mask) in enumerate([('control-cicero',cicero,0),('control-celsus',celsus,0),('control-cicero-masked',cicero,fraction)]):
        words=text.split();limit=min(24000,int(len(text)*.7));i=0;n=0
        while i<len(words) and n+len(words[i])+1<=limit:n+=len(words[i])+1;i+=1
        plain_train=chunk_text(' '.join(words[:i]));plain_test=chunk_text(' '.join(words[i:]))
        rng=np.random.default_rng(config['seed']+trial);encryption=list(map(int,rng.permutation(23)))+[SPACE,UNKNOWN]
        truekey=list(map(int,np.argsort(encryption[:23])))+[SPACE]
        def enc(rows):return [[UNKNOWN if c<23 and rng.random()<mask else encryption[c] for c in row] for row in rows]
        train=enc(plain_train);test=enc(plain_test)
        data=dict(name=name,train=train,test=test,plain_train=plain_train,plain_test=plain_test,true_key=truekey,mask_fraction=mask,word_split_index=i)
        dump(f'data/latin/{name}.json',data)
        fitted=solve(name,train,config['seed']+100+trial,config);key=fitted['selected']['key']
        accuracies=[]
        for r in fitted['restarts']:
            k=r['key'];pairs=[(a,b) for row,truth in zip(test,plain_test) for a,b in zip(row,truth) if a<23]
            accuracies.append(dict(restart=r['restart'],train_score=r['train_score'],test_token_accuracy=sum(k[a]==b for a,b in pairs)/len(pairs)))
        supported={c for row in train for c in row if c<23}
        record=dict(name=name,mask_fraction=mask,selected_restart=fitted['selected']['restart'],test_token_accuracy=accuracies[fitted['selected']['restart']]['test_token_accuracy'],supported_key_accuracy=sum(key[c]==truekey[c] for c in supported)/len(supported),supported_key_letters=len(supported),all_restarts=accuracies,test_score=score(test,key),correct_key_test_score=score(test,truekey),train_score=fitted['selected']['train_score'],data_sha256=digest(f'data/latin/{name}.json'),frozen_key_sha256=digest(f'results/substitution/{name}.frozen-key.json'))
        record['test_fullword_recovery']=word_accuracy(test,plain_test,key)
        out.append(record);dump('results/substitution/calibration.json',dict(status='IN_PROGRESS',controls=out))
        print('CONTROL',record,flush=True)
    passed=all(r['test_token_accuracy']>=.95 for r in out)
    dump('results/substitution/calibration.json',dict(status='PASS' if passed else 'FAIL',controls=out,criterion='each selected key >=95% unmasked test letter accuracy'))
    assert passed,'Preserve failures; tune on controls only before manuscript search.'
    paths=['src/substitution.py','src/substitution_search.cpp','data/substitution-config.json','data/substitution-manuscript.json','data/latin-model.json','data/latin/latin-4gram.bin','notes/experiment-010-design.md','results/substitution/calibration.json']
    dump('results/substitution/solver-freeze.json',dict(status='FROZEN_BEFORE_MANUSCRIPT_SEARCH',hashes=snap(paths)))

def attack():
    freeze=read('results/substitution/solver-freeze.json')
    for p,h in freeze['hashes'].items():assert digest(p)==h,(p,'changed since freeze')
    man=read('data/substitution-manuscript.json');config=read('data/substitution-config.json')
    train=[r['encoded'] for r in man['rows'] if r['split']=='train'];test=[r['encoded'] for r in man['rows'] if r['split']=='test']
    fitted=solve('voynich',train,config['seed']+1000,config);key=fitted['selected']['key']
    key_hash=digest('results/substitution/voynich.frozen-key.json');testscore=score(test,key)
    rng=np.random.default_rng(config['seed']+2000)
    randomscores=[score(test,list(map(int,rng.permutation(23)))+[SPACE]) for _ in range(999)]
    out=dict(status='NOT_A_TRANSLATION',glyphs=man['glyphs'],selected_key=key,frozen_key_sha256=key_hash,train_score=fitted['selected']['train_score'],test_score=testscore,random_key_scores=randomscores,random_key_p=(1+sum(x>=testscore for x in randomscores))/1000,nulls=[],restart_test_scores=[dict(restart=r['restart'],train_score=r['train_score'],test_score=score(test,r['key']),key=r['key']) for r in fitted['restarts']])
    def decode(row,k):return ''.join('□' if c==UNKNOWN else ' ' if c==SPACE else ALPH[k[c]] for c in row)
    dump('results/substitution/voynich-decoded.json',dict(status='NOT_A_TRANSLATION; mechanical output under one restricted candidate key',glyph_to_letter=dict(zip(man['glyphs'],[ALPH[x] for x in key[:23]])),rows=[dict(locus=r['locus'],split=r['split'],source=r['raw'],output=decode(r['encoded'],key)) for r in man['rows']]))
    dump('results/substitution/attack.json',out)
    for trial in range(9):
        nrng=np.random.default_rng(config['seed']+3000+trial)
        def shuffled(rows):
            result=[]
            for row in rows:
                idx=[i for i,c in enumerate(row) if c<23];arr=list(row);vals=nrng.permutation([row[i] for i in idx])
                for i,v in zip(idx,vals):arr[i]=int(v)
                result.append(arr)
            return result
        tr=shuffled(train);te=shuffled(test)
        f=solve(f'null-{trial}',tr,config['seed']+4000+trial,config)
        nk=f['selected']['key'];out['nulls'].append(dict(trial=trial,shuffle_seed=config['seed']+3000+trial,train_score=f['selected']['train_score'],test_score=score(te,nk),frozen_key_sha256=digest(f'results/substitution/null-{trial}.frozen-key.json')))
        dump('results/substitution/attack.json',out)
    out['null_test_p']=(1+sum(n['test_score']>=testscore for n in out['nulls']))/10
    out['complete']=True;dump('results/substitution/attack.json',out)
    print('ATTACK', {k:v for k,v in out.items() if k not in {'random_key_scores','restart_test_scores'}},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','controls','attack']);a=p.parse_args()
    {'prepare':prepare,'controls':controls,'attack':attack}[a.stage]()
