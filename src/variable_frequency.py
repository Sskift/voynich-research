"""Experiment016: exact necessary frequency bounds, not a decipherment."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from itertools import product
import argparse,hashlib,json,unicodedata
import numpy as np

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/variable-frequency';OUT=ROOT/'results/variable-frequency'
ALPH='abcdefghiklmnopqrstuxyz'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dump(p,x):
    p=ROOT/p;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def rat(x):return dict(numerator=x.numerator,denominator=x.denominator,value=float(x))
def frac(x):return F(x['numerator'],x['denominator'])
def verify(p):
    d=read(p);hs=d.get('hashes',d.get('files'))
    for f,h in hs.items():assert sha(f)==h,f
def snapshot(paths):return {p:sha(p) for p in paths}
def opaque(u):return any(unicodedata.category(c)=='Co' or c.isdigit() or c in '~&⁊' for c in u)
def blocks(words):
    result=[];start=size=0
    for i,w in enumerate(words):
        size+=len(w)
        if size>=5000:result.append((start,i+1));start=i+1;size=0
    return result,dict(start_word=start,words=len(words)-start,units=size)
def profile(words):
    units=Counter(c for w in words for c in w);N=sum(units.values());assert N>0
    clear={c:n for c,n in units.items() if not opaque(c)};unknown={c:n for c,n in units.items() if opaque(c)}
    pairs=Counter(tuple(w[i:i+2]) for w in words for i in range(len(w)-1))
    cp={p:n for p,n in pairs.items() if not any(opaque(c) for c in p)}
    U=sum(unknown.values());V=sum(n for p,n in pairs.items() if p not in cp)
    C=max(clear.values(),default=0);P=max(cp.values(),default=0);a=F(C+U,N);rawb=F(P+V,N);b=min(rawb,a)
    return dict(units=N,words=len(words),unit_counts=dict(sorted(units.items())),pair_counts=[dict(pair=list(p),count=n) for p,n in sorted(pairs.items())],opaque_counts=unknown,opaque_units=U,opaque_pair_positions=V,max_clear_letter_count=C,max_clear_pair_count=P,letter_maximizers=sorted(c for c,n in clear.items() if n==C),pair_maximizers=[list(p) for p,n in sorted(cp.items()) if n==P],a_upper=rat(a),b_raw_upper=rat(rawb),b_upper=rat(b))
def bound(counts,a,b):
    assert 0<a<=1 and 0<=b<=a and all(v>0 for v in counts.values())
    N=sum(counts.values());freq={c:F(n,N) for c,n in counts.items()};top=max(freq.values());lower=max(F(1),top/a);t=F(2);steps=[]
    for _ in range(len(counts)+2):
        forced=sorted(c for c,f in freq.items() if f>t*b);mass=sum((freq[c] for c in forced),F(0));new=min(t,2-mass)
        steps.append(dict(upper_before=rat(t),forced_monographs=forced,forced_mass=rat(mass),upper_after=rat(new)))
        if t==new:break
        t=new
    else:raise AssertionError('nonconvergence')
    return dict(cipher_units=N,cipher_types=len(counts),cipher_counts=dict(sorted(counts.items())),max_cipher_frequency=rat(top),maximizers=sorted(c for c,f in freq.items() if f==top),a=rat(a),b=rat(b),lower=rat(lower),upper=rat(t),rejects=lower>t,iterations=steps)
def envelope(sources,kind):
    entries=[]
    for name,obj in sources.items():
        if kind=='blocks':entries.extend((name,str(r['index']),r['profile']) for r in obj['blocks'])
        else:entries.append((name,'full',obj['full']))
    aa=max(frac(p['a_upper']) for _,_,p in entries);bb=max(frac(p['b_upper']) for _,_,p in entries)
    return dict(kind=kind,entries=len(entries),a=rat(aa),b=rat(bb),a_witnesses=[dict(source=n,scope=i) for n,i,p in entries if frac(p['a_upper'])==aa],b_witnesses=[dict(source=n,scope=i) for n,i,p in entries if frac(p['b_upper'])==bb])

def freeze():
    assert not (DATA/'input-freeze.json').exists()
    for p in ['data/language-screen/source-freeze.json','data/variable-units-v1/calibration-freeze.json','results/variable-units-v1/solver-freeze.json']:verify(p)
    paths=['notes/experiment-016-design.md','notes/experiment-016-rules.md','src/variable_frequency.py','data/language-screen/macer.json','data/language-screen/aquinas.json','data/homophonic/manuscript.json','data/latin/control-cicero.json','data/latin/control-celsus.json']
    paths += [f'data/latin/caesar-gallic-{n}.txt' for n in range(1,8)]
    paths += [f'data/variable-units-v1/{n}.json' for n in ['letter23','pair127','homophone127']]
    dump('data/variable-frequency/input-freeze.json',dict(status='FROZEN_BEFORE_REFERENCE_PAIR_COUNTS',hashes=snapshot(paths)))
    print('FROZEN',len(paths),flush=True)

def toys():
    verify('data/variable-frequency/input-freeze.json');cases=assignments=feasiblecases=falseadmit=0;witness=None
    for G in range(1,5):
        for cc in product(range(1,4),repeat=G):
            N=sum(cc);counts={str(i):n for i,n in enumerate(cc)}
            for ia in range(1,13):
                for ib in range(ia+1):
                    a=F(ia,12);b=F(ib,12);res=bound(counts,a,b);lo=frac(res['lower']);hi=frac(res['upper']);feasible=[];cases+=1
                    for lens in product([1,2],repeat=G):
                        length=sum(c*l for c,l in zip(cc,lens));assignments+=1
                        if all(c*12<=length*(ia if l==1 else ib) for c,l in zip(cc,lens)):
                            lam=F(length,N);assert lo<=lam<=hi;(feasible.append(lam))
                    if feasible:feasiblecases+=1;assert not res['rejects']
                    elif not res['rejects']:
                        falseadmit+=1
                        if witness is None:witness=dict(counts=list(cc),a=rat(a),b=rat(b),lower=rat(lo),upper=rat(hi))
    codebooks=[['a','b'],['a','b','aa','ab','ba','bb'],['a','b','aa','ab','ba','bb','a'],['a','b','ab']];ciphercases=0
    for units in codebooks:
        for g in range(1,5):
            for cipher in product(range(len(units)),repeat=g):
                w=list(''.join(units[c] for c in cipher));p=profile([w]);a=frac(p['a_upper']);b=frac(p['b_upper']);ct=Counter(map(str,cipher));r=bound(ct,a,b);lam=F(len(w),g)
                assert frac(r['lower'])<=lam<=frac(r['upper']) and not r['rejects']
                for c,n in ct.items():assert F(n,g)<=lam*(a if len(units[int(c)])==1 else b)
                ciphercases+=1
    result=dict(status='PASS',length_distribution_cases=cases,length_assignments=assignments,cases_with_feasible_assignment=feasiblecases,false_admissions=falseadmit,false_admission_witness=witness,exact_decoded_cipher_cases=ciphercases,note='False admissions are expected from a necessary relaxation. No false rejections observed in exhaustive finite cases; mathematical proof still supplies general validity.')
    dump('results/variable-frequency/toys.json',result);print(json.dumps(result,indent=2),flush=True)

def reference_words():
    caesar=[list(w) for n in range(1,8) for w in (ROOT/f'data/latin/caesar-gallic-{n}.txt').read_text().split()]
    return dict(caesar=caesar,**{name:[w for r in read(f'data/language-screen/{name}.json')['records'] for w in r['words']] for name in ['macer','aquinas']})
def references():
    verify('data/variable-frequency/input-freeze.json');assert read('results/variable-frequency/toys.json')['status']=='PASS';sources={}
    for name,words in reference_words().items():
        bs,tail=blocks(words);sources[name]=dict(full=profile(words),tail=tail,blocks=[dict(index=i,start_word=s,end_word=e,profile=profile(words[s:e])) for i,(s,e) in enumerate(bs)])
    envelopes={k:envelope(sources,k) for k in ['blocks','full']}
    dump('data/variable-frequency/references.json',dict(status='SOURCE_ONLY_ENVELOPES',sources=sources,envelopes=envelopes))
    paths=['data/variable-frequency/input-freeze.json','results/variable-frequency/toys.json','data/variable-frequency/references.json']
    dump('data/variable-frequency/envelope-freeze.json',dict(status='FROZEN_BEFORE_CONTROL_TRANSFER',hashes=snapshot(paths)))
    print('ENVELOPES',json.dumps(envelopes),flush=True)

def plainrows(rows):
    return [''.join(ALPH[c] if c<23 else ' ' for c in row).split() for row in rows]
def segment(w,pairs):
    out=[]
    while w:
        u=w[:2] if w[:2] in pairs else w[:1];out.append(u);w=w[len(u):]
    return out
def regenerate():
    caesar=[''.join(w) for w in reference_words()['caesar']];pairs=Counter(w[i:i+2] for w in caesar for i in range(len(w)-1));out={}
    for trial,(name,source,npairs,extras) in enumerate([('letter23','cicero',0,0),('pair127','cicero',104,0),('homophone127','celsus',40,64)]):
        chosen=sorted(pairs,key=lambda u:(-pairs[u],u))[:npairs];units=list(ALPH)+chosen;sf=Counter(u for w in caesar for u in segment(w,set(chosen)))
        f=np.array([sf[u] for u in units],dtype=float);f/=f.sum();ex=np.floor(extras*f).astype(int)
        for i in sorted(range(len(units)),key=lambda i:(-(extras*f[i]-ex[i]),units[i]))[:extras-int(ex.sum())]:ex[i]+=1
        rng=np.random.default_rng(20261025+trial);labels=list(map(int,rng.permutation(len(units)+extras)));groups={};off=0
        for i,u in enumerate(units):groups[u]=labels[off:off+1+ex[i]];off+=1+int(ex[i])
        truth={str(c):u for u,cs in groups.items() for c in cs};old=read(f'data/variable-units-v1/{name}.json');assert groups==old['groups'] and chosen==old['pairs']
        plain=read(f'data/latin/control-{source}.json');out[name]=dict(key=truth,source=source)
        for split in ['train','test']:
            words=[w for row in plainrows(plain['plain_'+split]) for w in row];rows=[];assert len(words)==len(old[split])
            for w,stored in zip(words,old[split]):
                us=segment(w,set(chosen));codes=[]
                for u in us:
                    gr=groups[u];codes.append(gr[0] if len(gr)==1 or rng.random()<.9 else int(rng.choice(gr[1:])))
                assert codes==stored['labels'] and w==stored['plain'] and us==stored['units']
                rows.append(dict(plain=list(w),cipher=list(map(str,codes))))
            out[name][split]=rows
    return out
def check_transfer(rows,key,env):
    cc=Counter(c for r in rows for c in r['cipher']);N=sum(cc.values());P=sum(len(r['plain']) for r in rows);lam=F(P,N);a=frac(env['a']);b=frac(env['b']);res=bound(cc,a,b)
    failures=[dict(code=c,output=key[c],cipher_count=n,required_frequency=rat(F(n,P)),allowed_frequency=rat(a if len(key[c])==1 else b)) for c,n in sorted(cc.items()) if F(n,N)>lam*(a if len(key[c])==1 else b)]
    length_pass=frac(res['lower'])<=lam<=frac(res['upper'])
    return dict(words=len(rows),plaintext_units=P,true_expansion=rat(lam),bound=res,true_expansion_admitted=length_pass,known_code_violations=failures,pass_transfer=length_pass and not failures)
def controls():
    verify('data/variable-frequency/input-freeze.json');verify('data/variable-frequency/envelope-freeze.json');envs=read('data/variable-frequency/references.json')['envelopes'];reg=regenerate();results=[];decodedwords=0
    for name,data in reg.items():
        own={};transfer={};sens=[]
        for split in ['train','test']:
            rows=data[split];decodedwords+=len(rows)
            for r in rows:assert list(''.join(data['key'][c] for c in r['cipher']))==r['plain']
            p=profile([r['plain'] for r in rows]);env=dict(a=p['a_upper'],b=p['b_upper']);own[split]=check_transfer(rows,data['key'],env);assert own[split]['pass_transfer']
        for kind,env in envs.items():transfer[kind]=check_transfer(data['test'],data['key'],env)
        bs,tail=blocks([r['plain'] for r in data['test']])
        for i,(s,e) in enumerate(bs):sens.append(dict(index=i,start_word=s,end_word=e,**check_transfer(data['test'][s:e],data['key'],envs['blocks'])))
        results.append(dict(name=name,source=data['source'],paired=own,full_test_transfer=transfer,test_block_transfer=sens,test_block_tail=tail))
    primary=all(r['full_test_transfer']['blocks']['pass_transfer'] for r in results)
    output=dict(status='PASS' if primary else 'FAIL',criterion='Every full test true length and known code must satisfy the block envelope; block sensitivities retained separately.',controls=results,regenerated_cipher_words=decodedwords,block_tests=sum(len(r['test_block_transfer']) for r in results),block_failures=sum(not s['pass_transfer'] for r in results for s in r['test_block_transfer']),full_source_sensitivity_pass=all(r['full_test_transfer']['full']['pass_transfer'] for r in results),independent_plaintext_works=2,note='Two conditions share Cicero; blocks share sources and are not independent. All known labels retained, including training-unseen ones.')
    dump('results/variable-frequency/controls.json',output)
    paths=['data/variable-frequency/envelope-freeze.json','results/variable-frequency/controls.json']
    dump('data/variable-frequency/control-freeze.json',dict(status='FROZEN_BEFORE_NATIVE_APPLICATION',hashes=snapshot(paths)))
    print('TRANSFER',output['status'],'blocks',output['block_tests'],'failures',output['block_failures'],'fullsource',output['full_source_sensitivity_pass'],flush=True)

def native():
    for p in ['input-freeze','envelope-freeze','control-freeze']:verify('data/variable-frequency/'+p+'.json')
    envs=read('data/variable-frequency/references.json')['envelopes'];cal=read('results/variable-frequency/controls.json');man=read('data/homophonic/manuscript.json');output={}
    for split in ['train','development','reserve']:
        words=[];excluded=0;loci=0
        for r in man['rows']:
            if r['split']!=split:continue
            loci+=1;w=[]
            for c in r['atoms']:
                if c==' ':
                    if w:
                        if '?' in w:excluded+=1
                        else:words.append(w)
                    w=[]
                else:w.append(c)
            assert not w
        cc=Counter(c for w in words for c in w)
        if split=='train':assert sum(cc.values())==23721
        output[split]=dict(loci=loci,known_words=len(words),excluded_unknown_words=excluded,cipher_units=sum(cc.values()),under_envelopes={kind:bound(cc,frac(env['a']),frac(env['b'])) for kind,env in envs.items()})
    result=dict(status='CONSTRAINT_NOT_A_DECIPHERMENT',transfer_status=cal['status'],block_transfer_failures=cal['block_failures'],native=output,scope='All used79H/16H/16P are development; no015key used, no new reserve consumed.',interpretation='Only compatibility with the stated reference envelopes. Neither rejection nor admission identifies language or supplies a reading.')
    dump('results/variable-frequency/native.json',result)
    for s,r in output.items():print('NATIVE',s,r['cipher_units'],[(k,x['lower']['value'],x['upper']['value'],x['rejects']) for k,x in r['under_envelopes'].items()],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['freeze','toys','references','controls','native']);a=p.parse_args();globals()[a.command]()
