"""Independent016 audit, using prefix fixed points instead of iterative pruning."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from itertools import product
import hashlib,json,unicodedata as U

ROOT=Path(__file__).resolve().parents[1];ALPH='abcdefghiklmnopqrstuxyz'
def read(p):return json.loads((ROOT/p).read_text())
def q(p):return F(p['numerator'],p['denominator'])
def checkq(saved,value):assert q(saved)==value and saved['value']==float(value),(saved,value)
def hidden(s):return any(U.category(c)=='Co' or c.isdigit() or c in {'~','&','⁊'} for c in s)
def ranges(ws):
    result=[];s=0;n=0
    for i,w in enumerate(ws):
        n+=len(w)
        if n>=5000:result.append((s,i+1));s=i+1;n=0
    return result,dict(start_word=s,words=len(ws)-s,units=n)
def measure(ws,stored=None):
    singles=Counter();pairs=Counter();obscure=Counter();uncertainpairs=0
    for word in ws:
        for c in word:
            singles[c]+=1
            if hidden(c):obscure[c]+=1
        for a,b in zip(word,word[1:]):
            pairs[a,b]+=1
            if hidden(a) or hidden(b):uncertainpairs+=1
    N=sum(singles.values());clear={s:n for s,n in singles.items() if s not in obscure};clear2={p:n for p,n in pairs.items() if not any(hidden(c) for c in p)}
    C=max(clear.values(),default=0);P=max(clear2.values(),default=0);a=F(C+sum(obscure.values()),N);b0=F(P+uncertainpairs,N);b=min(a,b0)
    if stored is not None:
        assert stored['units']==N and stored['words']==len(ws) and stored['unit_counts']==singles
        assert Counter({tuple(r['pair']):r['count'] for r in stored['pair_counts']})==pairs
        assert stored['opaque_counts']==obscure and stored['opaque_units']==sum(obscure.values()) and stored['opaque_pair_positions']==uncertainpairs
        assert stored['max_clear_letter_count']==C and stored['max_clear_pair_count']==P
        assert stored['letter_maximizers']==sorted(c for c,n in clear.items() if n==C)
        assert stored['pair_maximizers']==[list(p) for p,n in sorted(clear2.items()) if n==P]
        checkq(stored['a_upper'],a);checkq(stored['b_raw_upper'],b0);checkq(stored['b_upper'],b)
    return a,b
def endpoints(cc,a,b):
    N=sum(cc.values());freq=sorted((F(n,N) for n in cc.values()),reverse=True);cum=F(0);stable=[]
    for k in range(len(freq)+1):
        t=2-cum
        if sum(f>t*b for f in freq)==k:stable.append(t)
        if k<len(freq):cum+=freq[k]
    assert stable;hi=max(stable);lo=max(F(1),freq[0]/a)
    if lo<=hi:
        lens=[1 if f>hi*b else 2 for f in freq]
        assert sum(f*l for f,l in zip(freq,lens))==hi
        assert all(f<=hi*(a if l==1 else b) for f,l in zip(freq,lens))
    return lo,hi
def inspect_bound(cc,a,b,r):
    N=sum(cc.values());lo,hi=endpoints(cc,a,b)
    assert r['cipher_units']==N and r['cipher_types']==len(cc) and r['cipher_counts']==cc
    checkq(r['a'],a);checkq(r['b'],b);checkq(r['lower'],lo);checkq(r['upper'],hi);checkq(r['max_cipher_frequency'],F(max(cc.values()),N))
    assert r['maximizers']==sorted(c for c,n in cc.items() if n==max(cc.values())) and r['rejects']==(lo>hi)
    t=F(2)
    for step in r['iterations']:
        checkq(step['upper_before'],t);forced=sorted(c for c,n in cc.items() if F(n,N)>t*b);mass=F(sum(cc[c] for c in forced),N);nt=min(t,2-mass)
        assert forced==step['forced_monographs'];checkq(step['forced_mass'],mass);checkq(step['upper_after'],nt);t=nt
    assert t==hi

def main():
    hashes=0
    for p in ['data/variable-frequency/input-freeze.json','data/variable-frequency/envelope-freeze.json','data/variable-frequency/control-freeze.json','data/language-screen/source-freeze.json','data/variable-units-v1/calibration-freeze.json','results/variable-units-v1/solver-freeze.json']:
        ob=read(p)
        for f,h in ob.get('hashes',ob.get('files')).items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h,f;hashes+=1
    ref=read('data/variable-frequency/references.json');sources={'caesar':[list(w) for n in range(1,8) for w in (ROOT/f'data/latin/caesar-gallic-{n}.txt').read_text().split()]}
    for name in ['macer','aquinas']:sources[name]=[w for r in read(f'data/language-screen/{name}.json')['records'] for w in r['words']]
    profiles=0;collected={'blocks':[],'full':[]}
    for name,ws in sources.items():
        ob=ref['sources'][name];a,b=measure(ws,ob['full']);profiles+=1;collected['full'].append((name,'full',a,b));bs,tail=ranges(ws);assert tail==ob['tail'] and len(bs)==len(ob['blocks'])
        for i,((s,e),saved) in enumerate(zip(bs,ob['blocks'])):
            assert saved['index']==i and saved['start_word']==s and saved['end_word']==e;a,b=measure(ws[s:e],saved['profile']);profiles+=1;collected['blocks'].append((name,str(i),a,b))
    for kind,entries in collected.items():
        env=ref['envelopes'][kind];a=max(x[2] for x in entries);b=max(x[3] for x in entries);checkq(env['a'],a);checkq(env['b'],b);assert env['entries']==len(entries)
        assert env['a_witnesses']==[dict(source=n,scope=i) for n,i,x,y in entries if x==a]
        assert env['b_witnesses']==[dict(source=n,scope=i) for n,i,x,y in entries if y==b]
    toy=read('results/variable-frequency/toys.json');cases=assigns=admitted=0
    for G in range(1,5):
        for cc in product(range(1,4),repeat=G):
            N=sum(cc)
            for ia in range(1,13):
                for ib in range(ia+1):
                    lo,hi=endpoints(dict(enumerate(cc)),F(ia,12),F(ib,12));cases+=1;feasible=[]
                    for mask in range(2**G):
                        ls=[1+((mask>>i)&1) for i in range(G)];m=sum(n*l for n,l in zip(cc,ls));assigns+=1
                        if all(n*12<=m*(ia if l==1 else ib) for n,l in zip(cc,ls)):feasible.append(F(m,N))
                    assert bool(feasible)==(lo<=hi)
                    if feasible:admitted+=1;assert max(feasible)==hi and min(feasible)>=lo
    assert cases==toy['length_distribution_cases'] and assigns==toy['length_assignments'] and admitted==toy['cases_with_feasible_assignment'] and toy['false_admissions']==0
    decodedcases=0
    for book in [['a','b'],['a','b','aa','ab','ba','bb'],['a','b','aa','ab','ba','bb','a'],['a','b','ab']]:
        for length in range(1,5):
            for codes in product(range(len(book)),repeat=length):
                word=list(''.join(book[c] for c in codes));a,b=measure([word]);lo,hi=endpoints(Counter(codes),a,b);assert lo<=F(len(word),length)<=hi;decodedcases+=1
    assert decodedcases==toy['exact_decoded_cipher_cases']
    ctl=read('results/variable-frequency/controls.json');bounds=0;decodedwords=0;transferflags=[];blockflags=[];fullflags=[]
    def check(rows,key,env,saved):
        nonlocal bounds
        cc=Counter(c for r in rows for c in r['cipher']);P=sum(len(r['plain']) for r in rows);N=sum(cc.values());a,b=q(env['a']),q(env['b']);lam=F(P,N)
        inspect_bound(cc,a,b,saved['bound']);bounds+=1;checkq(saved['true_expansion'],lam);assert saved['words']==len(rows) and saved['plaintext_units']==P
        lo,hi=endpoints(cc,a,b);lp=lo<=lam<=hi;fail=[]
        for c,n in sorted(cc.items()):
            limit=a if len(key[c])==1 else b
            if F(n,P)>limit:fail.append((c,key[c],n,F(n,P),limit))
        assert len(fail)==len(saved['known_code_violations'])
        for (c,u,n,f,m),r in zip(fail,saved['known_code_violations']):
            assert (r['code'],r['output'],r['cipher_count'])==(c,u,n);checkq(r['required_frequency'],f);checkq(r['allowed_frequency'],m)
        assert saved['true_expansion_admitted']==lp and saved['pass_transfer']==(lp and not fail)
        return lp and not fail
    for rec in ctl['controls']:
        old=read(f"data/variable-units-v1/{rec['name']}.json");plain=read('data/latin/control-'+rec['source']+'.json');key={str(c):u for u,cs in old['groups'].items() for c in cs};rows={}
        for split in ['train','test']:
            ww=[list(w) for row in plain['plain_'+split] for w in ''.join(' ' if c==23 else ALPH[c] for c in row).split()];assert len(ww)==len(old[split]);rows[split]=[]
            for w,r in zip(ww,old[split]):
                cipher=list(map(str,r['labels']));assert ''.join(w)==r['plain']==''.join(key[c] for c in cipher);decodedwords+=1;rows[split].append(dict(plain=w,cipher=cipher))
            a,b=measure(ww);env=dict(a=dict(numerator=a.numerator,denominator=a.denominator),b=dict(numerator=b.numerator,denominator=b.denominator));assert check(rows[split],key,env,rec['paired'][split])
        transferflags.append(check(rows['test'],key,ref['envelopes']['blocks'],rec['full_test_transfer']['blocks']));fullflags.append(check(rows['test'],key,ref['envelopes']['full'],rec['full_test_transfer']['full']))
        bs,tail=ranges([r['plain'] for r in rows['test']]);assert tail==rec['test_block_tail'] and len(bs)==len(rec['test_block_transfer'])
        for i,((s,e),r) in enumerate(zip(bs,rec['test_block_transfer'])):
            assert (r['index'],r['start_word'],r['end_word'])==(i,s,e);blockflags.append(check(rows['test'][s:e],key,ref['envelopes']['blocks'],r))
    assert ctl['status']==('PASS' if all(transferflags) else 'FAIL') and ctl['full_source_sensitivity_pass']==all(fullflags)
    assert ctl['regenerated_cipher_words']==decodedwords and ctl['block_tests']==len(blockflags) and ctl['block_failures']==sum(not f for f in blockflags)
    native=read('results/variable-frequency/native.json');man=read('data/homophonic/manuscript.json');loci=0
    for split,saved in native['native'].items():
        ws=[];excluded=0;ls=0
        for r in man['rows']:
            if r['split']!=split:continue
            word=[];ls+=1
            for c in r['atoms']:
                if c!=' ':word.append(c)
                else:
                    if word:
                        if '?' in word:excluded+=1
                        else:ws.append(word)
                    word=[]
            assert not word
        counts=Counter(c for w in ws for c in w);assert saved['loci']==ls and saved['known_words']==len(ws) and saved['excluded_unknown_words']==excluded and saved['cipher_units']==sum(counts.values());loci+=ls
        for kind,r in saved['under_envelopes'].items():env=ref['envelopes'][kind];inspect_bound(counts,q(env['a']),q(env['b']),r);bounds+=1
    assert native['transfer_status']==ctl['status'] and native['block_transfer_failures']==ctl['block_failures']
    result=dict(status='PASS',scope='Exact necessary constraints relative to stated source envelopes, not a language theorem or decipherment.',frozen_hash_checks=hashes,reconstructed_source_profiles=profiles,exhaustive_length_cases=cases,exhaustive_length_assignments=assigns,endpoint_witnesses=admitted,decoded_toy_cases=decodedcases,verified_real_cipher_words=decodedwords,independently_recomputed_bounds=bounds,native_loci=loci,primary_transfer=ctl['status'],block_transfer_failures=ctl['block_failures'])
    (ROOT/'results/audit-012-frequency.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
