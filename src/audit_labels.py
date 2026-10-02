"""Audit source integrity, label alignment, and independently count exact hits."""
from collections import Counter, defaultdict
from pathlib import Path
import hashlib,json
from corpus import read

ROOT=Path(__file__).resolve().parents[1]


def main():
    checked=0
    for name in ['sources.json','round-002-sources.json','pharma-image-sources.json','herbal-image-sources.json']:
        for item in json.loads((ROOT/'data'/name).read_text()):
            data=(ROOT/item['file']).read_bytes()
            assert len(data)==item['bytes'],item['file']
            assert hashlib.sha256(data).hexdigest()==item['sha256'],item['file']
            checked+=1
    labels=json.loads((ROOT/'results/pharma-labels.json').read_text())
    primary=[x for x in labels if x['agreed_single_word'] and x['subtype']=='Lf']
    assert len(labels)==234 and len(primary)==131
    assert len({x['locus'] for x in labels})==234
    for x in labels:
        assert all((ROOT/p).is_file() for p in x['images']),x['locus']
        if x['agreed_single_word']:
            assert x['zl_words']==x['it_words']==[x['word']]
    stats=json.loads((ROOT/'results/pharma-label-page-association.json').read_text())
    herbal=json.loads((ROOT/'results/herbal-label-association.json').read_text())
    reports=[]
    for saved, hs in zip(stats['runs'],herbal['runs']):
        loci=read(ROOT/'data/raw'/saved['source'],saved['spacing'])
        body=defaultdict(set)
        for x in loci:
            if x.code[1]=='P' and x.code[0]!='!':
                body[x.page].update(w for w in x.words if w)
        stratum_words=defaultdict(Counter)
        def stratum(x): return x['subtype'],x['meta'].get('L','?'),x['meta'].get('H','?')
        for x in primary:
            stratum_words[stratum(x)][x['word']]+=1
        observed=sum(x['word'] in body[x['page']] for x in primary)
        expected=0.
        for x in primary:
            frequencies=stratum_words[stratum(x)]
            expected+=sum(n for w,n in frequencies.items() if w in body[x['page']])/sum(frequencies.values())
        test=saved['tests']['Lf']
        assert observed==test['observed_present']
        assert abs(expected-test['exact_null_expected_present'])<1e-10
        assert abs(expected-test['null_mc_mean'])<.1
        assert observed==len(test['matches'])
        assert hs['source']==saved['source'] and hs['spacing']==saved['spacing']
        assert hs['observed_exact_present']==sum(x['word'] in body[x['target']] for x in hs['pairs'])==0
        assert hs['permutations']==24
        assert hs['posthoc_currier_scribe_sensitivity']['permutations']==6
        reports.append(dict(source=saved['source'],spacing=saved['spacing'],
                            observed=observed,independent_expected=expected,
                            herbal_hits=hs['observed_exact_present']))
    result=dict(status='PASS',source_entries_hash_checked=checked,
        label_positions=234,primary_plant_labels=131,image_paths_exist=True,
        independent_exact_count_checks=reports,
        limitation='Automated checks verify hashes, counts and expectations; image identities require the documented manual review.')
    (ROOT/'results/audit-002.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
