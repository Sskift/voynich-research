"""Verify the frozen split and score using separate direct counting."""
from collections import Counter,defaultdict
from pathlib import Path
import hashlib,json,re
from corpus import read
from label_rules import transform,candidates

ROOT=Path(__file__).resolve().parents[1]


def main():
    fit_path=ROOT/'results/label-rule-fit.json'
    fit=json.loads(fit_path.read_text())
    ev=json.loads((ROOT/'results/label-rule-external-evaluation.json').read_text())
    assert hashlib.sha256(fit_path.read_bytes()).hexdigest()==ev['frozen_fit_sha256']
    for file,sha in fit['artifact_inputs'].items():assert hashlib.sha256((ROOT/file).read_bytes()).hexdigest()==sha
    proposals=json.loads((ROOT/'data/herbal-weak-image-pairs.json').read_text())
    extracted=set()
    for q in [15,19]:
        text=(ROOT/f'data/raw/q{q}-index.html').read_text()
        extracted|={(int(m[1]),m[2]) for m in re.finditer(r'Fragment (\d+) shows some similarity with the plant on\s*<A HREF="[^"]+">(f[^<]+)</A>',text)}
    assert extracted|{(150,'f55v')}=={(x['fragment'],x['target']) for x in proposals}
    assert len(proposals)==22 and sum(x['primary'] for x in proposals)==12
    inventory={x['locus']:x for x in json.loads((ROOT/'results/pharma-labels.json').read_text())}
    rows=[inventory[k] for k in fit['training_loci']]
    assert not({x['word'] for x in rows}&set(fit['excluded_word_types']))
    assert not(set(fit['fitting_body_pages'])&{x['target'] for x in proposals})
    rules=candidates([x['word'] for x in rows])
    assert len(rules)==fit['n_candidates']==670
    body=defaultdict(set)
    for x in read(ROOT/'data/raw/ZL3b-n.txt'):
        if x.page in fit['fitting_body_pages'] and x.code[1]=='P' and x.code[0]!='!':
            body[x.page].update(w for w in x.words if w)
    strata=defaultdict(list)
    for x in rows:strata[x['meta'].get('L'),x['meta'].get('H')].append(x)
    direct=[]
    for r in rules:
        score=0.
        for xs in strata.values():
            counts=Counter(x['page'] for x in xs)
            for x in xs:
                w=x['word'];v=transform(w,r)
                observed=int(v in body[x['page']])-int(w in body[x['page']])
                expected=sum(n*(int(v in body[p])-int(w in body[p])) for p,n in counts.items())/len(xs)
                score+=observed-expected
        direct.append(score/len(rows))
    maximum=max(direct)
    chosen=next(i for i,s in enumerate(direct) if s>=maximum-1e-12)
    assert rules[chosen]==fit['selected_rule']
    assert abs(direct[chosen]-fit['selected_centered_gain_per_label'])<1e-12
    for run in ev['runs']:
        text=defaultdict(set)
        for x in read(ROOT/'data/raw'/run['source'],run['spacing']):
            if x.code[1]=='P' and x.code[0]!='!':text[x.page].update(w for w in x.words if w)
        assert run['identity_hits']==sum(x['word'] in text[x['target']] for x in run['pairs'])==0
        assert run['rewritten_hits']==sum(transform(x['word'],fit['selected_rule']) in text[x['target']] for x in run['pairs'])==0
    profile=json.loads((ROOT/'results/label-position-profiles.json').read_text())
    for run in profile['runs']:
        for p in run['profiles']:
            assert len(p['matched'])==p['n_matched']
            assert p['n_matched']+len(p['excluded'])==131
            for x in p['matched']:
                assert x['pool_size']>0
                assert x['label_features']==[float(x['word'].startswith('o')),float(x['word'].startswith('qo')),float(x['word'].endswith('m'))]
    result=dict(status='PASS',frozen_fit_sha256_unchanged=True,training_input_hashes_unchanged=True,
                catalog_proposals_complete=22,primary_external_pairs=12,training_labels=len(rows),
                disjoint_training_word_types_and_herbal_pages=True,
                independently_scored_candidates=len(rules),selected_rule=rules[chosen]['id'],
                independently_recomputed_max_score=maximum,external_zero_hit_counts_confirmed=True,
                positional_profile_support_counts_checked=True)
    (ROOT/'results/audit-003.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
