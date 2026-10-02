"""Experiment013 score runner; v1 equality assertion compared int/string JSON keys.
Selection and all input hashes remain unchanged. No numerical method change.
"""
import json
from collections import Counter
from transcription_units import (FILES, MODES, SPLITS, read, sha, write, collect,
                                 corpus_result, screen)

def score():
    frozen=read('data/transcription-units/selection-freeze.json')
    for p,h in frozen['files'].items():assert sha(p)==h,p
    selection=read('data/transcription-units/selection.json');refs=read('results/language-screen/constraints.json');full={};paired={};coverage={}
    for mode in MODES:
        full[mode]={};paired[mode]={};coverage[mode]={}
        for split in SPLITS:
            full[mode][split]={};paired[mode][split]={};coverage[mode][split]={}
            subset=[r for r in selection['rows'] if r['split']==split]
            for s in FILES:
                fw=collect(selection,s,split,mode,False);pw=collect(selection,s,split,mode,True)
                f=corpus_result(fw,False);p=corpus_result(pw,False)
                f['language_screen']=screen(f['full'],refs);p['language_screen']=screen(p['full'],refs)
                full[mode][split][s]=f;paired[mode][split][s]=p
                eligible=[r['sources'][s] for r in subset if r['sources'][s]['eligible']]
                coverage[mode][split][s]=dict(scope_loci=len(subset),eligible_loci=len(eligible),missing_or_ineligible_loci=len(subset)-len(eligible),raw_words=sum(len(r['words'][mode]) for r in eligible),known_words=len(fw),unknown_words=sum(not w['known'] for r in eligible for w in r['words'][mode]),paired_known_words=len(pw),paired_word_fraction=len(pw)/len(fw),paired_unit_fraction=p['full']['units']/f['full']['units'])
    for split in SPLITS:
        assert json.loads(json.dumps(full['primary'][split]['GC']['full']))==refs['native'][split]['full']
    triplets={};ratios={}
    for mode in MODES:
        triplets[mode]={};ratios[mode]={}
        for split in SPLITS:
            ps=[p for p in selection['pairs'][mode] if p['split']==split]
            ct=Counter(tuple(''.join(p['words'][s]['atoms']) for s in FILES) for p in ps)
            triplets[mode][split]=[dict(GC=t[0],ZL=t[1],IT=t[2],tokens=n) for t,n in sorted(ct.items(),key=lambda v:(-v[1],v[0]))]
            ratios[mode][split]={}
            for s in ['ZL','IT']:
                base=paired[mode][split]['GC']['full'];dest=paired[mode][split][s]['full']
                ratios[mode][split][s]=dict(total_unit_ratio=dest['units']/base['units'],words=len(ps),shorter=sum(len(p['words'][s]['atoms'])<len(p['words']['GC']['atoms']) for p in ps),equal=sum(len(p['words'][s]['atoms'])==len(p['words']['GC']['atoms']) for p in ps),longer=sum(len(p['words'][s]['atoms'])>len(p['words']['GC']['atoms']) for p in ps))
    gate=full['primary']['train']['ZL']['language_screen'];passing=[la for la,m in gate.items() if m['primary']['screen_pass']]
    out=dict(status='TRANSCRIPTION_UNIT_SCREEN_NOT_A_DECIPHERMENT',selection_freeze_sha256=sha('data/transcription-units/selection-freeze.json'),primary_gate='full primary-space training ZL, reference primary apostrophe condition',primary_passing_languages=passing,reference_sha256=sha('results/language-screen/constraints.json'),full=full,paired=paired,coverage=coverage,paired_ratios=ratios,paired_triplets=triplets,limitations='Positional matches with equal source word counts can retain compensating split/join differences; not a verified character mapping. Three transcriptions concern one manuscript.')
    write('results/transcription-units/profiles.json',out)
    summary=dict(primary_passing_languages=passing,full_train={s:{k:full['primary']['train'][s]['full'][k] for k in ['units','words','mean_word_length','max_frequency','most_frequent','types']} for s in FILES},paired_train={s:{k:paired['primary']['train'][s]['full'][k] for k in ['units','words','mean_word_length','max_frequency','most_frequent','types']} for s in FILES},ratios=ratios['primary']['train'],coverage=coverage['primary']['train'])
    write('results/transcription-units/summary.json',summary);print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=="__main__":score()
