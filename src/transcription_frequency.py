"""Experiment017: frozen016frequency bounds on frozen013transcriptions."""
from collections import Counter
import argparse
import variable_frequency as v

ROOT=v.ROOT;SOURCES=['GC','ZL','IT'];SPLITS=['train','development','reserve'];MODES=['primary','comma_join'];SCOPES=['full','paired']

def freeze():
    assert not (ROOT/'data/transcription-frequency/input-freeze.json').exists()
    for p in ['data/variable-frequency/input-freeze.json','data/variable-frequency/envelope-freeze.json','data/variable-frequency/control-freeze.json','data/transcription-units/selection-freeze.json']:v.verify(p)
    paths=['notes/experiment-017-design.md','notes/experiment-016-report-correction.md','src/transcription_frequency.py','src/variable_frequency.py','data/variable-frequency/input-freeze.json','data/variable-frequency/envelope-freeze.json','data/variable-frequency/control-freeze.json','data/variable-frequency/references.json','results/variable-frequency/controls.json','results/variable-frequency/native.json','results/audit-012-frequency.json','data/transcription-units/selection-freeze.json','data/transcription-units/selection.json','data/raw/GC2a-n.txt','data/raw/ZL3b-n.txt','data/raw/IT2a-n.txt']
    v.dump('data/transcription-frequency/input-freeze.json',dict(status='FROZEN_BEFORE_017BOUNDS',hashes=v.snapshot(paths)))
    print('FROZEN',len(paths),flush=True)

def score():
    v.verify('data/transcription-frequency/input-freeze.json')
    assert v.read('results/variable-frequency/controls.json')['status']=='PASS'
    selection=v.read('data/transcription-units/selection.json');envs=v.read('data/variable-frequency/references.json')['envelopes'];results=[]
    previous=v.read('results/variable-frequency/native.json')['native']
    for mode in MODES:
        for split in SPLITS:
            rows=[r for r in selection['rows'] if r['split']==split]
            for source in SOURCES:
                eligible=[r['sources'][source] for r in rows if r['sources'][source]['eligible']]
                full=[w['atoms'] for r in eligible for w in r['words'][mode] if w['known']]
                paired=[p['words'][source]['atoms'] for p in selection['pairs'][mode] if p['split']==split]
                coverage=dict(selected_loci=len(rows),eligible_loci=len(eligible),missing_or_ineligible_loci=len(rows)-len(eligible),source_word_slots=sum(len(r['words'][mode]) for r in eligible),unreadable_words=sum(not w['known'] for r in eligible for w in r['words'][mode]),full_known_words=len(full),paired_known_words=len(paired),paired_word_fraction=len(paired)/len(full),full_known_units=sum(map(len,full)),paired_known_units=sum(map(len,paired)))
                for scope,words in [('full',full),('paired',paired)]:
                    counts=Counter(c for w in words for c in w)
                    bounds={kind:v.bound(counts,v.frac(env['a']),v.frac(env['b'])) for kind,env in envs.items()}
                    if mode=='primary' and source=='GC' and scope=='full':
                        assert bounds==previous[split]['under_envelopes'];assert len(words)==previous[split]['known_words']
                    results.append(dict(source=source,split=split,spacing=mode,scope=scope,words=len(words),coverage=coverage,bounds=bounds))
    primary=next(r for r in results if (r['source'],r['split'],r['spacing'],r['scope'])==('ZL','train','primary','full'))
    rejected=primary['bounds']['blocks']['rejects'];decision='INCOMPATIBLE_WITH_EMPIRICAL_ENVELOPE' if rejected else 'ADMITTED_NECESSARY_CONSTRAINT_ONLY'
    out=dict(status='COMPLETED_CONSTRAINT_NOT_A_DECIPHERMENT',primary_condition=dict(source='ZL',split='train',spacing='primary',scope='full',envelope='blocks'),primary_decision=decision,results=results,total_bounds=sum(len(r['bounds']) for r in results),primary_training_comparison=[r for r in results if r['split']=='train' and r['spacing']=='primary' and r['scope']=='full'],next_action='Do not run this bounded ZL Latin decoder under the current envelope assumptions.' if rejected else 'Freeze operational015solver reuse protocol, then actually fit ZL. Admission is not a reading.',limits='All source texts depict one manuscript; strict paired slots do not prove physical alignment. Existing016known controls are reused. No new reserve consumed. Not exclusion of all Latin.')
    v.dump('results/transcription-frequency/bounds.json',out)
    for r in out['primary_training_comparison']:
        b=r['bounds']['blocks'];print('PRIMARY COMPARISON',r['source'],b['cipher_units'],b['lower']['value'],b['upper']['value'],'rejects',b['rejects'],flush=True)
    print('DECISION',decision,'all bounds',out['total_bounds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['freeze','score']);a=p.parse_args();globals()[a.command]()
