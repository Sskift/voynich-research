"""Independent all-skeleton020 v2 infinite full/OOV completion-mass audit.

Run only after all eighteen v2 cases and their terminal result freeze exist.
The graph/linear-solve/prefix-exit implementation is the independent v1 audit,
unchanged; no production model or numerical implementation is imported.
"""
from pathlib import Path
import argparse
import datetime
import itertools
import json
import math
import sys
import traceback

from audit_omission_mass import GraphMass, log_add, omit
from audit_omission_recovery import IndependentModel, Integrity


ROOT = Path(__file__).resolve().parents[1]
DATA = 'data/omission-oracle'
RESULTS = 'results/omission-oracle-v2'
FREEZE = 'data/omission-oracle-v2/input-freeze.json'
LANGUAGES = ('la','de','it')
MECHANISMS = ('ordered','iid','copy')
VOWELS = frozenset('aeiou')


def need(condition,message):
    if not condition:
        raise AssertionError(message)


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def sha(path):
    import hashlib
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def mass_comparison(independent,archived,ordinary,where):
    if independent is None:
        need(archived is None and ordinary==0.,where+': absent mass mismatch')
        return dict(independent_log_mass=None,archived_log_mass=None,independent_mass=0.,
                    absolute_mass_error=0.,relative_mass_error=0.,absolute_log_error=0.,underflow=False)
    need(archived is not None and math.isfinite(archived),where+': missing finite log mass')
    error=independent-archived
    need(abs(error)<=1e-8,f'{where}: log disagreement {error}')
    observed=math.exp(independent)
    need(ordinary==math.exp(archived),where+': archived log/ordinary mismatch')
    return dict(independent_log_mass=independent,archived_log_mass=archived,
                independent_mass=observed,absolute_mass_error=abs(observed-ordinary),
                relative_mass_error=abs(math.expm1(error)),absolute_log_error=abs(error),
                underflow=observed==0.)


def audit_skeleton(model,skeleton,archived):
    need(archived['skeleton']==skeleton and archived['status']=='COMPLETE'
         and archived['numerical_complete'] is True and archived['search_complete'] is True,
         'Incomplete/mismatched archived certificate: '+skeleton)
    need(archived['numerical_certificate']['numerical_complete'] is True,
         'Numerical certificate completeness mismatch')
    lexical=sorted(model.lexical_images.get(skeleton,()))
    need(archived['lexicon_candidates']==lexical and
         set(archived['lexicon_log_char_probabilities'])==set(lexical), 'Lexical inventory mismatch: '+skeleton)
    singleton=all(c in VOWELS for c in skeleton)
    if singleton:
        need(lexical==([skeleton] if skeleton in model.words else []), 'All-vowel preimage is not singleton')
        full=model.char_log(skeleton)
        oov=None if lexical else full
        diagnostics=[]
        graph_states=graph_edges=0
        accounting=dict(kind='All-vowel singleton; no graph or infinite sum required',lexical_words=len(lexical))
        need([item['word'] for item in archived['top_oov']]==([] if lexical else [skeleton]),
             'All-vowel OOV singleton set')
    else:
        need(not any(c in VOWELS for c in skeleton),'Invalid mixed omission observation')
        graph=GraphMass(model,skeleton)
        full=graph.full_log_mass()
        oov,accounting=graph.oov_log_mass(lexical)
        diagnostics=graph.position_diagnostics
        graph_states=len(graph.states)
        graph_edges=sum(map(len,graph.edges.values()))
        need(len(archived['top_oov'])==2,'Nontrivial archived OOV candidate set has wrong size')
    lexical_logs={word:model.char_log(word) for word in lexical}
    terms=list(lexical_logs.values())+([] if oov is None else [oov])
    partition_error=abs(log_add(terms)-full)
    need(partition_error<=1e-10,'Independent full/OOV/lexicon partition mismatch: '+skeleton)
    word_checks=[]
    for word,independent in lexical_logs.items():
        saved=archived['lexicon_log_char_probabilities'][word]
        error=abs(independent-saved)
        need(error<=1e-10,'Lexical G(w) mismatch: '+word)
        word_checks.append(dict(kind='lexical',word=word,independent_log_probability=independent,
                                archived_log_probability=saved,absolute_log_error=error))
    seen=set()
    for item in archived['top_oov']:
        word=item['word']
        need(word not in seen and word not in model.words and omit(word)==skeleton,'OOV candidate invalid: '+word)
        seen.add(word)
        independent=model.char_log(word)
        error=abs(independent-item['log_char_probability'])
        need(error<=1e-10,'OOV G(w) mismatch: '+word)
        word_checks.append(dict(kind='archived_top_oov',word=word,independent_log_probability=independent,
                                archived_log_probability=item['log_char_probability'],absolute_log_error=error))
    return dict(singleton=singleton,full_mass=mass_comparison(full,archived['log_char_mass'],archived['char_mass'],skeleton+'.full'),
                oov_mass=mass_comparison(oov,archived['log_oov_char_mass'],archived['oov_char_mass'],skeleton+'.oov'),
                independent_partition_log_error=partition_error,lexical_candidates=lexical,
                concrete_word_probabilities=word_checks,reachable_graph_states=graph_states,
                reachable_graph_edges=graph_edges,position_equation_diagnostics=diagnostics,oov_accounting=accounting)


def audit():
    need((ROOT/f'{RESULTS}/summary.json').exists() and (ROOT/f'{RESULTS}/result-freeze.json').exists(),
         'Await complete v2 run and result freeze')
    need(not (ROOT/f'{RESULTS}/computation-failure.json').exists(),'v2 run is incomplete')
    integrity=Integrity()
    original=integrity.verify(f'{DATA}/input-freeze.json')
    prepared=integrity.verify(f'{DATA}/prepared-manifest.json')
    v2=integrity.verify(FREEZE)
    terminal=integrity.verify(f'{RESULTS}/result-freeze.json')
    original_hash=sha(f'{DATA}/input-freeze.json')
    need(prepared['input_hashes']==original['hashes'] and prepared['input_freeze_sha256']==original_hash,
         'Original prepared chain mismatch')
    need(v2['parent020_input_hashes']==original['hashes'] and v2['parent020_input_freeze_sha256']==original_hash,
         'v2 parent hash chain mismatch')
    need(FREEZE in terminal['hashes'] and terminal['hashes'][FREEZE]==sha(FREEZE),'v2 result freeze chain')
    need((set(original['hashes'])|set(prepared['hashes']))<=set(v2['hashes']), 'v2 missing old inputs/prepared artifacts')
    all_cases=[f'{l}/{m}-{rule}' for l,m,rule in itertools.product(LANGUAGES,MECHANISMS,('identity','omit-aeiou'))]
    expected_files={f'{RESULTS}/{kind}/{case}.json' for kind in ('predictions','scores') for case in all_cases}
    expected_files|={f'{RESULTS}/{l}-completion-certificates.json' for l in LANGUAGES}
    expected_files|={f'{RESULTS}/summary.json',f'{RESULTS}/result-freeze.json'}
    need({str(p.relative_to(ROOT)) for p in (ROOT/RESULTS).rglob('*.json')}==expected_files,'Incomplete v2 output inventory')
    need(expected_files-{f'{RESULTS}/result-freeze.json'}<=set(terminal['hashes']),'Unfrozen v2 result file')
    summary=read(f'{RESULTS}/summary.json')
    need(summary['status']=='COMPLETE_ORACLE_CALIBRATION_NO_NATIVE_SCORES','v2 completion marker')
    result={}; total_unique=total_singletons=total_word_checks=0
    maxima=dict(full_log_error=0.,oov_log_error=0.,full_relative_error=0.,oov_relative_error=0.,
                partition_log_error=0.,concrete_word_log_error=0.,scaled_equation_residual=0.,
                componentwise_relative_equation_residual=0.)
    for language in LANGUAGES:
        model=IndependentModel(language,integrity)
        certificate_path=f'{RESULTS}/{language}-completion-certificates.json'
        certificates=read(certificate_path)
        omission={key.split('|',1)[1]:value for key,value in certificates.items() if key.startswith('omit-aeiou|')}
        observed=set()
        for mechanism in MECHANISMS:
            inputs=read(f'{DATA}/tests/{language}/{mechanism}-omit-aeiou.json')
            need(inputs['rule']=='omit-aeiou' and inputs['language']==language,'Observation metadata mismatch')
            observed.update(''.join(t['units']) for r in inputs['rows'] for t in r['tokens'] if t['units'] is not None)
        need(set(omission)==observed,'Archived/observed omission skeleton inventory mismatch')
        audited={}
        for number,skeleton in enumerate(sorted(omission),1):
            current=audit_skeleton(model,skeleton,omission[skeleton])
            audited[skeleton]=current
            total_unique+=1; total_singletons+=int(current['singleton'])
            total_word_checks+=len(current['concrete_word_probabilities'])
            maxima['full_log_error']=max(maxima['full_log_error'],current['full_mass']['absolute_log_error'])
            maxima['oov_log_error']=max(maxima['oov_log_error'],current['oov_mass']['absolute_log_error'])
            maxima['full_relative_error']=max(maxima['full_relative_error'],current['full_mass']['relative_mass_error'])
            maxima['oov_relative_error']=max(maxima['oov_relative_error'],current['oov_mass']['relative_mass_error'])
            maxima['partition_log_error']=max(maxima['partition_log_error'],current['independent_partition_log_error'])
            for word in current['concrete_word_probabilities']:
                maxima['concrete_word_log_error']=max(maxima['concrete_word_log_error'],word['absolute_log_error'])
            for position in current['position_equation_diagnostics']:
                maxima['scaled_equation_residual']=max(maxima['scaled_equation_residual'],position['maximum_scaled_all_equation_residual'])
                maxima['componentwise_relative_equation_residual']=max(maxima['componentwise_relative_equation_residual'],position['maximum_componentwise_relative_equation_residual'])
            if number%100==0:
                print(json.dumps(dict(progress='INDEPENDENT_MASS_AUDIT',language=language,complete=number,total=len(omission))),flush=True)
        result[language]=dict(unique_omission_skeletons=len(audited),all_vowel_singletons=sum(x['singleton'] for x in audited.values()),
                              completion_certificates_sha256=sha(certificate_path),skeletons=audited)
        print(json.dumps(dict(progress='LANGUAGE_MASS_AUDIT_COMPLETE',language=language,unique_omission_skeletons=len(audited))),flush=True)
    return dict(status='V2_ALL_UNIQUE_OMISSION_MASS_AUDIT_PASS_NOT_A_READING',
                languages=list(LANGUAGES),unique_omission_skeletons=total_unique,all_vowel_singletons=total_singletons,
                concrete_word_probability_checks=total_word_checks,maximum_errors=maxima,results=result,
                scope=dict(independent_algorithm='Existing independent v1 explicit graph; per-position numpy.linalg.solve; iterative noncore fill; nonnegative forward prefix-exit OOV sum.',
                    character_word_probabilities='Independently recomputed from frozen018 tables normalized in memory.',
                    infinite_mass='Every unique omission skeleton in all three full completion-certificate archives; no word-length truncation.',
                    residuals='Represented scaled equation residuals only; not interval or rigorous forward-error certificates.',
                    top2_search='Only archived concrete word probabilities and membership checked; global top2 search completeness is not claimed.',
                    production_modules_imported=False,production_v2_numerical_implementation_used=False,
                    manuscript_scoring=False,calibration_claim=False),
                tolerances=dict(mass_log_agreement=1e-8,partition_log_agreement=1e-10,concrete_word_log_agreement=1e-10),
                hash_manifests=integrity.manifests,distinct_frozen_files=len(integrity.files),
                v2_input_freeze_sha256=sha(FREEZE),parent020_input_freeze_sha256=original_hash,
                result_freeze_sha256=sha(f'{RESULTS}/result-freeze.json'),audit_source_sha256=sha('src/audit_omission_mass_v2.py'),
                independent_helpers_sha256={path:sha(path) for path in ('src/audit_omission_mass.py','src/audit_omission_recovery.py')})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='results/audit-016-omission-mass-v2.json')
    args=parser.parse_args()
    output=ROOT/args.output
    need(not output.exists(),'Preserve existing audit; choose another --output')
    try:
        report=audit()
    except Exception as error:
        report=dict(status='V2_MASS_AUDIT_FAILED',error_type=type(error).__name__,error=str(error),
                    traceback=traceback.format_exc(),audit_source_sha256=sha('src/audit_omission_mass_v2.py'))
        stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        failed=output.with_name(output.stem+'.failure-'+stamp+'.json')
        failed.parent.mkdir(parents=True,exist_ok=True)
        with failed.open('x',encoding='utf-8') as handle:json.dump(report,handle,ensure_ascii=False,indent=2)
        print(json.dumps(report,ensure_ascii=False,indent=2),file=sys.stderr)
        raise
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8') as handle:
        json.dump(report,handle,ensure_ascii=False,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({key:report[key] for key in ('status','unique_omission_skeletons','all_vowel_singletons',
                     'concrete_word_probability_checks','maximum_errors')},indent=2),flush=True)


if __name__=='__main__':
    main()
