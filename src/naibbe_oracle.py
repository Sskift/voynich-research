"""Experiment021 staged runner: freeze -> prepare -> predict -> score.

Every output is exclusive-create. No manuscript paths or unknown-key search.
Prediction functions never parse truth, generator traces, or source texts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys

from naibbe_channel import normalize, load_table, encode_line, emission_model
from naibbe_decode import Decoder, train_model, edit_distance, SEGMENT_LENGTH
from naibbe_inverse import CodebookInverse, build_base_relation, SUPPORT_SCOPE
from naibbe_sources import extract_sources, INPUTS as SOURCE_INPUTS

ROOT = Path(__file__).resolve().parents[1]
DATA = 'data/naibbe-oracle'
OUT = 'results/naibbe-oracle'
AUTHOR = 'data/naibbe-review/author-f2675ec5'
TABLE = AUTHOR+'/naibbe_tables.csv'
LANGUAGES = ('la','it')
MODES = ('preserved','deleted')
RUN_CONTEXT = {}
INPUTS = (
    'notes/experiment-021-design.md','src/naibbe_channel.py','src/naibbe_decode.py',
    'src/naibbe_inverse.py','src/naibbe_sources.py','src/naibbe_oracle.py',
    *SOURCE_INPUTS,'data/raw/language-screen/dante-1012.txt',
    'data/raw/language-screen/boccaccio-21424.html','src/language_sources.py',
    'data/latin-sources.json',TABLE,AUTHOR+'/naibbe.py',AUTHOR+'/README.md',
    AUTHOR+'/LICENSE',AUTHOR+'/archive-inventory.json',
    'notes/naibbe-channel-review.md','notes/naibbe-and-units-primary-review.md',
)


def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def read(path): return json.loads((ROOT/path).read_text())
def write(path,obj):
    target = ROOT/path; target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('x') as handle:
        json.dump(obj,handle,ensure_ascii=False,indent=2,allow_nan=False);handle.write('\n')
def hashes(paths): return {path:sha(path) for path in sorted(set(paths))}
def verify(path):
    manifest = read(path)
    for name,digest in manifest['hashes'].items(): assert sha(name)==digest,(path,name)
    return manifest


def expected_cases():
    cases=[]
    for language in LANGUAGES:
        for variant in ('ordered-r0','ordered-r1','ordered-r2','shuffled-r0'):
            for mode in MODES:
                name=f'{language}-{variant}-{mode}'
                cases.append(dict(name=name,language=language,variant=variant,mode=mode,
                    test=f'{DATA}/tests/{name}.json',truth=f'{DATA}/truth/{name}.json',
                    model=f'{DATA}/{language}/model.json'))
    return cases


def verify_inputs():
    frozen=verify(DATA+'/input-freeze.json')
    assert frozen['status']=='FROZEN_BEFORE_SOURCE_PREPARATION_OR_SCORING'
    assert set(frozen['hashes'])==set(INPUTS) and frozen['config']==configuration()
    return frozen


def verify_prepared():
    verify_inputs()
    frozen=verify(DATA+'/prepared-freeze.json')
    assert frozen['status']=='PREPARED_NOT_PREDICTED'
    assert frozen['input_freeze_sha256']==sha(DATA+'/input-freeze.json')
    assert frozen['cases']==expected_cases()
    paths={DATA+'/config.json',DATA+'/source-provenance.json'}
    for case in expected_cases():
        paths.update((case['test'],case['truth'],case['model'],
                      f'{DATA}/generation/{case["language"]}-{case["variant"]}.json'))
    assert len(paths)==44 and set(frozen['hashes'])==paths
    assert read(DATA+'/config.json')==configuration()
    return frozen


def verify_prediction_freeze():
    frozen=verify(OUT+'/prediction-freeze.json')
    expected={f'{OUT}/predictions/{c["name"]}.json' for c in expected_cases()}
    assert len(expected)==16 and set(frozen['hashes'])==expected
    assert frozen['status']=='ALL_PREDICTIONS_SAVED_BEFORE_SCORING'
    assert frozen['prepared_freeze_sha256']==sha(DATA+'/prepared-freeze.json')
    assert frozen['input_freeze_sha256']==sha(DATA+'/input-freeze.json')
    return frozen


def validate_rows(rows,language):
    assert len(rows)==64
    assert [(r['id'],r['block']) for r in rows]==[
        (f'{language}-{i:03d}','A' if i<32 else 'B') for i in range(64)]


def configuration():
    return dict(experiment=21,channel='BASE52_AUTHOR_f2675ec5',source_chunk_length=128,
                test_chunks=64,blocks={'A':[0,32],'B':[32,64]},languages=list(LANGUAGES),
                modes=list(MODES),ordered_replicates=3,shuffle_replicates=1,
                ordered_seed='202102100+1000*language_index+replicate',
                shuffle_encoder_seed='202102190+1000*language_index',
                shuffle_source_seed='202102199+1000*language_index',
                known_key=True,matching_language_prior=True,score='APPROXIMATE_PATH_MAX_NOT_POSTERIOR',
                support_scope=SUPPORT_SCOPE,model_order=2,baseline_order=0,tau=20,additive=.5,
                single_probability='17/36 except final forced singleton',gap_drop=.03,
                max_dfa_states=100000,max_dfa_edges=1000000,max_viterbi_states=1000000,
                max_encoder_draws_per_chunk=1000000,
                primary_gate=dict(coverage=1.,character_edit_accuracy=.99,exact_chunk_accuracy=.8),
                manuscripts_read=False,original_word_boundaries_recovered=False)


def freeze():
    from naibbe_channel import self_check as channel_check
    from naibbe_decode import self_check as decode_check
    from naibbe_sources import self_test as source_check
    from naibbe_inverse import self_check as inverse_check
    checks = dict(channel=channel_check(),decoder=decode_check(),sources=source_check())
    inverse_check()
    checks['inverse'] = {'status':'PASS','pure_synthetic_checks':7222}
    write(DATA+'/input-freeze.json',dict(status='FROZEN_BEFORE_SOURCE_PREPARATION_OR_SCORING',
          config=configuration(),hashes=hashes(INPUTS),synthetic_checks=checks))
    print(json.dumps({'status':'INPUT_FROZEN','hashes':len(INPUTS)}),flush=True)


def prepare():
    verify_inputs()
    raw_sources = extract_sources(ROOT)
    _, table = load_table(ROOT/TABLE)
    outputs, test_cases, source_info = [], [], {}
    config_path = DATA+'/config.json'; write(config_path,configuration());outputs.append(config_path)
    for language_index,language in enumerate(LANGUAGES):
        RUN_CONTEXT.update(language=language)
        source = raw_sources[language]
        # Verify the previously selected Italian records still bind their raw
        # download; the new extractor uses text, never old normalized words.
        if language == 'it':
            for role in ('train','test'):
                record = read(source['provenance'][role]['path'])
                assert sha(record['raw_file']) == record['raw_sha256']
        train,test = normalize(source['train']),normalize(source['test'])
        assert len(test)>=64*SEGMENT_LENGTH
        chunks = [test[i*128:(i+1)*128] for i in range(64)]
        model_path = f'{DATA}/{language}/model.json'
        write(model_path,train_model(train));outputs.append(model_path)
        source_info[language] = dict(provenance=source['provenance'],
                normalized_train_sha256=hashlib.sha256(train.encode()).hexdigest(),
                normalized_test_sha256=hashlib.sha256(test.encode()).hexdigest(),
                normalized_train_letters=len(train),normalized_test_letters=len(test),
                test_selected_interval=[0,64*128],orthography='AUTHOR_NORMALIZATION_NO_UV_MERGE')
        variants = [(f'ordered-r{rep}',202102100+1000*language_index+rep,chunks) for rep in range(3)]
        shuffle_rng = random.Random(202102199+1000*language_index)
        shuffled = []
        for chunk in chunks:
            letters = list(chunk);shuffle_rng.shuffle(letters);shuffled.append(''.join(letters))
        variants.append(('shuffled-r0',202102190+1000*language_index,shuffled))
        for variant,seed,plain_chunks in variants:
            RUN_CONTEXT.update(variant=variant)
            rng = random.Random(seed)
            cases = {mode:dict(experiment=21,language=language,mode=mode,rows=[]) for mode in MODES}
            truths = {mode:dict(rows=[]) for mode in MODES}
            trace_rows = []
            for row_index,plain in enumerate(plain_chunks):
                RUN_CONTEXT.update(row_index=row_index)
                trace = encode_line(plain,table,rng)
                row_id = f'{language}-{row_index:03d}'
                block = 'A' if row_index<32 else 'B'
                trace_rows.append(dict(id=row_id,plaintext=plain,**trace))
                for mode in MODES:
                    observed = [u['cipher'] for u in trace['units']] if mode=='preserved' else trace['observed']
                    payloads = [u['plaintext'] for u in trace['units']] if mode=='preserved' else trace['payloads']
                    cases[mode]['rows'].append(dict(id=row_id,block=block,observations=observed))
                    truths[mode]['rows'].append(dict(id=row_id,block=block,plaintext=plain,payloads=payloads))
            trace_path = f'{DATA}/generation/{language}-{variant}.json'
            write(trace_path,dict(language=language,variant=variant,seed=seed,rows=trace_rows));outputs.append(trace_path)
            for mode in MODES:
                name = f'{language}-{variant}-{mode}'
                test_path,truth_path = f'{DATA}/tests/{name}.json',f'{DATA}/truth/{name}.json'
                write(test_path,cases[mode]);write(truth_path,truths[mode]);outputs.extend((test_path,truth_path))
                test_cases.append(dict(name=name,language=language,variant=variant,mode=mode,
                                       test=test_path,truth=truth_path,model=model_path))
            print(json.dumps(dict(prepared=language+'-'+variant,rows=len(plain_chunks))),flush=True)
    provenance_path=DATA+'/source-provenance.json';write(provenance_path,source_info);outputs.append(provenance_path)
    write(DATA+'/prepared-freeze.json',dict(status='PREPARED_NOT_PREDICTED',
          input_freeze_sha256=sha(DATA+'/input-freeze.json'),hashes=hashes(outputs),cases=test_cases))
    verify_prepared()
    print(json.dumps(dict(status='PREPARATION_COMPLETE',outputs=len(outputs),cases=len(test_cases))),flush=True)


def predict():
    prepared = verify_prepared()
    rows,table = load_table(ROOT/TABLE)
    relation = build_base_relation(rows)
    inverse = CodebookInverse(relation)
    emissions,emission_checks = emission_model(table)
    assert set(emissions)==set(relation)
    decoders = {}
    outputs = []
    for case in prepared['cases']:
        RUN_CONTEXT.update(case=case['name'])
        observations = read(case['test'])
        assert set(observations)=={'experiment','language','mode','rows'}
        assert observations['experiment']==21
        assert observations['language']==case['language'] and observations['mode']==case['mode']
        validate_rows(observations['rows'],case['language'])
        language=observations['language'];mode=observations['mode']
        if language not in decoders:
            model = read(case['model'])
            decoders[language] = {order:Decoder(inverse,emissions,model,order=order) for order in (2,0)}
        predictions=[]
        for row in observations['rows']:
            RUN_CONTEXT.update(row_id=row['id'])
            assert set(row)=={'id','block','observations'}
            methods = {label:decoders[language][order].decode(row['observations'],mode)
                       for label,order in (('trigram',2),('unigram',0))}
            predictions.append(dict(id=row['id'],block=row['block'],**methods))
        path = f'{OUT}/predictions/{case["name"]}.json'
        write(path,dict(case=case['name'],test_sha256=sha(case['test']),model_sha256=sha(case['model']),
                        rows=predictions,truth_parsed=False));outputs.append(path)
        print(json.dumps(dict(predicted=case['name'],rows=len(predictions))),flush=True)
    write(OUT+'/prediction-freeze.json',dict(status='ALL_PREDICTIONS_SAVED_BEFORE_SCORING',
          hashes=hashes(outputs),prepared_freeze_sha256=sha(DATA+'/prepared-freeze.json'),
          input_freeze_sha256=sha(DATA+'/input-freeze.json'),emission_rows_checked=emission_checks,
          truth_access='Only integrity hashing; no truth/source/trace parsed by prediction stage'))
    print(json.dumps(dict(status='PREDICTIONS_COMPLETE',files=len(outputs))),flush=True)


def score():
    prepared=verify_prepared()
    verify_prediction_freeze()
    rows,_=load_table(ROOT/TABLE);inverse=CodebookInverse(build_base_relation(rows))
    inverse_cache={};summary=[];all_details=[]
    for case in prepared['cases']:
        RUN_CONTEXT.update(case=case['name'])
        tests=read(case['test'])['rows'];truth=read(case['truth'])['rows']
        prediction=read(f'{OUT}/predictions/{case["name"]}.json')
        assert prediction['case']==case['name'] and prediction['test_sha256']==sha(case['test'])
        assert prediction['model_sha256']==sha(case['model']) and prediction['truth_parsed'] is False
        predicted=prediction['rows']
        for collection in (tests,truth,predicted):validate_rows(collection,case['language'])
        details=[]
        for obs,answer,pred in zip(tests,truth,predicted):
            RUN_CONTEXT.update(row_id=obs['id'])
            assert obs['id']==answer['id']==pred['id'] and obs['block']==answer['block']==pred['block']
            assert ''.join(answer['payloads'])==answer['plaintext'] and len(answer['plaintext'])==128
            assert len(answer['payloads'])==len(obs['observations'])
            inverse_info=[]
            for cipher,plain in zip(obs['observations'],answer['payloads']):
                key=(case['mode'],cipher)
                if key not in inverse_cache:inverse_cache[key]=inverse.analyze(cipher,gap_mode=case['mode'])
                info=inverse_cache[key]
                assert info['complete']
                contains=inverse._membership(info['dag'],plain)
                inverse_info.append(dict(contains_truth=contains,distinct=info['distinct_plaintext_count']))
            base=dict(id=obs['id'],block=obs['block'],support=all(x['contains_truth'] for x in inverse_info),
                      ambiguous_tokens=sum(x['distinct']>1 for x in inverse_info),tokens=len(inverse_info))
            for method in ('trigram','unigram'):
                result=pred[method]
                assert len(result['payloads'])==len(answer['payloads'])
                assert all(result['payloads']) and ''.join(result['payloads'])==result['plaintext']
                assert len(result['token_inverses'])==len(obs['observations'])
                assert result['all_local_payloads_unique']==all(i['distinct']==1 for i in inverse_info)
                correct=[a==b for a,b in zip(result['payloads'],answer['payloads'])]
                for cipher,payload,a,b in zip(obs['observations'],result['payloads'],result['token_inverses'],inverse_info):
                    info=inverse_cache[case['mode'],cipher]
                    assert a['observation']==cipher and a['path_count']==info['path_count']
                    assert a['distinct_plaintext_count']==b['distinct']
                    assert inverse._membership(info['dag'],payload)
                details.append(dict(**base,method=method,edit_distance=edit_distance(result['plaintext'],answer['plaintext']),
                                    true_characters=len(answer['plaintext']),exact=result['plaintext']==answer['plaintext'],
                                    correct_tokens=sum(correct),ambiguous_correct=sum(c for c,i in zip(correct,inverse_info) if i['distinct']>1),
                                    states=result['states'],transitions=result['transitions']))
        for method in ('trigram','unigram'):
            for block in ('A','B'):
                values=[d for d in details if d['method']==method and d['block']==block]
                assert len(values)==32
                edits=sum(d['edit_distance'] for d in values);characters=sum(d['true_characters'] for d in values)
                exact=sum(d['exact'] for d in values);supported=sum(d['support'] for d in values)
                ambiguous=[d for d in values if d['ambiguous_tokens']]
                atokens=sum(d['ambiguous_tokens'] for d in values)
                gate=supported==len(values) and 100*edits<=characters and 5*exact>=4*len(values)
                summary.append(dict(case=case['name'],language=case['language'],variant=case['variant'],mode=case['mode'],
                        method=method,block=block,rows=len(values),support_rows=supported,
                        edit_distance=edits,true_characters=characters,character_accuracy=1-edits/characters,
                        exact_rows=exact,exact_row_accuracy=exact/len(values),
                        tokens=sum(d['tokens'] for d in values),correct_tokens=sum(d['correct_tokens'] for d in values),
                        ambiguous_tokens=atokens,ambiguous_correct_tokens=sum(d['ambiguous_correct'] for d in values),
                        ambiguous_rows=len(ambiguous),ambiguous_exact_rows=sum(d['exact'] for d in ambiguous),
                        primary=method=='trigram' and case['variant'].startswith('ordered'),gate_pass=gate,
                        max_states=max(d['states'] for d in values),max_transitions=max(d['transitions'] for d in values)))
        all_details.append(dict(case=case['name'],rows=details))
    primary=[s for s in summary if s['primary']]
    assert len(summary)==64 and len(primary)==24 and len(all_details)==16
    final=dict(status='KNOWN_KEY_CALIBRATION_COMPLETE_NOT_A_TRANSLATION',
               all_primary_gates_pass=all(s['gate_pass'] for s in primary),primary_gates=len(primary),
               primary_gates_passed=sum(s['gate_pass'] for s in primary),summary=summary,
               observed_unique_mode_token_pairs=len(inverse_cache),
               original_word_boundaries_recovered=False,manuscript_or_unknown_key_run=False)
    write(OUT+'/row-scores.json',all_details);write(OUT+'/summary.json',final)
    write(OUT+'/score-freeze.json',dict(status='SCORED',hashes=hashes([OUT+'/row-scores.json',OUT+'/summary.json']),
              prediction_freeze_sha256=sha(OUT+'/prediction-freeze.json')))
    print(json.dumps({k:v for k,v in final.items() if k!='summary'}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('freeze','prepare','predict','score'))
    args=parser.parse_args()
    try:globals()[args.stage]()
    except Exception as exc:
        if args.stage!='freeze':
            failure=ROOT/OUT/f'{args.stage}-failure.json'
            if not failure.exists():write(str(failure.relative_to(ROOT)),dict(stage=args.stage,type=type(exc).__name__,message=str(exc),context=RUN_CONTEXT))
        raise
