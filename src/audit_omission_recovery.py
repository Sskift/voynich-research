"""Independent020 inference/evaluation audit, including explicitly partial runs.

No production module is imported. Character word probabilities are recomputed
from frozen018 tables. For non-singleton omission observations, forward/backward
uses ARCHIVED OOV TOTAL MASS and the archived OOV top-two candidate set. Agreement
therefore does not independently certify infinite completion mass or search.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import argparse
import ast
import datetime
import hashlib
import itertools
import json
import math
import re
import sys
import traceback

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = 'data/omission-oracle'
RESULTS = 'results/omission-oracle'
LANGUAGES = ('la', 'de', 'it')
MECHANISMS = ('ordered', 'iid', 'copy')
RULES = ('identity', 'omit-aeiou')
SPLITS = ('development', 'reserve')
VOWELS = frozenset('aeiou')
COUNTS = Counter()
MAX_ERROR = defaultdict(float)


class AuditError(Exception):
    pass


def need(condition, message):
    if not condition:
        raise AuditError(message)


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def close(expected, actual, where, category='probability'):
    need(isinstance(actual, (int, float)) and math.isfinite(actual), where+': nonfinite value')
    error = abs(float(expected)-actual)
    MAX_ERROR[category] = max(MAX_ERROR[category], error)
    tolerance = (2e-8+2e-12*abs(expected)) if category=='log' else 3e-9
    need(error <= tolerance, f'{where}: expected {expected!r}, archived {actual!r}, error {error}')
    COUNTS['numeric_comparisons'] += 1


def compare(expected, actual, where):
    if isinstance(expected, dict):
        need(isinstance(actual, dict) and expected.keys()==actual.keys(), where+': fields differ')
        for key in expected:
            compare(expected[key], actual[key], where+'.'+key)
    elif isinstance(expected, list):
        need(isinstance(actual, list) and len(expected)==len(actual), where+': length differs')
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f'{where}[{i}]')
    elif isinstance(expected, float):
        close(expected, actual, where)
    else:
        need(expected==actual, f'{where}: expected {expected!r}, archived {actual!r}')


def logsum(values, axis=None):
    values = np.asarray(values, dtype=float)
    peak = np.max(values, axis=axis, keepdims=True)
    answer = peak + np.log(np.exp(values-peak).sum(axis=axis, keepdims=True))
    return float(answer.ravel()[0]) if axis is None else np.squeeze(answer, axis=axis)


def scalar_sum(*values):
    peak = max(values)
    return peak+math.log(math.fsum(math.exp(v-peak) for v in values))


def probability(logp):
    need(math.isfinite(logp) and logp<=2e-8, 'Invalid posterior log probability')
    return min(1., math.exp(logp))


def image(word, rule):
    need(bool(word), 'Empty reference/candidate word')
    return word if rule=='identity' else (''.join(c for c in word if c not in VOWELS) or word)


class Integrity:
    def __init__(self):
        self.files = {}
        self.manifests = []

    def verify(self, filename):
        document = read(filename)
        need(isinstance(document.get('hashes'), dict) and document['hashes'], filename+': no hashes')
        for path, expected in document['hashes'].items():
            need(digest(path)==expected, 'Frozen hash mismatch: '+path)
            if path in self.files:
                need(self.files[path]==expected, 'Conflicting freezes: '+path)
            self.files[path] = expected
            COUNTS['hash_checks'] += 1
        self.manifests.append(filename)
        return document


def audit_truth_access_order():
    """Static frozen-code constraint plus output hashes; not a runtime access log."""
    tree = ast.parse((ROOT/'src/omission_oracle.py').read_text())
    functions = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    run = functions['run']
    calls = [n for n in ast.walk(run) if isinstance(n, ast.Call)]
    writes = [n.lineno for n in calls if isinstance(n.func, ast.Name) and n.func.id=='write'
              and n.args and isinstance(n.args[0], ast.Name) and n.args[0].id=='predpath']
    hashes = [n.lineno for n in calls if isinstance(n.func, ast.Name) and n.func.id=='sha'
              and n.args and isinstance(n.args[0], ast.Name) and n.args[0].id=='predpath']
    truth_reads = [n.lineno for n in calls if isinstance(n.func, ast.Name) and n.func.id=='read'
                   and n.args and isinstance(n.args[0], ast.Name) and n.args[0].id=='truthpath']
    need(len(writes)==len(hashes)==len(truth_reads)==1 and writes[0]<hashes[0]<truth_reads[0],
         'Prediction write/hash does not precede truth read in frozen run()')
    for name in ('decode', 'segment_inference', 'fixed_span_posterior'):
        for node in ast.walk(functions[name]):
            need(not (isinstance(node, ast.Name) and 'truth' in node.id.lower()),
                 'Truth variable in inference: '+name)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                need(node.func.id not in ('read','open','load_model'), 'File read within inference: '+name)
    return dict(frozen_production_sha256=digest('src/omission_oracle.py'),
                prediction_write_line=writes[0], prediction_hash_line=hashes[0], truth_read_line=truth_reads[0],
                evidence='Frozen code order and saved prediction/score hashes; not an independent runtime file-access trace.')


class IndependentModel:
    def __init__(self, language, integrity):
        self.info = read(f'{DATA}/{language}/model.json')
        info = self.info
        need(info['language']==language and info['complete'] and info['NOT_A_TRANSLATION'], 'Invalid model metadata')
        self.words, self.pairs, self.context = info['words'], info['pairs'], info['context']
        self.total = info['lexical_N']
        need(sum(self.words.values())==self.total>0, 'Lexical count total')
        need(info['mixture_lexical']==.9 and info['mixture_character']==.1 and info['bigram_tau']==100,
             'Changed prior weights')
        need({u:sum(vs.values()) for u,vs in self.pairs.items()}==self.context, 'Pair/context count mismatch')
        need(info['pair_N']==sum(self.context.values()), 'Pair total mismatch')
        need(all(u in self.words and all(v in self.words and type(c) is int and c>0 for v,c in vs.items())
                 for u,vs in self.pairs.items()), 'Invalid lexical pair table')
        self.alphabet = info['alphabet']
        need(self.alphabet.endswith(' ') and len(set(self.alphabet))==len(self.alphabet), 'Alphabet')
        self.index = {letter:i for i,letter in enumerate(self.alphabet)}
        self.space = self.index[' ']
        self.size = len(self.alphabet)
        old = read(f'data/gap-inference/{language}/language.json')
        need(old['alphabet']==self.alphabet and old['word_counts']==self.words, '018 lexical/alphabet mismatch')
        path = f'data/gap-inference/{language}/conditional.bin'
        need(info['char_model']['conditional_file']==path and
             info['char_model']['conditional_sha256']==integrity.files[path], 'Character-model provenance')
        self.tables = []
        raw = np.fromfile(ROOT/path, dtype='<f8')
        need(len(raw)==sum(self.size**n for n in range(1,5)) and np.isfinite(raw).all(), '018 table dimensions')
        offset = 0
        for order in range(1,5):
            length = self.size**order
            table = raw[offset:offset+length].copy().reshape((-1,self.size))
            maximum = table.max(axis=1, keepdims=True)
            normalizer = maximum+np.log(np.exp(table-maximum).sum(axis=1,keepdims=True))
            table -= normalizer
            self.tables.append(table.reshape((self.size,)*order))
            offset += length
        first = self.tables[1][self.space]
        self.log_nonempty = math.log(math.fsum(math.exp(first[i]) for i in range(self.size) if i!=self.space))
        self.lexical_images = defaultdict(list)
        for word in sorted(self.words):
            self.lexical_images[image(word,'omit-aeiou')].append(word)

    @lru_cache(maxsize=None)
    def char_log(self, word):
        need(word and all(c in self.index and c!=' ' for c in word), 'Invalid concrete candidate')
        history = [self.space]
        terms = []
        for letter in word+' ':
            value = self.index[letter]
            terms.append(float(self.tables[len(history)][tuple(history+[value])]))
            history = (history+[value])[-3:]
        return math.fsum(terms)-self.log_nonempty

    @lru_cache(maxsize=None)
    def q0(self, word):
        char = math.log(.1)+self.char_log(word)
        count = self.words.get(word,0)
        return scalar_sum(math.log(.9)+math.log(count)-math.log(self.total),char) if count else char

    @lru_cache(maxsize=None)
    def q1(self, previous, word):
        if previous is None:
            return self.q0(word)
        count = self.context.get(previous,0)
        pair = self.pairs.get(previous,{}).get(word,0)
        numerator = math.log(100)+self.q0(word)
        if pair:
            numerator = scalar_sum(math.log(pair),numerator)
        return numerator-math.log(count+100)

    def path_score(self, words):
        return self.q0(words[0])+math.fsum(self.q1(u,v) for u,v in zip(words,words[1:]))

    def inspect_analysis(self, skeleton, rule, analysis, require_complete=True):
        known = ([skeleton] if skeleton in self.words else []) if rule=='identity' else self.lexical_images.get(skeleton,[])
        need(analysis['lexicon_candidates']==known, 'Incomplete/extra lexical candidates: '+skeleton)
        if require_complete:
            need(analysis['numerical_complete'] is True and analysis['search_complete'] is True,
                 'Incomplete archived completion: '+skeleton)
        singleton = rule=='identity' or all(c in VOWELS for c in skeleton)
        if singleton:
            expected_oov = None if known else self.char_log(skeleton)
            if expected_oov is None:
                need(analysis['log_oov_char_mass'] is None, 'Singleton known word has OOV mass')
            else:
                close(expected_oov, analysis['log_oov_char_mass'], 'Singleton OOV mass', 'log')
            need([x['word'] for x in analysis['top_oov']]==([] if known else [skeleton]), 'Singleton candidate set')
        else:
            need(not any(c in VOWELS for c in skeleton), 'Mixed omission observation')
            need(analysis['log_oov_char_mass'] is not None and len(analysis['top_oov'])==2,
                 'Infinite OOV support lost')
        seen = set()
        for item in analysis['top_oov']:
            word = item['word']
            need(word not in self.words and word not in seen and image(word,rule)==skeleton, 'Invalid OOV candidate')
            close(self.char_log(word), item['log_char_probability'], 'OOV character score', 'log')
            seen.add(word)
        if rule!='identity':
            need(analysis['skeleton']==skeleton, 'Certificate skeleton mismatch')
            for word,value in analysis['lexicon_log_char_probabilities'].items():
                close(self.char_log(word),value,'Lexical character score','log')
            need(set(analysis['lexicon_log_char_probabilities'])==set(known), 'Lexical probability inventory')
            pieces = [self.char_log(w) for w in known]
            if analysis['log_oov_char_mass'] is not None:pieces.append(analysis['log_oov_char_mass'])
            close(logsum(pieces), analysis['log_char_mass'], 'Archived partition consistency','log')
            for logkey,valuekey in [('log_char_mass','char_mass'),('log_oov_char_mass','oov_char_mass')]:
                expected = 0. if analysis[logkey] is None else math.exp(analysis[logkey])
                close(expected, analysis[valuekey], 'Ordinary/log mass consistency')
        COUNTS['completion_analyses'] += 1

    def lattice(self, skeleton, rule, certificates):
        key = rule+'|'+skeleton
        if rule=='identity':
            known = skeleton in self.words
            analysis = dict(lexicon_candidates=[skeleton] if known else [],
                            log_oov_char_mass=None if known else self.char_log(skeleton),
                            top_oov=[] if known else [dict(word=skeleton,log_char_probability=self.char_log(skeleton))],
                            numerical_complete=True,search_complete=True)
            if key in certificates:
                self.inspect_analysis(skeleton,rule,certificates[key])
        else:
            need(key in certificates, 'Missing completion certificate: '+key)
            analysis = certificates[key]
            self.inspect_analysis(skeleton,rule,analysis)
            COUNTS['lattices_using_archived_oov_mass'] += 1
        labels = list(analysis['lexicon_candidates'])
        weights = [self.q0(w) for w in labels]
        if analysis['log_oov_char_mass'] is not None:
            labels.append(None); weights.append(math.log(.1)+analysis['log_oov_char_mass'])
        candidates = sorted(analysis['lexicon_candidates']+[x['word'] for x in analysis['top_oov']])
        return dict(labels=labels,logweights=np.array(weights),candidates=candidates,
                    group={w:labels.index(w if w in self.words else None) for w in candidates},analysis=analysis)

    def transition_matrix(self, left, right):
        matrix = np.empty((len(left['labels']),len(right['labels'])))
        for i,u in enumerate(left['labels']):
            denominator = math.log(self.context.get(u,0)+100)
            for j,v in enumerate(right['labels']):
                numerator = math.log(100)+right['logweights'][j]
                pair = self.pairs.get(u,{}).get(v,0)
                if pair:numerator = scalar_sum(math.log(pair),numerator)
                matrix[i,j] = numerator-denominator
        return matrix


def layout_segments(rows):
    segments, current = [], []
    prior_header = None
    for rno,row in enumerate(rows):
        need(row['split'] in SPLITS and row['page'] not in ('f49v','f58r','f58v'), 'Unexpected scope')
        need(row['tokens'], 'Empty physical row')
        header = (row['page'],row['split'])
        if prior_header!=header and current:
            segments.append(current); current=[]
        for tno,token in enumerate(row['tokens']):
            need(type(token['reset_before']) is bool, 'Invalid reset flag')
            if token['reset_before'] or token['units'] is None:
                if current:segments.append(current); current=[]
            if token['units'] is not None:
                need(isinstance(token['units'],list) and token['units'], 'Invalid units')
                current.append((rno,tno,''.join(token['units'])))
        prior_header=header
    if current:segments.append(current)
    return segments


def forward_backward(lattices, matrices):
    # Scale at every position, independently of the production unscaled logs.
    normalizers = [logsum(lattices[0]['logweights'])]
    alpha = [lattices[0]['logweights']-normalizers[0]]
    for matrix in matrices:
        new = logsum(alpha[-1][:,None]+matrix,axis=0)
        scale = logsum(new)
        normalizers.append(scale); alpha.append(new-scale)
    beta = [None]*len(lattices)
    beta[-1] = np.zeros(len(lattices[-1]['labels']))
    for i in range(len(lattices)-2,-1,-1):
        beta[i] = logsum(matrices[i]+beta[i+1][None,:],axis=1)-normalizers[i+1]
    return alpha,beta,math.fsum(normalizers)


def constrained_evidence(model,lattices,matrices,clamps):
    # Generic evidence under concrete word constraints, instead of the production
    # alpha/internal-transition/beta formula for a special contiguous span.
    current=lattices[0]['logweights'].copy()
    total=0.
    for i,lat in enumerate(lattices):
        if i:current=logsum(current[:,None]+matrices[i-1],axis=0)
        if i in clamps:
            word=clamps[i]; group=lat['group'][word]
            selected=current[group]+model.q0(word)-lat['logweights'][group]
            current=np.full(len(current),-np.inf); current[group]=selected
        scale=logsum(current)
        total+=scale; current-=scale
    return total


@lru_cache(maxsize=None)
def distance(left,right):
    grid=[[0]*(len(right)+1) for _ in range(len(left)+1)]
    for i in range(len(left)+1):grid[i][0]=i
    for j in range(len(right)+1):grid[0][j]=j
    for i,a in enumerate(left,1):
        for j,b in enumerate(right,1):
            grid[i][j]=min(grid[i-1][j]+1,grid[i][j-1]+1,grid[i-1][j-1]+(a!=b))
    return grid[-1][-1]


def ratio(numerator,denominator):
    return numerator/denominator if denominator else None


def evaluation(prediction,truth,observations):
    output={}
    for split in SPLITS:
        records=[]; eligible=[]; excluded=0; strata={}
        for p,t,o in zip(prediction['rows'],truth['rows'],observations['rows']):
            if t['split']!=split:continue
            row_results=[]
            for pp,tt,oo in zip(p['tokens'],t['tokens'],o['tokens']):
                if tt['units'] is None:continue
                true=''.join(tt['units']); seen=''.join(oo['units']); guessed=pp['prediction']
                good=true==guessed; hard=bool(set(seen)-VOWELS)
                row_results.append(good)
                gaps=re.split('[^aeiou]',true); guessed_gaps=re.split('[^aeiou]',guessed)
                missing_count=missing_error=0
                if prediction['rule']=='omit-aeiou' and hard:
                    need(len(gaps)==len(guessed_gaps), 'Consonant gaps changed')
                    missing_count=sum(map(len,gaps))
                    missing_error=sum(distance(a,b) for a,b in zip(gaps,guessed_gaps))
                records.append(dict(good=good,length=len(true),edit=distance(true,guessed),
                                    baseline_edit=distance(true,seen),baseline_good=true==seen,
                                    changed=true!=seen,hard=hard,posterior=pp['posterior'],
                                    missing_count=missing_count,missing_error=missing_error))
            count=len(t['tokens'])
            if count>=4 and len(row_results)==count:
                exact=all(row_results); eligible.append((count,exact))
                label='4-7' if count<=7 else '8-15' if count<=15 else '16+'
                if label not in strata:strata[label]=dict(rows=0,exact=0)
                strata[label]['rows']+=1; strata[label]['exact']+=int(exact)
            else:excluded+=1
        changed=[r for r in records if r['changed']]
        hard=[r for r in records if r['hard']]
        confident=[r for r in hard if r['posterior']>=.9]
        def total(field,items=records):return sum(r[field] for r in items)
        def accuracy(items):
            denominator=total('length',items)
            return 1-total('edit',items)/denominator if denominator else None
        def bins(items):
            out=[dict(words=0,correct=0) for _ in range(4)]
            for item in items:
                index=sum(item['posterior']>=threshold for threshold in (.5,.9,.99))
                out[index]['words']+=1; out[index]['correct']+=int(item['good'])
            return out
        exact_rows=sum(x[1] for x in eligible)
        values=dict(known_words=len(records),exact_words=total('good'),reference_characters=total('length'),
                    word_edit_distance=total('edit'),char_edit_accuracy=accuracy(records),
                    exact_word_accuracy=ratio(total('good'),len(records)),lossy_words=len(changed),
                    exact_lossy_words=total('good',changed),exact_lossy_word_accuracy=ratio(total('good',changed),len(changed)),
                    nontrivial_words=len(hard),lossy_char_edit_accuracy=accuracy(changed),
                    nontrivial_exact_word_accuracy=ratio(total('good',hard),len(hard)),
                    nontrivial_char_edit_accuracy=accuracy(hard),nontrivial_confident_words=len(confident),
                    nontrivial_confident_correct=total('good',confident),
                    nontrivial_confident_accuracy=ratio(total('good',confident),len(confident)),
                    nontrivial_confident_coverage=ratio(len(confident),len(hard)),confidence_bins=bins(records),
                    nontrivial_confidence_bins=bins(hard),eligible_complete_rows=len(eligible),
                    exact_complete_rows=exact_rows,excluded_rows=excluded,
                    eligible_row_word_coverage=ratio(sum(x[0] for x in eligible),len(records)),row_length_strata=strata,
                    exact_complete_row_accuracy=ratio(exact_rows,len(eligible)),omitted_characters=total('missing_count'),
                    omitted_gap_edit_distance=total('missing_error'),
                    omitted_symbol_edit_accuracy=1-total('missing_error')/total('missing_count') if total('missing_count') else None,
                    baseline_char_edit_accuracy=1-total('baseline_edit')/total('length') if total('length') else None,
                    baseline_exact_word_accuracy=ratio(total('baseline_good'),len(records)))
        denominators=bool(total('length') and records and eligible and hard and confident)
        values['ordered_recovery_thresholds_pass']=bool(denominators and accuracy(records)>=.95
            and values['exact_word_accuracy']>=.9 and values['exact_complete_row_accuracy']>=.5
            and values['nontrivial_confident_accuracy']>=.95 and values['nontrivial_confident_coverage']>=.25)
        values['high_confidence_error_flag']=total('good',confident)/len(confident)<.95 if confident else None
        values['identity_exact']=bool(records and eligible and total('good')==len(records)
                                      and exact_rows==len(eligible) and total('edit')==0)
        output[split]=values
    return output


def audit_case(name,model,certificates,input_hash):
    language,tail=name.split('/')
    mechanism=next(m for m in MECHANISMS if tail.startswith(m+'-'))
    rule=tail[len(mechanism)+1:]
    paths={kind:f'{DATA}/{kind}/{name}.json' for kind in ('tests','truth')}
    predpath=f'{RESULTS}/predictions/{name}.json'; scorepath=f'{RESULTS}/scores/{name}.json'
    observations=read(paths['tests']); prediction=read(predpath); score=read(scorepath)
    need(prediction['complete'] and prediction['status']=='ORACLE_KEY_PREDICTIONS_NOT_A_TRANSLATION', 'Prediction incomplete')
    need(prediction['rule']==rule and observations['rule']==rule, 'Rule mismatch')
    need(prediction['observation_sha256']==digest(paths['tests']) and prediction['model_sha256']==digest(f'{DATA}/{language}/model.json')
         and prediction['input_freeze_sha256']==input_hash, 'Prediction provenance')
    need(score['prediction_sha256']==digest(predpath) and score['truth_sha256']==digest(paths['truth']), 'Score provenance')
    need(score['case']==name and score['truth_used_only_after_prediction'] is True
         and score['status']=='ORACLE_RECOVERY_NOT_A_MANUSCRIPT_READING', 'Score status')
    need(len(prediction['rows'])==len(observations['rows']), 'Prediction row count')
    for p,o in zip(prediction['rows'],observations['rows']):
        need(all(p[k]==o[k] for k in ('locus','page','split')), 'Prediction row identity')
        need(len(p['tokens'])==len(o['tokens']), 'Prediction token count')
        for pp,oo in zip(p['tokens'],o['tokens']):
            need((pp is None)==(oo['units'] is None), 'Known/None mismatch')
            if pp is not None:
                need(image(pp['prediction'],rule)==''.join(oo['units']), 'Prediction outside true channel')
    segments=layout_segments(observations['rows'])
    need(len(segments)==len(prediction['segments']), 'Segment count')
    lattice_cache={}; row_logs=defaultdict(float); recomputed_segments=[]
    for sn,(positions,archived) in enumerate(zip(segments,prediction['segments'])):
        need(archived['index']==sn and archived['positions']==[[r,t] for r,t,_ in positions], 'Segment layout')
        lattices=[]
        for _,_,skeleton in positions:
            if skeleton not in lattice_cache:lattice_cache[skeleton]=model.lattice(skeleton,rule,certificates)
            lattices.append(lattice_cache[skeleton])
        matrices=[model.transition_matrix(left,right) for left,right in zip(lattices,lattices[1:])]
        alpha,beta,evidence=forward_backward(lattices,matrices)
        close(evidence,archived['log_evidence'], 'Segment evidence','log')
        paths_found=archived['paths']
        need(1<=len(paths_found)<=2, 'Missing/extra top paths')
        expected_number=1 if all(len(lat['candidates'])==1 for lat in lattices) else 2
        need(len(paths_found)==expected_number, 'Top-two path count')
        need(len({tuple(p['words']) for p in paths_found})==len(paths_found), 'Duplicate top paths')
        path_scores=[]
        for path in paths_found:
            need(len(path['words'])==len(positions), 'Top path length')
            need(all(w in lat['candidates'] for w,lat in zip(path['words'],lattices)), 'Top path candidate missing')
            value=model.path_score(path['words']); path_scores.append(value)
            close(value,path['log_prior'],'Concrete top path prior','log')
            close(probability(value-evidence),path['posterior'],'Concrete top path posterior')
            COUNTS['top_path_scores']+=1
        if len(paths_found)==2:
            need(path_scores[0]>=path_scores[1]-2e-8, 'Top paths not ordered')
            saved_gap=paths_found[0]['log_prior']-paths_found[1]['log_prior']
            close(path_scores[0]-path_scores[1],archived['log_gap'],'Top path gap','log')
            need(archived['top2_numerical_tie']==(saved_gap<=1e-12), 'Numerical tie flag')
        else:
            need(archived['log_gap'] is None and archived['top2_numerical_tie'] is False, 'Singleton gap')
        chosen=paths_found[0]['words']; row_clamps=defaultdict(dict)
        for i,((r,t,skeleton),lat) in enumerate(zip(positions,lattices)):
            item=prediction['rows'][r]['tokens'][t]
            need(item['prediction']==chosen[i] and item['segment']==sn and item['segment_offset']==i, 'Chosen path/layout mismatch')
            need(item['observation']==skeleton and item['completion_certificate_key']==rule+'|'+skeleton, 'Certificate reference')
            marginals=np.exp(alpha[i]+beta[i]); close(1.,float(marginals.sum()),'Marginal normalization')
            candidates=item['candidates']
            need(sorted(x['word'] for x in candidates)==lat['candidates'], 'Candidate inventory')
            need(candidates==sorted(candidates,key=lambda x:(-x['posterior'],x['word'])), 'Candidate order')
            shown_oov=0.; chosen_posterior=None
            for candidate in candidates:
                word=candidate['word']; group=lat['group'][word]
                post=probability(float(alpha[i][group]+beta[i][group]+model.q0(word)-lat['logweights'][group]))
                close(post,candidate['posterior'],'Concrete candidate posterior')
                need(candidate['lexicon']==(word in model.words), 'Candidate lexical flag')
                if word not in model.words:shown_oov+=post
                if word==chosen[i]:chosen_posterior=post
                COUNTS['candidate_posteriors']+=1
            need(chosen_posterior is not None, 'Chosen word omitted')
            close(chosen_posterior,item['posterior'],'Chosen posterior')
            need(item['prediction_in_lexicon']==(chosen[i] in model.words), 'Chosen lexical flag')
            need(item['lexical_candidates']==len(lat['analysis']['lexicon_candidates']), 'Lexical ambiguity count')
            group_post=float(marginals[lat['labels'].index(None)]) if None in lat['labels'] else 0.
            close(group_post,item['oov_group_posterior'],'OOV group posterior')
            close(max(0.,group_post-shown_oov),item['unlisted_oov_posterior'],'Unlisted OOV posterior')
            row_clamps[r][i]=chosen[i]
            COUNTS['word_posteriors']+=1
        for r,clamps in row_clamps.items():
            row_logs[r]+=constrained_evidence(model,lattices,matrices,clamps)-evidence
        recomputed_segments.append(dict(index=sn,evidence=evidence,map_score=path_scores[0],chosen=chosen))
        COUNTS['segments']+=1
    for rno,row in enumerate(prediction['rows']):
        complete=all(t is not None for t in row['tokens'])
        need(row['complete_known_row']==complete, 'Known row flag')
        if complete:close(probability(row_logs[rno]),row['exact_output_row_posterior'],'Whole physical row posterior')
        else:need(row['exact_output_row_posterior'] is None,'Unknown row posterior claimed')
        COUNTS['physical_rows']+=1
    # Truth is read only after the independent inference checks above.
    truth=read(paths['truth'])
    need(len(truth['rows'])==len(observations['rows']), 'Truth row count')
    for t,o in zip(truth['rows'],observations['rows']):
        need({k:v for k,v in t.items() if k!='tokens'}=={k:v for k,v in o.items() if k!='tokens'}, 'Truth header/layout')
        need(len(t['tokens'])==len(o['tokens']), 'Truth token count')
        for tt,oo in zip(t['tokens'],o['tokens']):
            need(tt['reset_before']==oo['reset_before'] and (tt['units'] is None)==(oo['units'] is None), 'Truth boundary/None')
            if tt['units'] is not None:need(image(''.join(tt['units']),rule)==''.join(oo['units']), 'Truth channel mismatch')
    tests=evaluation(prediction,truth,observations)
    compare(tests,score['tests'],name+'.metrics')
    need(len(score['true_path_diagnostics'])==len(segments), 'True-path diagnostic count')
    for positions,recomputed,archived in zip(segments,recomputed_segments,score['true_path_diagnostics']):
        words=[''.join(truth['rows'][r]['tokens'][t]['units']) for r,t,_ in positions]
        value=model.path_score(words)
        need(value<=recomputed['map_score']+2e-8, 'True path outranks archived MAP')
        need(archived['segment']==recomputed['index'] and archived['exact_map']==(words==recomputed['chosen']), 'True path label')
        close(value,archived['log_prior'],'True path prior','log')
        close(probability(value-recomputed['evidence']),archived['posterior'],'True path posterior')
        COUNTS['true_path_scores']+=1
    COUNTS['completed_cases']+=1
    return dict(status='PASS',rule=rule,tests=tests,unique_observations=len(lattice_cache),
                prediction_sha256=digest(predpath),score_sha256=digest(scorepath),
                certificate_keys=[rule+'|'+s for s in lattice_cache])


def audit_failure(document,models,input_hash):
    need(document['status']=='COMPUTATION_INCOMPLETE_NO_GATE_OR_NATIVE', 'Unexpected failure status')
    need(document['input_freeze_sha256']==input_hash, 'Failure input hash')
    context=document['run_context']; name=context['case']; language=context['language']
    need(name==language+'/ordered-omit-aeiou' and context['phase']=='decode-observations'
         and context['rule']=='omit-aeiou', 'Unexpected partial-run failure location')
    need(not (ROOT/f'{RESULTS}/predictions/{name}.json').exists()
         and not (ROOT/f'{RESULTS}/scores/{name}.json').exists(), 'Failed case has a completed prediction/score')
    skeleton=document['observation']; analysis=document['analysis']
    need(analysis['status']=='NUMERICAL_INCOMPLETE' and analysis['search_complete'] is True
         and analysis['numerical_complete'] is False, 'Failure completeness flags')
    certificate=analysis['numerical_certificate']
    need(certificate['numerical_complete'] is False and certificate['numerical_log_tolerance']==1e-8
         and certificate['conservative_log_mass_error_bound']>certificate['numerical_log_tolerance'], 'Failure threshold decision')
    need(certificate['vowel_escape_lower_bound']>0 and certificate['max_scaled_linear_residual']>=0, 'Failure certificate validity')
    search=analysis['search_counters']
    need(search['expansion_cap']==100000 and search['heap_cap']==1000000
         and search['expanded_prefixes']<100000 and search['heap_peak']<=1000000, 'Search budget changed')
    models[language].inspect_analysis(skeleton,'omit-aeiou',analysis,require_complete=False)
    observations=read(f'{DATA}/tests/{name}.json')
    need(any(t['units'] is not None and ''.join(t['units'])==skeleton for r in observations['rows'] for t in r['tokens']),
         'Failure observation absent from case')
    return dict(case=name,observation=skeleton,numerical_complete=False,search_complete=True,
                archived_conservative_log_mass_error_bound=certificate['conservative_log_mass_error_bound'],
                fixed_tolerance=certificate['numerical_log_tolerance'],
                independently_checked='Frozen stopping decision, lexical/OOV concrete word probabilities and partition consistency only.',
                not_independently_certified='Infinite completion mass, the residual error-bound derivation and global OOV search completeness.',
                failure_sha256=digest(f'{RESULTS}/computation-failure.json'))


def audit():
    has_failure=(ROOT/f'{RESULTS}/computation-failure.json').exists()
    has_summary=(ROOT/f'{RESULTS}/summary.json').exists()
    need(has_failure!=has_summary,'Need exactly one terminal failure or complete summary')
    integrity=Integrity()
    inputs=integrity.verify(f'{DATA}/input-freeze.json')
    prepared=integrity.verify(f'{DATA}/prepared-manifest.json')
    input_hash=digest(f'{DATA}/input-freeze.json')
    need(prepared['input_hashes']==inputs['hashes'] and prepared['input_freeze_sha256']==input_hash, 'Prepared hash chain')
    need(prepared['complete'] is True and prepared['NOT_A_TRANSLATION'] is True
         and prepared['status']=='PREPARED_NOT_SCORED', 'Prepared completion status')
    all_cases=[f'{l}/{m}-{r}' for l,m,r in itertools.product(LANGUAGES,MECHANISMS,RULES)]
    expected={f'{DATA}/config.json'} | {f'{DATA}/{l}/model.json' for l in LANGUAGES}
    expected|={f'{DATA}/{k}/{name}.json' for k in ('tests','truth','cipher') for name in all_cases}
    need(len(expected)==58 and set(prepared['hashes'])==expected, 'Prepared output inventory')
    need(set(inputs['mandatory_inputs'])<=set(inputs['hashes']), 'Mandatory input hashes missing')
    config=read(f'{DATA}/config.json')
    need(config['cases']==all_cases and config['lexical_weight']==.9 and config['character_weight']==.1
         and config['bigram_tau']==100 and config['search_budget']==dict(expanded_prefixes=100000,heap_entries=1000000), 'Config')
    code_order=audit_truth_access_order()
    complete_names=[]
    for name in all_cases:
        pred=(ROOT/f'{RESULTS}/predictions/{name}.json').exists()
        score=(ROOT/f'{RESULTS}/scores/{name}.json').exists()
        need(pred==score, 'Unpaired prediction/score: '+name)
        if pred:complete_names.append(name)
    if has_failure:
        need(complete_names==['la/ordered-identity'], 'Unexpected partial-run inventory; review scope before auditing')
    else:
        need(complete_names==all_cases, 'Complete summary has missing cases')
        integrity.verify(f'{RESULTS}/result-freeze.json')
    # Optional terminal-failure freeze, written by the parent after the stopped run.
    for path in sorted((ROOT/RESULTS).glob('*freeze*.json')):
        relative=str(path.relative_to(ROOT))
        if relative not in integrity.manifests:integrity.verify(relative)
    models={}; results={}; observed_keys=defaultdict(set)
    for name in complete_names:
        language=name.split('/')[0]
        if language not in models:models[language]=IndependentModel(language,integrity)
        certpath=f'{RESULTS}/{language}-completion-certificates.json'
        certificates=read(certpath) if (ROOT/certpath).exists() else {}
        results[name]=audit_case(name,models[language],certificates,input_hash)
        observed_keys[language].update(results[name].pop('certificate_keys'))
    summary=None; failure=None
    if has_failure:
        failure=audit_failure(read(f'{RESULTS}/computation-failure.json'),models,input_hash)
        need(not (ROOT/f'{RESULTS}/result-freeze.json').exists(), 'Failed run has successful result freeze')
    else:
        for language,keys in observed_keys.items():
            need(set(read(f'{RESULTS}/{language}-completion-certificates.json'))==keys, 'Completion certificate inventory')
        stored=read(f'{RESULTS}/summary.json')
        ordered={l:{s:results[f'{l}/ordered-omit-aeiou']['tests'][s]['ordered_recovery_thresholds_pass'] for s in SPLITS}
                 for l in LANGUAGES}
        identity=all(results[f'{l}/{m}-identity']['tests'][s]['identity_exact'] for l,m,s in itertools.product(LANGUAGES,MECHANISMS,SPLITS))
        summary=dict(ordered_by_language_split=ordered,identity_all_exact=identity,
                     ordered_oracle_gate_pass=bool(identity and all(v for lang in ordered.values() for v in lang.values())),
                     confidence_flags={name:{s:r['tests'][s]['high_confidence_error_flag'] for s in SPLITS} for name,r in results.items()})
        need(stored['status']=='COMPLETE_ORACLE_CALIBRATION_NO_NATIVE_SCORES', 'Summary completion marker')
        for field,value in summary.items():compare(value,stored[field],'summary.'+field)
    return dict(status='PARTIAL_RUN_AUDIT_PASS_NOT_CALIBRATION' if has_failure else 'COMPLETE_RUN_INFERENCE_EVALUATION_AUDIT_PASS',
                complete_run=not has_failure,completed_cases=len(complete_names),planned_cases=len(all_cases),
                uncompleted_cases=[n for n in all_cases if n not in complete_names],
                manuscript_scoring_performed=False,calibration_claim=False if has_failure else summary['ordered_oracle_gate_pass'],
                dependency_scope=dict(character_word_probabilities='Independently recomputed from018 normalized P1..P4.',
                    non_singleton_oov_mass='Archived values only; infinite mass must be independently audited separately.',
                    oov_candidate_search='Archived top2 candidate set; concrete word scores independently recomputed, global optimality not re-searched.',
                    actual_completed_case_uses_archived_oov_mass=bool(COUNTS['lattices_using_archived_oov_mass'])),
                code_truth_separation=code_order,hash_manifests=integrity.manifests,distinct_frozen_files=len(integrity.files),
                checks=dict(COUNTS),maximum_absolute_errors=dict(MAX_ERROR),cases=results,summary=summary,failure=failure,
                input_freeze_sha256=input_hash,prepared_manifest_sha256=digest(f'{DATA}/prepared-manifest.json'),
                audit_source_sha256=digest('src/audit_omission_recovery.py'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='results/audit-016-omission-recovery.json')
    arguments=parser.parse_args()
    output=ROOT/arguments.output
    need(not output.exists(),'Audit output exists; preserve it and choose a new --output')
    try:
        report=audit()
    except Exception as error:
        report=dict(status='AUDIT_FAILED',error_type=type(error).__name__,error=str(error),
                    checks=dict(COUNTS),maximum_absolute_errors=dict(MAX_ERROR),traceback=traceback.format_exc(),
                    audit_source_sha256=digest('src/audit_omission_recovery.py'))
        stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        failed=output.with_name(output.stem+'.failure-'+stamp+'.json')
        failed.parent.mkdir(parents=True,exist_ok=True)
        with failed.open('x',encoding='utf-8') as handle:json.dump(report,handle,ensure_ascii=False,indent=2,allow_nan=False)
        print(json.dumps(report,ensure_ascii=False,indent=2),file=sys.stderr)
        raise
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8') as handle:
        json.dump(report,handle,ensure_ascii=False,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(status=report['status'],completed_cases=report['completed_cases'],
                          planned_cases=report['planned_cases'],checks=report['checks'],
                          maximum_absolute_errors=report['maximum_absolute_errors'],output=str(output)),indent=2))


if __name__=='__main__':
    main()
