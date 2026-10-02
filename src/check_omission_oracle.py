"""Independent finite enumeration checks for020 latent-word inference.

Only synthetic objects are used. No prepared/source/manuscript file is read.
Character automata are checked separately; this enumerates every concrete word
sequence to check the OOV aggregation, top-two paths, marginals and row events.
"""
from collections import Counter
import copy
import itertools
import json
import math

import numpy as np
import omission_oracle as production

COUNTS = Counter()
MAX_ERROR = 0.


def near(actual, expected):
    global MAX_ERROR
    error=abs(actual-expected)
    MAX_ERROR=max(MAX_ERROR,error)
    assert error<1e-9,(actual,expected,error)
    COUNTS['numerical_comparisons']+=1


class FiniteWords:
    def __init__(self, values, lexicon):
        assert abs(sum(values.values())-1)<1e-12
        self.values,self.lexicon=values,lexicon

    def log_word_probability(self, word):
        return math.log(self.values[word])

    def analyze(self, skeleton):
        matched=[w for w in self.values if production.transform(w)==skeleton]
        known=sorted(w for w in matched if w in self.lexicon)
        unknown=sorted((w for w in matched if w not in self.lexicon),
                       key=lambda w:(-self.values[w],w))
        mass=sum(self.values[w] for w in unknown)
        return dict(lexicon_candidates=known,log_oov_char_mass=math.log(mass) if mass else None,
                    top_oov=[dict(word=w,log_char_probability=self.log_word_probability(w)) for w in unknown[:2]],
                    numerical_complete=True,search_complete=True)


def model_fixture():
    variants=['b'+a+'t'+b for a in 'aeiou' for b in 'aeiou'][:20]
    values={w:.89/len(variants) for w in variants}
    values.update(bt=.01,ct=.08,a=.01,ae=.01)
    info=dict(words={'bt':1,'ct':99},lexical_N=100,
              context={'bt':1,'ct':5},pairs={'bt':{'ct':1},'ct':{'bt':5}},
              mixture_lexical=.9,mixture_character=.1,bigram_tau=100)
    return production.SourceModel(info,FiniteWords(values,set(info['words'])))


def independent_p0(model, word):
    return .9*model.words.get(word,0)/model.N+.1*model.character.values[word]


def independent_path_probability(model, words):
    result=independent_p0(model,words[0])
    for u,v in zip(words,words[1:]):
        result*=(model.pairs.get(u,{}).get(v,0)+100*independent_p0(model,v))/(model.context.get(u,0)+100)
    return result


def enumerate_paths(model, observed, rule):
    possibilities=[[w for w in model.character.values
                    if (w if rule=='identity' else production.transform(w))==s] for s in observed]
    paths=[(tuple(ws),independent_path_probability(model,ws)) for ws in itertools.product(*possibilities)]
    paths.sort(key=lambda x:(-x[1],x[0]))
    return paths


def check_segment(model, observed, rule='omit-aeiou'):
    actual=production.segment_inference(model,observed,rule)
    paths=enumerate_paths(model,observed,rule)
    evidence=sum(p for _,p in paths)
    near(actual['log_evidence'],math.log(evidence))
    for got,(words,p) in zip(actual['paths'],paths[:2]):
        near(got['log_prior'],math.log(p))
        near(got['posterior'],p/evidence)
        # Exact tied OOV spellings can differ at machine-roundoff in arithmetic;
        # their prior and posterior must still attain the enumerated ranks.
        returned=next(value for ws,value in paths if ws==tuple(got['words']))
        near(returned,p)
    for i,word in enumerate(actual['words']):
        for cand in word['candidates']:
            p=sum(value for ws,value in paths if ws[i]==cand['word'])/evidence
            near(cand['posterior'],p)
        oov=sum(value for ws,value in paths if ws[i] not in model.words)/evidence
        near(word['oov_group_posterior'],oov)
        assert word['posterior']<=word['oov_group_posterior']+1e-12 or word['prediction_in_lexicon']
    chosen=tuple(actual['paths'][0]['words'])
    for start in range(len(observed)):
        for stop in range(start+1,len(observed)+1):
            expected=sum(p for ws,p in paths if ws[start:stop]==chosen[start:stop])/evidence
            got,_=production.fixed_span_posterior(model,actual,start,stop)
            near(got,expected)
    COUNTS['complete_segments']+=1
    COUNTS['enumerated_concrete_paths']+=len(paths)
    return actual


def check_concrete_oov_winner():
    # The reverse of the aggregate-mass counterexample: a particular OOV word
    # must win against an available lexical word, even inside a latent sequence.
    values=dict(bt=.01,bata=.70,beta=.08,ct=.19,a=.01,ae=.01)
    info=dict(words={'bt':1,'ct':99},lexical_N=100,
              context={'bt':1,'ct':5},pairs={'bt':{'ct':1},'ct':{'bt':5}},
              mixture_lexical=.9,mixture_character=.1,bigram_tau=100)
    model=production.SourceModel(info,FiniteWords(values,set(info['words'])))
    for observed in [['bt'],['bt','ct'],['ct','bt'],['bt','bt','ct']]:
        result=check_segment(model,observed)
        expected=[('bata' if s=='bt' else s) for s in observed]
        assert result['paths'][0]['words']==expected
        for skeleton,word in zip(observed,result['words']):
            if skeleton!='bt':continue
            assert not word['prediction_in_lexicon']
            assert word['lexical_candidates']==1
            assert 0<word['posterior']<word['oov_group_posterior']<1
            assert any(c['word']=='bt' and c['lexicon'] for c in word['candidates'])
        COUNTS['concrete_oov_winner_segments']+=1


def check_exact_float_lexicographic_ties():
    # Equal probabilities and identical arithmetic give exactly equal represented
    # scores; unlike near-ties, these must obey the promised lexical path order.
    fixtures=[
        (dict(bat=.25,bet=.25,bit=.25,ct=.25),{'ct':1},False),
        (dict(bat=.125,bet=.125,ct=.75),{'bat':1,'bet':1},True),
    ]
    for values,words,known in fixtures:
        info=dict(words=words,lexical_N=sum(words.values()),context={},pairs={},
                  mixture_lexical=.9,mixture_character=.1,bigram_tau=100)
        model=production.SourceModel(info,FiniteWords(values,set(words)))
        lattice=model.lattice('bt','omit-aeiou')
        assert lattice['candidates']==['bat','bet']
        assert lattice['candidate_logs'][0]==lattice['candidate_logs'][1]
        for length in (1,2,3):
            observations=['bt']*length
            result=check_segment(model,observations)
            paths=enumerate_paths(model,observations,'omit-aeiou')
            assert result['paths'][0]['log_prior']==result['paths'][1]['log_prior']
            assert result['log_gap']==0
            assert [p['words'] for p in result['paths']]==[list(ws) for ws,_ in paths[:2]]
            assert result['paths'][0]['words']==['bat']*length
            assert result['paths'][1]['words']==['bat']*(length-1)+['bet']
            assert all(w['prediction_in_lexicon']==known for w in result['words'])
            COUNTS['exact_float_lexicographic_ties']+=1


def row(locus, words, resets=(), split='development', page='toy'):
    return dict(locus=locus,page=page,split=split,paragraph_first_line=False,
                tokens=[dict(units=list(w) if w is not None else None,reset_before=i in resets)
                        for i,w in enumerate(words)])


def check_layout(model):
    data=dict(rule='omit-aeiou',rows=[
        row('toy.1',['bt','ct'],(0,)),
        row('toy.2',['bt',None,'ct','bt'],(3,)),
        row('next.1',['bt','ct'],(),page='next'),
        row('r.1',['bt','ct','bt','ct'],(0,),split='reserve',page='r'),
    ])
    result=production.decode(model,data)
    expected=[[[0,0],[0,1],[1,0]],[[1,2]],[[1,3]],[[2,0],[2,1]],[[3,0],[3,1],[3,2],[3,3]]]
    assert [s['positions'] for s in result['segments']]==expected
    assert result['rows'][1]['tokens'][1] is None
    assert result['rows'][1]['exact_output_row_posterior'] is None
    altered=copy.deepcopy(data)
    altered['truth_rows']=[{'hidden_word':'a completely different answer'}]
    altered['mechanism']='different-truth-tag'
    assert production.decode(model,altered)==result
    COUNTS['same_observation_different_truth_checks']+=1
    for rno,r in enumerate(result['rows']):
        if not r['complete_known_row']:continue
        expected_probability=1.
        for s in result['segments']:
            covered=[i for i,(rr,_) in enumerate(s['positions']) if rr==rno]
            if not covered:continue
            observed=[''.join(data['rows'][rr]['tokens'][tt]['units']) for rr,tt in s['positions']]
            paths=enumerate_paths(model,observed,'omit-aeiou')
            desired={i:r['tokens'][s['positions'][i][1]]['prediction'] for i in covered}
            mass=sum(p for ws,p in paths if all(ws[i]==w for i,w in desired.items()))
            expected_probability*=mass/sum(p for _,p in paths)
        near(r['exact_output_row_posterior'],expected_probability)
    COUNTS['layout_rows']+=len(result['rows'])


def check_internal_reset_full_row(model):
    observations=dict(rule='omit-aeiou',rows=[
        row('toy.reset',['bt','ct','bt','ct'],(0,2)),
        row('r.reset',['bt','ct','bt','ct'],(0,2),split='reserve',page='r'),
    ])
    prediction=production.decode(model,observations)
    assert [s['positions'] for s in prediction['segments']]==[
        [[0,0],[0,1]],[[0,2],[0,3]],[[1,0],[1,1]],[[1,2],[1,3]]]
    paths=enumerate_paths(model,['bt','ct'],'omit-aeiou')
    assert paths[0][0]==('bt','ct')
    segment_probability=paths[0][1]/sum(p for _,p in paths)
    expected_row_probability=segment_probability*segment_probability
    assert 0<expected_row_probability<segment_probability<1
    for r in prediction['rows']:
        assert r['complete_known_row'] and len(r['tokens'])==4
        assert [t['prediction'] for t in r['tokens']]==['bt','ct','bt','ct']
        near(r['exact_output_row_posterior'],expected_row_probability)

    correct_truth=dict(rows=copy.deepcopy(observations['rows']))
    correct=production.score(model,prediction,correct_truth,observations)
    for split in production.SPLITS:
        metrics=correct['tests'][split]
        assert metrics['eligible_complete_rows']==metrics['exact_complete_rows']==1
        assert metrics['exact_complete_row_accuracy']==1
        assert metrics['excluded_rows']==0 and metrics['eligible_row_word_coverage']==1

    # Only the second fragment is wrong. A correct first fragment must not be
    # substituted for the whole physical row or counted as a second eligible row.
    changed_truth=copy.deepcopy(correct_truth)
    changed_truth['rows'][0]['tokens'][2]['units']=list('bata')
    changed=production.score(model,prediction,changed_truth,observations)
    development=changed['tests']['development']
    assert development['known_words']==4 and development['exact_words']==3
    assert development['eligible_complete_rows']==1 and development['exact_complete_rows']==0
    assert development['exact_complete_row_accuracy']==0 and development['excluded_rows']==0
    assert development['row_length_strata']=={'4-7':dict(rows=1,exact=0)}
    assert changed['tests']['reserve']['exact_complete_row_accuracy']==1
    # Neither truth variant may alter a subsequent decode, including model caches.
    assert production.decode(model,observations)==prediction
    COUNTS['same_observation_different_truth_checks']+=1
    COUNTS['internal_reset_complete_rows']+=len(prediction['rows'])
    COUNTS['evaluation_cases']+=2


def check_evaluation(model):
    truth=dict(rows=[row('toy.1',['bt','ct','ae','bata'],(0,)),
                     row('toy.2',['bt',None,'ct'],(0,)),
                     row('r.1',['ct','bt','a','beta'],(0,),split='reserve',page='r')])
    identity=dict(rule='identity',rows=copy.deepcopy(truth['rows']))
    prediction=production.decode(model,identity)
    scores=production.score(model,prediction,truth,identity)
    for split,t in scores['tests'].items():
        assert t['identity_exact'] and t['exact_word_accuracy']==1 and t['char_edit_accuracy']==1
        assert t['exact_complete_row_accuracy']==1 and t['eligible_complete_rows']==1
    omitted=dict(rule='omit-aeiou',rows=copy.deepcopy(truth['rows']))
    for r in omitted['rows']:
        for t in r['tokens']:
            if t['units'] is not None:t['units']=list(production.transform(''.join(t['units'])))
    prediction=production.decode(model,omitted)
    scores=production.score(model,prediction,truth,omitted)
    assert scores['tests']['development']['excluded_rows']==1
    assert scores['tests']['development']['known_words']==6
    assert scores['tests']['development']['nontrivial_words']==5
    assert scores['tests']['development']['baseline_exact_word_accuracy']==5/6
    assert scores['tests']['reserve']['nontrivial_words']==3
    for s in scores['true_path_diagnostics']:assert 0<=s['posterior']<=1
    COUNTS['evaluation_cases']+=2


def check_real_character_integration():
    """Connect both production models using only an analytic synthetic prior."""
    from omission_channel import CompletionModel

    alphabet='ab '
    probabilities=np.array([.25,.5,.25])
    tables=[np.tile(np.log(probabilities),len(alphabet)**(n-1)).reshape((len(alphabet),)*n)
            for n in range(1,5)]
    words={'a':4,'ab':3,'bb':1}
    info=dict(words=words,lexical_N=8,context={'a':2},pairs={'a':{'ab':2}},
              mixture_lexical=.9,mixture_character=.1,bigram_tau=100)
    character=CompletionModel(alphabet,tables,set(words))
    model=production.SourceModel(info,character)

    # Independent closed forms: conditioned nonempty memoryless words and the
    # geometric sum of arbitrary a insertions in each gap around retained b's.
    def skeleton(word):
        return word.replace('a','') or word

    def g(word):
        return (.25/.75)*.25**word.count('a')*.5**word.count('b')

    def gmass(observation,rule):
        if rule=='identity' or 'b' not in observation:return g(observation)
        return (.25/.75)*.5**len(observation)/.75**(len(observation)+1)

    def q0(word):
        return .9*words.get(word,0)/8+.1*g(word)

    def groups(observation,rule):
        if rule=='identity':return [(observation,q0(observation))]
        known=[w for w in words if skeleton(w)==observation]
        result=[(w,q0(w)) for w in known]
        oov=gmass(observation,rule)-sum(g(w) for w in known)
        if oov>0:result.append((None,.1*oov))
        return result

    def q1(left,right,right_mass):
        return (info['pairs'].get(left,{}).get(right,0)+100*right_mass)/(info['context'].get(left,0)+100)

    def enumerate_groups(observed,rule):
        output=[]
        for states in itertools.product(*(groups(s,rule) for s in observed)):
            value=states[0][1]
            for left,right in zip(states,states[1:]):value*=q1(left[0],right[0],right[1])
            output.append((states,value))
        return output

    for observed in ('a','aa','b','bb'):
        analysis=character.analyze(observed)
        assert analysis['search_complete'] and analysis['numerical_complete']
        known=sorted(w for w in words if skeleton(w)==observed)
        assert analysis['lexicon_candidates']==known
        near(analysis['char_mass'],gmass(observed,'omit-aeiou'))
        near(analysis['oov_char_mass'],gmass(observed,'omit-aeiou')-sum(g(w) for w in known))
        for w in known:
            near(math.exp(character.log_word_probability(w)),g(w))
        if observed=='a':assert analysis['top_oov']==[]
        if observed=='aa':assert [v['word'] for v in analysis['top_oov']]==['aa']
        if observed=='b':assert [v['word'] for v in analysis['top_oov']]==['b','ba']
        if observed=='bb':assert [v['word'] for v in analysis['top_oov']]==['abb','bab']

    truth=dict(rows=[row('int.1',['a','ab','aa','bb'],(0,2)),
                     row('int.2',['ab',None,'ba','a']),
                     row('int.r',['bb','ab','a','aa'],(0,2),split='reserve',page='r')])
    expected_positions=[[[0,0],[0,1]],[[0,2],[0,3],[1,0]],[[1,2],[1,3]],
                        [[2,0],[2,1]],[[2,2],[2,3]]]
    for rule in ('identity','omit-aeiou'):
        observations=dict(rule=rule,rows=copy.deepcopy(truth['rows']))
        for r in observations['rows']:
            for t in r['tokens']:
                if t['units'] is not None and rule=='omit-aeiou':
                    t['units']=list(skeleton(''.join(t['units'])))
        prediction=production.decode(model,observations)
        assert [s['positions'] for s in prediction['segments']]==expected_positions
        assert prediction['rows'][1]['tokens'][1] is None
        assert prediction['rows'][1]['exact_output_row_posterior'] is None
        for segment in prediction['segments']:
            positions=segment['positions']
            observed=[''.join(observations['rows'][r]['tokens'][t]['units']) for r,t in positions]
            enumerated=enumerate_groups(observed,rule)
            evidence=math.fsum(value for _,value in enumerated)
            near(segment['log_evidence'],math.log(evidence))
            for path in segment['paths']:
                prior=q0(path['words'][0])
                for left,right in zip(path['words'],path['words'][1:]):prior*=q1(left,right,q0(right))
                near(path['log_prior'],math.log(prior))
                near(path['posterior'],prior/evidence)
            for index,(r,t) in enumerate(positions):
                detail=prediction['rows'][r]['tokens'][t]
                for candidate in detail['candidates']:
                    target=candidate['word'] if rule=='identity' or candidate['word'] in words else None
                    groupmass=dict(groups(observed[index],rule))[target]
                    marginal=math.fsum(value for states,value in enumerated if states[index][0]==target)/evidence
                    near(candidate['posterior'],marginal*q0(candidate['word'])/groupmass)
                if observed[index] in ('a','aa'):near(detail['posterior'],1.)
            COUNTS['integration_analytic_group_paths']+=len(enumerated)
        scores=production.score(model,prediction,truth,observations)
        development=scores['tests']['development']
        reserve=scores['tests']['reserve']
        assert development['known_words']==7 and development['eligible_complete_rows']==1
        assert development['excluded_rows']==1 and development['exact_complete_row_accuracy']==1
        assert reserve['known_words']==4 and reserve['exact_complete_row_accuracy']==1
        if rule=='identity':
            assert development['identity_exact'] and reserve['identity_exact']
        else:
            assert development['exact_words']==6 and development['word_edit_distance']==2
            near(development['char_edit_accuracy'],10/12)
            assert development['omitted_characters']==3 and development['omitted_gap_edit_distance']==2
            assert development['nontrivial_words']==4 and reserve['nontrivial_words']==2
            assert prediction['rows'][1]['tokens'][2]['prediction']=='ab'
        assert production.decode(model,observations)==prediction
        COUNTS['real_character_integration_cases']+=1
        COUNTS['evaluation_cases']+=1


def main():
    model=model_fixture()
    for observations in [['bt'],['bt','ct'],['ct','bt','ct'],['bt','bt','ct'],['a','ae']]:
        check_segment(model,observations)
    check_segment(model,['bata','ct','ae'],'identity')
    a=production.segment_inference(model,['bt'],'omit-aeiou')
    assert a['words'][0]['oov_group_posterior']>.8
    assert a['paths'][0]['words']==['bt'], 'OOV aggregate mass was mistaken for a word'
    check_concrete_oov_winner()
    check_exact_float_lexicographic_ties()
    check_layout(model)
    check_internal_reset_full_row(model)
    check_evaluation(model)
    check_real_character_integration()
    print(json.dumps(dict(status='PASS',scope='SYNTHETIC_ONLY_NO_SOURCE_OR_MANUSCRIPT_READS',
                          checks=dict(COUNTS),maximum_absolute_error=MAX_ERROR),indent=2))


if __name__=='__main__':main()
