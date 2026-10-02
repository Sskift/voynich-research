"""Finite, unpruned Viterbi under a declared approximate base52 score.

Input interface is observations, codebook inverse, emission scores, and a
separately trained LM. No seed, true text, original units, lengths, or paired
preserved-boundary observations are accepted or loaded here.
"""
from collections import Counter, defaultdict
import itertools
import math

from naibbe_channel import ALPHABET, SINGLE, DROP

EOS = '$'
BOS = '^'
TAU = 20
SEGMENT_LENGTH = 128
MAX_VITERBI_STATES = 1000000


def train_model(stream):
    assert stream and not set(stream)-set(ALPHABET)
    counts = [Counter(), defaultdict(Counter), defaultdict(Counter)]
    chunks = [stream[i:i+SEGMENT_LENGTH] for i in range(0,len(stream),SEGMENT_LENGTH)]
    for chunk in chunks:
        history = BOS+BOS
        for letter in chunk+EOS:
            counts[0][letter] += 1
            counts[1][history[-1]][letter] += 1
            counts[2][history][letter] += 1
            history = (history+letter)[-2:]
    alphabet = ALPHABET+EOS
    p0 = {c:(counts[0][c]+.5)/(sum(counts[0].values())+.5*len(alphabet)) for c in alphabet}
    p1 = {}
    for context in BOS+ALPHABET:
        row = counts[1][context]; total = sum(row.values())
        p1[context] = {c:(row[c]+TAU*p0[c])/(total+TAU) for c in alphabet}
    p2 = {}
    contexts = [BOS+BOS]+[BOS+c for c in ALPHABET]+[a+b for a in ALPHABET for b in ALPHABET]
    for context in contexts:
        row = counts[2][context]; total = sum(row.values())
        p2[context] = {c:(row[c]+TAU*p1[context[-1]][c])/(total+TAU) for c in alphabet}
    for row in [p0]+list(p1.values())+list(p2.values()):
        assert abs(sum(row.values())-1) < 1e-12 and min(row.values()) > 0
    return dict(alphabet=ALPHABET,order=2,tau=TAU,p0=p0,p1=p1,p2=p2,
                training_letters=len(stream),training_chunks=len(chunks),
                training_chunk_length=SEGMENT_LENGTH,original_word_spaces_encoded=False)


class Decoder:
    def __init__(self, inverse, emissions, model, order=2, max_states=MAX_VITERBI_STATES):
        assert order in (0,2)
        self.inverse, self.emissions, self.model = inverse, emissions, model
        self.order, self.max_states = order, max_states
        self.cache = {}
        self.char_log = {h:{c:math.log(p) for c,p in row.items()} for h,row in model['p2'].items()}
        self.unigram_log = {c:math.log(p) for c,p in model['p0'].items()}
        self.unit_cache = {}

    def _letters(self, history, letters):
        key = (history,letters)
        if key not in self.unit_cache:
            score, end = 0., history
            for c in letters:
                score += self.char_log[end][c] if self.order == 2 else self.unigram_log[c]
                end = (end+c)[-2:]
            self.unit_cache[key] = score, end
        return self.unit_cache[key]

    def analyze(self, observation, mode):
        key = (mode, observation)
        if key not in self.cache:
            self.cache[key] = self.inverse.analyze(observation, gap_mode=mode, truth=None)
        answer = self.cache[key]
        if not answer['complete']:
            raise RuntimeError('INVERSE_BUDGET_EXCEEDED '+repr((observation,answer['budget'])))
        if not answer['path_count']:
            raise RuntimeError('NO_LEGAL_CODEBOOK_PATH '+repr(observation))
        return answer

    def decode(self, observations, mode):
        assert observations and mode in ('preserved','deleted')
        total = sum(map(len,observations))
        outgoing = [[] for _ in range(total+1)]
        summaries, offset = [], 0
        for token_index, observation in enumerate(observations):
            answer = self.analyze(observation,mode)
            summaries.append(dict(observation=observation,path_count=answer['path_count'],
                                  distinct_plaintext_count=answer['distinct_plaintext_count']))
            for start, edges in enumerate(answer['dag']['edges']):
                for edge in edges:
                    end, plain = edge['to'], edge['output']
                    score = self.emissions[observation[start:end],plain]
                    # Only the very last one-letter chunk is forced singleton.
                    score += 0. if len(plain) == 1 and offset+end == total else math.log(SINGLE if len(plain) == 1 else 1-SINGLE)
                    if mode == 'deleted' and offset+end != total:
                        score += math.log(1-DROP if end == len(observation) else DROP)
                    outgoing[offset+start].append((offset+end,plain,score,token_index))
            offset += len(observation)
        states = [{} for _ in range(total+1)]
        states[0][BOS+BOS] = (0.,'',None)
        state_count, transitions = 1, 0
        for start in range(total):
            for history, (value,text,path) in states[start].items():
                for end,plain,channel,token_index in outgoing[start]:
                    language,new_history = self._letters(history,plain)
                    candidate = (value+channel+language,text+plain,(path,start,end,plain,token_index))
                    old = states[end].get(new_history)
                    if old is None:
                        state_count += 1
                        if state_count > self.max_states:
                            raise RuntimeError('VITERBI_BUDGET_EXCEEDED '+repr(dict(
                                cipher_position=end,states_used=state_count,limit=self.max_states)))
                    if old is None or candidate[0] > old[0] or (candidate[0] == old[0] and candidate[1] < old[1]):
                        states[end][new_history] = candidate
                    transitions += 1
        endings = []
        for history,(value,text,path) in states[-1].items():
            eos = self.char_log[history][EOS] if self.order == 2 else self.unigram_log[EOS]
            endings.append((value+eos,text,path))
        assert endings, 'No terminal Viterbi state'
        value,text,path = min(endings,key=lambda x:(-x[0],x[1]))
        units = []
        while path is not None:
            parent,start,end,plain,token_index = path
            units.append(dict(start=start,end=end,plaintext=plain,observed_token=token_index))
            path = parent
        units.reverse()
        payloads = ['']*len(observations)
        for unit in units: payloads[unit['observed_token']] += unit['plaintext']
        assert ''.join(payloads) == text and all(payloads)
        return dict(plaintext=text,payloads=payloads,units=units,score=value,
                    score_scope='UNPRUNED_FLOAT_VITERBI_APPROXIMATE_WITH_REPLACEMENT_CHANNEL_NOT_POSTERIOR',
                    states=state_count,transitions=transitions,token_inverses=summaries,
                    all_local_payloads_unique=all(x['distinct_plaintext_count']==1 for x in summaries))


def edit_distance(a,b):
    row = list(range(len(b)+1))
    for i,x in enumerate(a,1):
        new = [i]
        for j,y in enumerate(b,1):
            new.append(min(row[j]+1,new[-1]+1,row[j-1]+(x!=y)))
        row = new
    return row[-1]


def self_check():
    from naibbe_inverse import CodebookInverse, _reference_paths
    relation = [('x','a'),('y','b'),('xy','ab'),('xy','ba'),('xx','b')]
    inverse = CodebookInverse(relation)
    emissions = {pair:math.log(.2) for pair in relation}
    model = train_model('abbaababaababba'*8)
    checks = 0
    for order in (0,2):
        decoder = Decoder(inverse,emissions,model,order)
        for mode in ('preserved','deleted'):
            for observations in (['x'],['xy'],['x','y'],['xy','xx'],['xx','xy','y']):
                result = decoder.decode(observations,mode)
                # Independent brute force over semantic paths and a direct
                # complete-string score, including final forced singleton/EOS.
                def unit_paths(obs,offset=0):
                    if offset == len(obs): return [[]]
                    return [[(c,p)]+tail for c,p in relation if obs.startswith(c,offset)
                            and (mode=='deleted' or offset==0 and len(c)==len(obs))
                            for tail in unit_paths(obs,offset+len(c))]
                best = -math.inf
                for pieces in itertools.product(*(unit_paths(o) for o in observations)):
                    units = [item for piece in pieces for item in piece]
                    text = ''.join(p for _,p in units)
                    score,history = 0.,BOS+BOS
                    for letter in text+EOS:
                        score += math.log(model['p2'][history][letter] if order==2 else model['p0'][letter])
                        history = (history+letter)[-2:]
                    index = 0
                    for piece in pieces:
                        for j,(cipher,plain) in enumerate(piece):
                            score += emissions[cipher,plain]
                            final = index == len(units)-1
                            if not (len(plain)==1 and final): score += math.log(SINGLE if len(plain)==1 else 1-SINGLE)
                            if mode=='deleted' and not final: score += math.log(1-DROP if j==len(piece)-1 else DROP)
                            index += 1
                    best = max(best,score)
                assert abs(best-result['score']) < 1e-10,(best,result)
                checks += 1
    try: Decoder(inverse,emissions,model,max_states=1).decode(['x'],'preserved')
    except RuntimeError as exc: assert 'BUDGET' in str(exc)
    else: raise AssertionError('Budget failure missing')
    assert edit_distance('abc','ac') == 1 and edit_distance('ab','ba') == 2
    ties = Decoder(CodebookInverse([('x','a'),('x','b')]),
                   {('x','a'):math.log(.5),('x','b'):math.log(.5)},train_model('ab'),order=0)
    assert ties.decode(['x'],'preserved')['plaintext']=='a'
    return dict(status='PASS',brute_force_viterbi_cases=checks,budget_failure_checked=True,
                equal_score_local_tie_checked=True,real_codebook_or_source_read=False)


if __name__ == '__main__':
    import json
    print(json.dumps(self_check()))
