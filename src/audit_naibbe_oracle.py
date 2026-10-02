"""Independent021 input, finite-inverse, path-score and evaluation audit.

The reference calculations import no production/author module. An explicitly
separate isolation check replays the frozen production Decoder in a subprocess
after denying ALL file opens, using only preloaded public tables/models/tests.
That replay is not substituted for independent mathematical verification.
"""
import argparse
from collections import Counter, defaultdict, deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import traceback
import unicodedata
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
D = 'data/naibbe-oracle'
R = 'results/naibbe-oracle'
A = 'data/naibbe-review/author-f2675ec5'
LETTERS = 'abcdefghilmnopqrstuvxyz'
LETTERS = ''.join(c for c in 'abcdefghijklmnopqrstuvwxyz' if c not in 'jkw')
WEIGHTS = [('alpha',20),('beta1',8),('beta2',8),('beta3',8),('gamma1',4),('gamma2',4)]
MODES = ('preserved','deleted')
VARIANTS = ('ordered-r0','ordered-r1','ordered-r2','shuffled-r0')
STRING_CAP = 200000
COUNTS = Counter()
ERRORS = defaultdict(float)


def require(condition, reason):
    if not condition:
        raise AssertionError(reason)


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def compare(expected, actual, where='root'):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and expected.keys() == actual.keys(), where+': object keys')
        for key in expected:
            compare(expected[key], actual[key], where+'.'+str(key))
    elif isinstance(expected, (list, tuple)):
        require(isinstance(actual, (list, tuple)) and len(expected) == len(actual), where+': list length')
        for i, (a,b) in enumerate(zip(expected, actual)):
            compare(a,b,where+f'[{i}]')
    elif isinstance(expected, float):
        require(isinstance(actual, (float,int)) and math.isfinite(actual), where+': nonfinite')
        error = abs(expected-actual)
        ERRORS['general_numeric'] = max(ERRORS['general_numeric'], error)
        require(error <= 3e-12, where+f': {expected!r} != {actual!r}')
        COUNTS['numeric_comparisons'] += 1
    else:
        require(expected == actual, where+f': {expected!r} != {actual!r}')


class Integrity:
    def __init__(self):
        self.files = {}

    def check(self, path):
        data = read(path)
        require(isinstance(data.get('hashes'), dict) and data['hashes'], path+': empty manifest')
        for name, digest in data['hashes'].items():
            require(sha(name) == digest, 'Frozen SHA mismatch: '+name)
            require(name not in self.files or self.files[name] == digest, 'Conflicting SHA: '+name)
            self.files[name] = digest
            COUNTS['hash_assertions'] += 1
        return data


def cases():
    answer = []
    for language, variant, mode in itertools.product(('la','it'), VARIANTS, MODES):
        name = f'{language}-{variant}-{mode}'
        answer.append(dict(name=name,language=language,variant=variant,mode=mode,
                           test=f'{D}/tests/{name}.json',truth=f'{D}/truth/{name}.json',
                           model=f'{D}/{language}/model.json'))
    return answer


def configuration():
    return dict(experiment=21,channel='BASE52_AUTHOR_f2675ec5',source_chunk_length=128,
                test_chunks=64,blocks={'A':[0,32],'B':[32,64]},languages=['la','it'],
                modes=list(MODES),ordered_replicates=3,shuffle_replicates=1,
                ordered_seed='202102100+1000*language_index+replicate',
                shuffle_encoder_seed='202102190+1000*language_index',
                shuffle_source_seed='202102199+1000*language_index',known_key=True,
                matching_language_prior=True,score='APPROXIMATE_PATH_MAX_NOT_POSTERIOR',
                support_scope='EXACT_LOCAL_CODEBOOK_LANGUAGE_OVERAPPROXIMATES_STOCHASTIC_ENCODER_SUPPORT',
                model_order=2,baseline_order=0,tau=20,additive=.5,
                single_probability='17/36 except final forced singleton',gap_drop=.03,
                max_dfa_states=100000,max_dfa_edges=1000000,max_viterbi_states=1000000,
                max_encoder_draws_per_chunk=1000000,
                primary_gate=dict(coverage=1.,character_edit_accuracy=.99,exact_chunk_accuracy=.8),
                manuscripts_read=False,original_word_boundaries_recovered=False)


def input_paths():
    return {'notes/experiment-021-design.md',*[f'src/naibbe_{name}.py' for name in
            ('channel','decode','inverse','sources','oracle')],
            'data/raw/latin/caesar-gallic.xml','data/raw/latin/cicero-catilinam.xml',
            'data/language-screen/dante.json','data/language-screen/boccaccio.json',
            'data/raw/language-screen/dante-1012.txt','data/raw/language-screen/boccaccio-21424.html',
            'src/language_sources.py','data/latin-sources.json',A+'/naibbe_tables.csv',
            A+'/naibbe.py',A+'/README.md',A+'/LICENSE',A+'/archive-inventory.json',
            'notes/naibbe-channel-review.md','notes/naibbe-and-units-primary-review.md'}


def normalized(text):
    replacements = {'æ':'ae','Æ':'ae','œ':'oe','Œ':'oe','ð':'d','Ð':'d',
                    'þ':'th','Þ':'th','ł':'l','Ł':'l','ß':'ss','ø':'o','Ø':'o'}
    expanded = ''.join(replacements.get(c,c) for c in unicodedata.normalize('NFD',text)
                       if unicodedata.category(c) != 'Mn')
    upper = ''.join(c for c in expanded if c.isalpha()).upper()
    output = upper.translate(str.maketrans({'W':'UU','J':'I','K':'C'})).lower()
    require(set(output) <= set(LETTERS), 'Unsupported normalized source')
    return output


def tei_source(raw, subtype, numbers):
    """Independent tree walk: retain selected paragraphs and alternative paths."""
    tree = ET.fromstring(raw)
    local = lambda e: e.tag.rsplit('}',1)[-1] if isinstance(e.tag,str) else ''
    ignored = {'note','head','foreign','rdg','del','num'}
    def selected_children(node):
        children = list(node)
        if local(node) == 'app':
            selected = [c for c in children if local(c) == 'lem']
            require(selected, 'Apparatus without direct lemma')
            return selected[:1]
        if local(node) == 'choice':
            for tag in ('expan','reg'):
                selected = [c for c in children if local(c) == tag]
                if selected:
                    return selected[:1]
            require(children, 'Empty choice')
            return children[:1]
        return children
    def text(node):
        tag = local(node)
        if tag in ignored or not tag:
            return ''
        if tag in ('choice','app'):
            return text(selected_children(node)[0])
        terms = [node.text or '']
        for child in node:
            terms.extend((text(child),child.tail or ''))
        return ''.join(terms)
    paths = {}
    stack = [(tree,'/'+local(tree)+'[1]')]
    while stack:
        node,path = stack.pop(); paths[node] = path
        seen = Counter()
        for child in node:
            tag = local(child); seen[tag] += 1
            stack.append((child,path+f'/{tag}[{seen[tag]}]'))
    bodies = [e for e in tree.iter() if local(e) == 'body']
    require(len(bodies) == 1, 'Expected single body')
    sections = [(int(e.attrib['n']),e) for outer in bodies[0] if local(outer)=='div'
                for e in outer if local(e)=='div' and e.get('subtype')==subtype and int(e.attrib['n']) in numbers]
    require([n for n,_ in sections] == list(numbers), 'Unexpected source divisions')
    paragraphs, locators = [], []
    for number, division in sections:
        pending = [division]; selected = []
        while pending:
            node = pending.pop(); tag = local(node)
            if tag in ignored or not tag:
                continue
            if tag == 'p':
                selected.append(node)
            else:
                pending.extend(reversed(selected_children(node)))
        require(selected, 'Selected division has no paragraphs')
        for node in selected:
            paragraphs.append(text(node))
            locators.append(dict(locator=paths[node],division_subtype=subtype,division_n=number))
    return '\n'.join(paragraphs),locators


def rebuild_sources():
    result = {}
    specs = [('la','train','data/raw/latin/caesar-gallic.xml','book',tuple(range(1,8))),
             ('la','test','data/raw/latin/cicero-catilinam.xml','speech',(1,2,3,4)),
             ('it','train','data/language-screen/dante.json',None,None),
             ('it','test','data/language-screen/boccaccio.json',None,None)]
    for language, role, path, subtype, numbers in specs:
        raw = (ROOT/path).read_bytes()
        if language == 'la':
            text, locations = tei_source(raw,subtype,numbers)
        else:
            data = json.loads(raw)
            require(sha(data['raw_file']) == data['raw_sha256'], 'Italian raw provenance')
            COUNTS['raw_provenance_hashes'] += 1
            text = '\n'.join(r['text'] for r in data['records'])
            locations = []
            for i,row in enumerate(data['records']):
                location = dict(locator=f'/records/{i}/text')
                if 'locator' in row:
                    location['source_locator'] = row['locator']
                locations.append(location)
        require(text.strip(), 'Empty source')
        result[language,role] = dict(text=normalized(text),locations=locations,path=path,
                                     raw_sha256=hashlib.sha256(raw).hexdigest())
    return result


def language_model(stream):
    unigram, bigram, trigram = Counter(), Counter(), Counter()
    chunks = [stream[i:i+128] for i in range(0,len(stream),128)]
    for chunk in chunks:
        padded = '^^'+chunk+'$'
        unigram.update(chunk+'$')
        bigram.update(padded[i-1:i+1] for i in range(2,len(padded)))
        trigram.update(padded[i-2:i+1] for i in range(2,len(padded)))
    letters = LETTERS+'$'
    p0 = {c:Fraction(2*unigram[c]+1,2*sum(unigram.values())+len(letters)) for c in letters}
    p1 = {}
    for history in '^'+LETTERS:
        total = sum(bigram[history+c] for c in letters)
        p1[history] = {c:(bigram[history+c]+20*p0[c])/(total+20) for c in letters}
    p2 = {}
    for history in ['^^']+['^'+c for c in LETTERS]+[a+b for a in LETTERS for b in LETTERS]:
        total = sum(trigram[history+c] for c in letters)
        p2[history] = {c:float((trigram[history+c]+20*p1[history[-1]][c])/(total+20)) for c in letters}
    COUNTS['lm_training_characters'] += len(stream)
    COUNTS['lm_conditional_rows'] += 1+len(p1)+len(p2)
    return dict(alphabet=LETTERS,order=2,tau=20,p0={c:float(p) for c,p in p0.items()},
                p1={h:{c:float(p) for c,p in row.items()} for h,row in p1.items()},p2=p2,
                training_letters=len(stream),training_chunks=len(chunks),training_chunk_length=128,
                original_word_spaces_encoded=False)


def table_relation():
    import csv
    rows = list(csv.DictReader((ROOT/(A+'/naibbe_tables.csv')).read_text(encoding='utf-8-sig').splitlines()))
    table = {r['code']:r['glyphs'] for r in rows}
    expected = {f'{role}_{name}_{letter}' for role in ('unigram','prefix','suffix')
                for name,_ in WEIGHTS for letter in LETTERS}
    require(len(rows) == len(table) == 414 and set(table) == expected, 'Table inventory')
    singles = {value for code,value in table.items() if code.startswith('unigram_')}
    probabilities = {}
    for a in LETTERS:
        marginal = defaultdict(Fraction)
        for name,w in WEIGHTS:
            marginal[table[f'unigram_{name}_{a}']] += Fraction(w,52)
        require(sum(marginal.values()) == 1, 'Unigram emission row')
        for cipher,p in marginal.items():
            probabilities[cipher,a] = p
    for a,b in itertools.product(LETTERS,repeat=2):
        marginal = defaultdict(Fraction)
        for (left,w),(right,v) in itertools.product(WEIGHTS,repeat=2):
            cipher = table[f'prefix_{left}_{a}']+table[f'suffix_{right}_{b}']
            if cipher not in singles:
                marginal[cipher] += w*v
        total = sum(marginal.values())
        require(total > 0, 'No allowed emission for source bigram')
        for cipher,p in marginal.items():
            probabilities[cipher,a+b] = p/total
    COUNTS['emission_rows'] = len(LETTERS)+len(LETTERS)**2
    return rows,table,{key:math.log(float(value)) for key,value in probabilities.items()}


def replay(plaintext, table, generator):
    remaining, chunks = deque(plaintext), []
    while remaining:
        take = 1 if len(remaining)==1 else (1 if generator.random() < 17/36 else 2)
        chunks.append(''.join(remaining.popleft() for _ in range(take)))
    def shuffled():
        cards = sum(([name]*number for name,number in WEIGHTS),[])
        generator.shuffle(cards)
        return deque(cards)
    cards = shuffled()
    draws = []
    def card():
        nonlocal cards
        if not cards:
            cards = shuffled()
        require(len(draws) < 1000000, 'Reference encoding draw cap')
        name = cards.popleft(); draws.append(name)
        return name
    singles = {g for c,g in table.items() if c.startswith('unigram_')}
    units, retries = [], 0
    for chunk in chunks:
        if len(chunk)==1:
            codes = [f'unigram_{card()}_{chunk}']
        else:
            while True:
                codes = [f'prefix_{card()}_{chunk[0]}',f'suffix_{card()}_{chunk[1]}']
                if table[codes[0]]+table[codes[1]] not in singles:
                    break
                retries += 1
        units.append(dict(plaintext=chunk,cipher=''.join(table[c] for c in codes),codes=codes))
    deleted = [generator.random() < .03 for _ in units[1:]]
    grouped = [[]]
    for i, unit in enumerate(units):
        if i and not deleted[i-1]:
            grouped.append([])
        grouped[-1].append(unit)
    return dict(units=units,deleted=deleted,observed=[''.join(u['cipher'] for u in group) for group in grouped],
                payloads=[''.join(u['plaintext'] for u in group) for group in grouped],
                draws=draws,ambiguity_retries=retries)


def inspect_prepared(report):
    integrity = Integrity()
    frozen = integrity.check(D+'/input-freeze.json')
    require(set(frozen['hashes']) == input_paths(), 'Input manifest scope')
    require(frozen['status'] == 'FROZEN_BEFORE_SOURCE_PREPARATION_OR_SCORING', 'Input status')
    compare(configuration(),frozen['config'],'frozen config')
    prepared = integrity.check(D+'/prepared-freeze.json')
    require(prepared['status']=='PREPARED_NOT_PREDICTED' and
            prepared['input_freeze_sha256']==sha(D+'/input-freeze.json'), 'Preparation binding')
    compare(cases(),prepared['cases'],'cases')
    inventory = {D+'/config.json',D+'/source-provenance.json'}
    for case in cases():
        inventory.update((case['test'],case['truth'],case['model'],
                          f'{D}/generation/{case["language"]}-{case["variant"]}.json'))
    require(len(inventory)==44 and set(prepared['hashes'])==inventory, 'Prepared44 inventory')
    compare(configuration(),read(D+'/config.json'),'prepared config')
    sources, provenance = rebuild_sources(),read(D+'/source-provenance.json')
    rows,table,emissions = table_relation()
    models = {}
    for language_index, language in enumerate(('la','it')):
        train,test = sources[language,'train']['text'],sources[language,'test']['text']
        info = provenance[language]
        for role, text in (('train',train),('test',test)):
            source = sources[language,role]
            require(info['normalized_'+role+'_sha256']==hashlib.sha256(text.encode()).hexdigest(), 'Normalized source hash')
            require(info['normalized_'+role+'_letters']==len(text), 'Normalized source length')
            saved = info['provenance'][role]
            require(saved['path']==source['path'] and saved['sha256']==source['raw_sha256'], 'Source identity')
            compare(source['locations'],saved['locators'],'source locators')
            require(isinstance(saved['selection_rule'],str) and saved['selection_rule'], 'No source-selection rule')
        require(info['test_selected_interval']==[0,8192] and info['orthography']=='AUTHOR_NORMALIZATION_NO_UV_MERGE',
                'Source interval/orthography')
        models[language] = language_model(train)
        compare(models[language],read(f'{D}/{language}/model.json'),'LM '+language)
        require(len(test)>=8192, 'Short test source')
        chunks = [test[i:i+128] for i in range(0,8192,128)]
        mixed = []
        shuffler = random.Random(202102199+1000*language_index)
        for chunk in chunks:
            letters = list(chunk); shuffler.shuffle(letters); mixed.append(''.join(letters))
        for variant in VARIANTS:
            seed = (202102190 if variant=='shuffled-r0' else 202102100+int(variant[-1]))+1000*language_index
            generator = random.Random(seed)
            traces, tests, truths = [], {m:[] for m in MODES}, {m:[] for m in MODES}
            for i, plain in enumerate(mixed if variant=='shuffled-r0' else chunks):
                generated = replay(plain,table,generator)
                row_id, block = f'{language}-{i:03d}', 'A' if i<32 else 'B'
                traces.append(dict(id=row_id,plaintext=plain,**generated))
                for mode in MODES:
                    obs = [u['cipher'] for u in generated['units']] if mode=='preserved' else generated['observed']
                    pay = [u['plaintext'] for u in generated['units']] if mode=='preserved' else generated['payloads']
                    tests[mode].append(dict(id=row_id,block=block,observations=obs))
                    truths[mode].append(dict(id=row_id,block=block,plaintext=plain,payloads=pay))
                    COUNTS['observed_tokens'] += len(obs)
                COUNTS['generation_rows'] += 1
                COUNTS['plaintext_letters_encoded'] += len(plain)
                COUNTS['card_draws'] += len(generated['draws'])
                COUNTS['unigram_collision_retries'] += generated['ambiguity_retries']
                COUNTS['deleted_gaps'] += sum(generated['deleted'])
            compare(dict(language=language,variant=variant,seed=seed,rows=traces),
                    read(f'{D}/generation/{language}-{variant}.json'),'generation '+language+'/'+variant)
            for mode in MODES:
                name = f'{language}-{variant}-{mode}'
                compare(dict(experiment=21,language=language,mode=mode,rows=tests[mode]),read(f'{D}/tests/{name}.json'),'test '+name)
                compare(dict(rows=truths[mode]),read(f'{D}/truth/{name}.json'),'truth '+name)
            print(json.dumps(dict(stage='inputs',language=language,variant=variant,status='PASS')),file=sys.stderr,flush=True)
    report['prepared'] = dict(status='PASS',input_hashes=len(frozen['hashes']),prepared_hashes=len(prepared['hashes']),
                              direct_scope_only=True,reference_lm_arithmetic='Exact Fraction backoff then float conversion',
                              italian_scope='Frozen records.text and raw provenance hash; no reinterpretation of older HTML extraction')
    return integrity,rows,table,emissions,models


class InverseReference:
    """Full finite string-set enumeration, independent of production's DFA."""
    def __init__(self, emissions):
        mapping = defaultdict(set)
        for cipher,plain in emissions:
            mapping[cipher].add(plain)
        self.mapping = {c:sorted(ps) for c,ps in mapping.items()}
        self.longest = max(map(len,mapping))
        self.cache = {}

    def analyze(self, observation, mode):
        key = mode,observation
        if key in self.cache:
            return self.cache[key]
        n = len(observation)
        edges = [[] for _ in range(n+1)]
        for start in range(n):
            if mode=='preserved' and start:
                continue
            stops = [n] if mode=='preserved' else range(start+1,min(n,start+self.longest)+1)
            for end in stops:
                cipher = observation[start:end]
                edges[start].extend((end,plain,cipher) for plain in self.mapping.get(cipher,()))
        strings, paths = [set() for _ in range(n+1)],[0]*(n+1)
        strings[n],paths[n] = {''},1
        for start in range(n-1,-1,-1):
            for end,plain,_ in edges[start]:
                paths[start] += paths[end]
                strings[start].update(plain+tail for tail in strings[end])
                require(len(strings[start]) <= STRING_CAP, 'AUDIT_STRING_ENUMERATION_INCOMPLETE; no truncation allowed')
        require(paths[0] and strings[0], 'No codebook inverse for saved observation')
        result = dict(edges=edges,path_count=paths[0],strings=frozenset(strings[0]),distinct=len(strings[0]))
        self.cache[key] = result
        COUNTS['unique_inverse_analyses'] += 1
        COUNTS['inverse_candidate_strings'] += len(strings[0])
        COUNTS['maximum_local_candidate_count'] = max(COUNTS['maximum_local_candidate_count'],len(strings[0]))
        return result


def distance(a,b):
    # Full matrix is intentionally independent of production's rolling row.
    matrix = [[0]*(len(b)+1) for _ in range(len(a)+1)]
    for i in range(len(a)+1): matrix[i][0] = i
    for j in range(len(b)+1): matrix[0][j] = j
    for i in range(1,len(a)+1):
        for j in range(1,len(b)+1):
            matrix[i][j] = min(matrix[i-1][j]+1,matrix[i][j-1]+1,matrix[i-1][j-1]+(a[i-1]!=b[j-1]))
    return matrix[-1][-1]


def independent_path_check(observations, mode, inverse, emissions, model, order, archived):
    total = sum(map(len,observations))
    edges = [[] for _ in range(total+1)]
    expected_tokens,offset = [],0
    for index,observation in enumerate(observations):
        info = inverse.analyze(observation,mode)
        expected_tokens.append(dict(observation=observation,path_count=info['path_count'],distinct_plaintext_count=info['distinct']))
        for start,outgoing in enumerate(info['edges']):
            for end,plain,cipher in outgoing:
                terms = [emissions[cipher,plain]]
                final = offset+end==total
                if not (len(plain)==1 and final):
                    terms.append(math.log(17/36 if len(plain)==1 else 19/36))
                if mode=='deleted' and not final:
                    terms.append(math.log(.97 if end==len(observation) else .03))
                edges[offset+start].append((offset+end,plain,math.fsum(terms),index))
        offset += len(observation)
    def language(history,plain):
        terms = []
        for c in plain:
            terms.append(math.log(model['p2'][history][c] if order==2 else model['p0'][c]))
            history = (history+c)[-2:]
        return math.fsum(terms),history
    scores = [dict() for _ in range(total+1)]
    scores[0]['^^'] = 0.
    transitions = 0
    for start,outgoing in enumerate(edges):
        for history,value in scores[start].items():
            for end,plain,channel,_ in outgoing:
                lm,new = language(history,plain)
                candidate = math.fsum((value,channel,lm))
                if new not in scores[end] or candidate > scores[end][new]:
                    scores[end][new] = candidate
                transitions += 1
    optimum = max(value+language(h,'$')[0] for h,value in scores[-1].items())
    chosen_terms,history,position,payloads = [],'^^',0,['']*len(observations)
    for unit in archived['units']:
        require(unit['start']==position, 'Saved unit positions do not form a complete path')
        matching = [edge for edge in edges[position] if
                    edge[0]==unit['end'] and edge[1]==unit['plaintext'] and edge[3]==unit['observed_token']]
        require(len(matching)==1, 'Saved unit is absent/ambiguous in semantic DAG')
        end,plain,channel,index = matching[0]
        lm,history = language(history,plain)
        chosen_terms.extend((channel,lm));payloads[index] += plain;position = end
    require(position==total, 'Saved path stopped before end')
    chosen_terms.append(language(history,'$')[0])
    chosen = math.fsum(chosen_terms)
    error = max(abs(chosen-optimum),abs(chosen-archived['score']))
    ERRORS['path_score'] = max(ERRORS['path_score'],error)
    require(error <= 1e-8, f'Saved path is not maximum score: gap={error}')
    compare(payloads,archived['payloads'],'predicted payload ownership')
    require(''.join(payloads)==archived['plaintext'], 'Predicted plaintext/unit mismatch')
    compare(expected_tokens,archived['token_inverses'],'inverse counts')
    require(archived['states']==sum(map(len,scores)) and archived['transitions']==transitions, 'DP graph/state counters')
    require(archived['score_scope']=='UNPRUNED_FLOAT_VITERBI_APPROXIMATE_WITH_REPLACEMENT_CHANNEL_NOT_POSTERIOR','Score scope')
    require(archived['all_local_payloads_unique']==all(t['distinct_plaintext_count']==1 for t in expected_tokens),'Unique flag')
    COUNTS['prediction_method_rows'] += 1


def isolation_replay(rows,models,tests,predictions):
    """Additional production-only replay; public in-memory inputs, ALL opens denied."""
    program = r'''
import hashlib,json,sys
from naibbe_channel import emission_model
from naibbe_inverse import CodebookInverse,build_base_relation
from naibbe_decode import Decoder
data=json.load(sys.stdin)
attempts=[]
def guard(event,args):
    if event=='open':
        attempts.append(repr(args))
        raise RuntimeError('FORBIDDEN_FILE_OPEN_DURING_DECODER_REPLAY '+repr(args))
sys.addaudithook(guard)
table={r['code']:r['glyphs'] for r in data['rows']}
emissions,_=emission_model(table)
inverse=CodebookInverse(build_base_relation(data['rows']))
decoders={l:{o:Decoder(inverse,emissions,m,order=o) for o in (2,0)} for l,m in data['models'].items()}
answers=[]
for case in data['tests']:
    for row in case['rows']:
        for method,order in (('trigram',2),('unigram',0)):
            answer=decoders[case['language']][order].decode(row['observations'],case['mode'])
            digest=hashlib.sha256(json.dumps(answer,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
            answers.append([case['name'],row['id'],method,digest])
print(json.dumps({'digests':answers,'file_open_attempts':attempts}))
'''
    # Use saved, independently verified float models here to reproduce exact
    # production ties. The mathematical audit uses its Fraction-derived models.
    # One condition per fresh subprocess: even the paired observations and
    # earlier predictions are absent from that subprocess's memory.
    all_digests = []
    for case in tests:
        language = case['language']
        inputs = dict(rows=rows,models={language:read(f'{D}/{language}/model.json')},tests=[case])
        result = subprocess.run([sys.executable,'-B','-c',program],input=json.dumps(inputs),
                                text=True,capture_output=True,cwd=ROOT/'src')
        require(result.returncode==0,'Isolation subprocess failed: '+result.stderr[-4000:])
        output = json.loads(result.stdout)
        require(not output['file_open_attempts'],'Decoder attempted a forbidden open')
        all_digests.extend(output['digests'])
        print(json.dumps(dict(stage='isolation',case=case['name'],status='PASS')),file=sys.stderr,flush=True)
    expected = []
    for case in tests:
        for row in predictions[case['name']]['rows']:
            for method in ('trigram','unigram'):
                digest = hashlib.sha256(json.dumps(row[method],sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
                expected.append([case['name'],row['id'],method,digest])
    compare(expected,all_digests,'isolated exact production replay')
    return dict(status='PASS',rows=len(expected),file_opens_attempted=0,
                scope='Additional frozen production replay only; all reference mathematics separately implemented',
                denied='All open audit events after module import/input transfer; each fresh subprocess gets only one condition, with no source, truth, traces, paired observations or saved predictions')


def inspect_recovery(report,prepared_result):
    integrity,rows,table,emissions,models = prepared_result
    prediction_freeze = integrity.check(R+'/prediction-freeze.json')
    require(prediction_freeze['status']=='ALL_PREDICTIONS_SAVED_BEFORE_SCORING' and
            prediction_freeze['input_freeze_sha256']==sha(D+'/input-freeze.json') and
            prediction_freeze['prepared_freeze_sha256']==sha(D+'/prepared-freeze.json'),'Prediction freeze chain')
    require(set(prediction_freeze['hashes'])=={f'{R}/predictions/{c["name"]}.json' for c in cases()},'Prediction inventory')
    require(prediction_freeze['emission_rows_checked']==552,'Emission row count')
    scored = integrity.check(R+'/score-freeze.json')
    require(scored['status']=='SCORED' and scored['prediction_freeze_sha256']==sha(R+'/prediction-freeze.json'),'Score freeze chain')
    require(set(scored['hashes'])=={R+'/row-scores.json',R+'/summary.json'},'Score inventory')
    inverse = InverseReference(emissions)
    details,summary,tests_for_guard,predictions = [],[],[],{}
    for case in cases():
        observation = read(case['test']); truth = read(case['truth'])['rows']
        prediction = read(f'{R}/predictions/{case["name"]}.json')
        predictions[case['name']] = prediction
        tests_for_guard.append(dict(name=case['name'],language=case['language'],mode=case['mode'],rows=observation['rows']))
        require(prediction['case']==case['name'] and prediction['test_sha256']==sha(case['test']) and
                prediction['model_sha256']==sha(case['model']) and prediction['truth_parsed'] is False,'Prediction header')
        require(len(prediction['rows'])==len(truth)==len(observation['rows'])==64,'Prediction row count')
        case_details = []
        for i,(test,actual,pred) in enumerate(zip(observation['rows'],truth,prediction['rows'])):
            row_id,block = f'{case["language"]}-{i:03d}','A' if i<32 else 'B'
            require(all(r['id']==row_id and r['block']==block for r in (test,actual,pred)),'Row identity')
            infos = [inverse.analyze(word,case['mode']) for word in test['observations']]
            support = all(plain in info['strings'] for plain,info in zip(actual['payloads'],infos))
            require(len(actual['payloads'])==len(infos),'Truth payload count')
            ambiguous = [info['distinct']>1 for info in infos]
            for method,order in (('trigram',2),('unigram',0)):
                answer = pred[method]
                independent_path_check(test['observations'],case['mode'],inverse,emissions,models[case['language']],order,answer)
                require(len(answer['payloads'])==len(infos),'Predicted payload count')
                correct = [a==b for a,b in zip(answer['payloads'],actual['payloads'])]
                case_details.append(dict(id=row_id,block=block,support=support,ambiguous_tokens=sum(ambiguous),tokens=len(infos),
                    method=method,edit_distance=distance(answer['plaintext'],actual['plaintext']),true_characters=len(actual['plaintext']),
                    exact=answer['plaintext']==actual['plaintext'],correct_tokens=sum(correct),
                    ambiguous_correct=sum(c and a for c,a in zip(correct,ambiguous)),states=answer['states'],transitions=answer['transitions']))
        for method,block in itertools.product(('trigram','unigram'),('A','B')):
            group = [d for d in case_details if d['method']==method and d['block']==block]
            require(len(group)==32,'Score block size')
            count = lambda field:sum(row[field] for row in group)
            edits,characters,exact,supported = map(count,('edit_distance','true_characters','exact','support'))
            primary = method=='trigram' and case['variant'].startswith('ordered')
            summary.append(dict(case=case['name'],language=case['language'],variant=case['variant'],mode=case['mode'],
                method=method,block=block,rows=32,support_rows=supported,edit_distance=edits,true_characters=characters,
                character_accuracy=1-edits/characters,exact_rows=exact,exact_row_accuracy=exact/32,
                tokens=count('tokens'),correct_tokens=count('correct_tokens'),ambiguous_tokens=count('ambiguous_tokens'),
                ambiguous_correct_tokens=count('ambiguous_correct'),ambiguous_rows=sum(bool(d['ambiguous_tokens']) for d in group),
                ambiguous_exact_rows=sum(d['exact'] for d in group if d['ambiguous_tokens']),primary=primary,
                gate_pass=supported==32 and edits<=40 and exact>=26,
                max_states=max(d['states'] for d in group),max_transitions=max(d['transitions'] for d in group)))
        details.append(dict(case=case['name'],rows=case_details))
        print(json.dumps(dict(stage='recovery',case=case['name'],status='PASS')),file=sys.stderr,flush=True)
    compare(details,read(R+'/row-scores.json'),'all row scores')
    primary = [s for s in summary if s['primary']]
    require(len(summary)==64 and len(primary)==24,'Gate inventory')
    final = dict(status='KNOWN_KEY_CALIBRATION_COMPLETE_NOT_A_TRANSLATION',
                 all_primary_gates_pass=all(s['gate_pass'] for s in primary),primary_gates=24,
                 primary_gates_passed=sum(s['gate_pass'] for s in primary),summary=summary,
                 observed_unique_mode_token_pairs=len(inverse.cache),original_word_boundaries_recovered=False,
                 manuscript_or_unknown_key_run=False)
    compare(final,read(R+'/summary.json'),'summary and gates')
    report['recovery'] = dict(status='PASS',prediction_conditions=16,prediction_method_rows=2048,summary_rows=64,
                              primary_gates=24,primary_gates_passed=final['primary_gates_passed'],
                              inverse_method='Complete finite string-set enumeration; no production DFA imported',
                              candidate_string_cap=STRING_CAP,path_score_tolerance=1e-8,
                              isolated_production_replay=isolation_replay(rows,models,tests_for_guard,predictions))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('prepared','complete'),required=True)
    parser.add_argument('--output')
    args = parser.parse_args()
    destination = ROOT/(args.output or ('results/audit-017-naibbe-inputs.json' if args.stage=='prepared' else
                                        'results/audit-017-naibbe-oracle.json'))
    require(not destination.exists(),'Refusing to overwrite an audit artifact')
    report = dict(status='RUNNING_NOT_COMPLETE',stage=args.stage,audit_sha256=sha('src/audit_naibbe_oracle.py'),
                  independent_reference_imports_production=False,
                  limits='Local inverse language is a deck-ignorant support superset. Viterbi is a path maximum, not a plaintext posterior.')
    code = 0
    try:
        prepared_result = inspect_prepared(report)
        if args.stage=='complete':
            inspect_recovery(report,prepared_result)
        report.update(status='PASS',counts=dict(COUNTS),maximum_errors=dict(ERRORS),
                      unique_frozen_files_checked=len(prepared_result[0].files))
    except Exception as error:
        report.update(status='FAIL_OR_INCOMPLETE',exception=type(error).__name__,reason=str(error),
                      counts=dict(COUNTS),maximum_errors=dict(ERRORS))
        traceback.print_exc();code = 1
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('x',encoding='utf-8') as handle:
        json.dump(report,handle,ensure_ascii=False,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:report[k] for k in ('status','stage','counts','reason') if k in report}),flush=True)
    sys.exit(code)


if __name__=='__main__':
    main()
