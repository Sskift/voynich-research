"""Independent infinite completion-mass audit for the stopped020 v1 case.

Only independent audit helpers are imported. A reachable product graph is built
explicitly; NumPy solves each vowel-core block and an iterative acyclic fill
handles other histories. OOV mass is a nonnegative first-prefix-exit sum, never
full mass minus lexicon mass. This does not certify the global OOV top-two search.
"""
from collections import defaultdict, deque
from pathlib import Path
import argparse
import datetime
import json
import math
import sys
import traceback

import numpy as np

from audit_omission_recovery import IndependentModel, Integrity


ROOT = Path(__file__).resolve().parents[1]
VOWELS = frozenset('aeiou')


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def sha(path):
    import hashlib
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def log_add(terms):
    terms = list(terms)
    if not terms:
        return -math.inf
    peak = max(terms)
    return peak+math.log(math.fsum(math.exp(term-peak) for term in terms))


def omit(word):
    return ''.join(c for c in word if c not in VOWELS) or word


class GraphMass:
    """Unconditioned remaining-completion mass at every reachable graph node."""
    def __init__(self, model, skeleton):
        need(skeleton and not any(c in VOWELS for c in skeleton), 'Expected nontrivial skeleton')
        self.model = model
        self.skeleton = skeleton
        self.vowels = tuple(c for c in model.alphabet if c in VOWELS)
        self.start = (0, (' ',))
        self.states = {self.start}
        self.by_position = defaultdict(list)
        self.edges = {}
        self.terminal = {}
        pending = deque([self.start])
        while pending:
            node = pending.popleft()
            position, history = node
            self.by_position[position].append(history)
            self.terminal[node] = self.log_next(history,' ') if position==len(skeleton) else None
            edges = []
            letters = self.vowels + ((skeleton[position],) if position<len(skeleton) else ())
            for letter in letters:
                target = (position+(letter not in VOWELS), self.shift(history,letter))
                edges.append((letter,target,self.log_next(history,letter)))
                if target not in self.states:
                    self.states.add(target); pending.append(target)
            self.edges[node] = edges
        self.logs = {}
        self.position_diagnostics = []
        self.solve()

    @staticmethod
    def shift(history, letter):
        return (history+(letter,))[-3:]

    def log_next(self, history, letter):
        coordinates = tuple(self.model.index[c] for c in history+(letter,))
        return float(self.model.tables[len(history)][coordinates])

    @staticmethod
    def trailing_vowels(history):
        count = 0
        for letter in reversed(history):
            if letter not in VOWELS:break
            count += 1
        return count

    def solve(self):
        for position in range(len(self.skeleton),-1,-1):
            histories = sorted(self.by_position[position])
            core = [h for h in histories if len(h)==3 and all(c in VOWELS for c in h)]
            need(len(core)==len(self.vowels)**3, 'Incomplete reachable vowel core')
            core_index = {h:i for i,h in enumerate(core)}
            vowel_matrix = np.zeros((len(core),len(core)))
            forcing = {}
            for history in histories:
                node = (position,history)
                if position==len(self.skeleton):
                    forcing[history] = self.terminal[node]
                else:
                    consonant_edges = [edge for edge in self.edges[node] if edge[0] not in VOWELS]
                    need(len(consonant_edges)==1, 'Incorrect consonant transition')
                    _,target,logweight = consonant_edges[0]
                    forcing[history] = logweight+self.logs[target]
                if history in core_index:
                    for letter,target,logweight in self.edges[node]:
                        if letter in VOWELS:
                            need(target[1] in core_index, 'Core vowel edge escaped')
                            vowel_matrix[core_index[history],core_index[target[1]]] += math.exp(logweight)
            # Position-specific scale protects arbitrarily small completion mass.
            scale = max(forcing.values())
            scaled_forcing = {h:math.exp(value-scale) for h,value in forcing.items()}
            need(all(value>0 and math.isfinite(value) for value in scaled_forcing.values()),
                 'Within-position forcing underflow; no mass claim')
            system = np.eye(len(core))-vowel_matrix
            rhs = np.array([scaled_forcing[h] for h in core])
            solved = np.linalg.solve(system,rhs)
            need(np.isfinite(solved).all() and np.all(solved>0), 'Invalid core solution')
            values = {h:float(solved[i]) for h,i in core_index.items()}
            outside = [h for h in histories if h not in core_index]
            # Each vowel edge raises this count until it enters the solved core.
            outside.sort(key=lambda h:(-self.trailing_vowels(h),h))
            for history in outside:
                terms = [scaled_forcing[history]]
                for letter,target,logweight in self.edges[(position,history)]:
                    if letter in VOWELS:
                        need(target[1] in values, 'Noncore vowel fill is not acyclic')
                        terms.append(math.exp(logweight)*values[target[1]])
                values[history] = math.fsum(terms)
            absolute_residual = relative_residual = 0.
            for history in histories:
                terms = [scaled_forcing[history]]
                for letter,target,logweight in self.edges[(position,history)]:
                    if letter in VOWELS:terms.append(math.exp(logweight)*values[target[1]])
                residual = abs(values[history]-math.fsum(terms))
                absolute_residual = max(absolute_residual,residual)
                relative_residual = max(relative_residual,residual/values[history])
                self.logs[(position,history)] = math.log(values[history])+scale
                need(self.logs[(position,history)]<=2e-10, 'Remaining event probability exceeds one')
            core_residual = float(np.max(np.abs(system@solved-rhs)))
            self.position_diagnostics.append(dict(position=position,reachable_histories=len(histories),
                vowel_core_states=len(core),noncore_states=len(outside),log_scale=scale,
                maximum_vowel_row_sum=float(vowel_matrix.sum(axis=1).max()),
                maximum_scaled_core_equation_residual=core_residual,
                maximum_scaled_all_equation_residual=absolute_residual,
                maximum_componentwise_relative_equation_residual=relative_residual,
                minimum_scaled_solution=min(values.values()),maximum_scaled_solution=max(values.values())))

    def full_log_mass(self):
        return self.logs[self.start]-self.model.log_nonempty

    def oov_log_mass(self, lexical):
        prefixes = {word[:n] for word in lexical for n in range(len(word)+1)} | {''}
        lexical = set(lexical)
        pending = [('',self.start,0.)]
        contributions = []
        visited = eos_terms = exit_terms = 0
        while pending:
            prefix,node,log_prefix = pending.pop()
            position,history = node
            visited += 1
            if position==len(self.skeleton) and prefix not in lexical:
                need(prefix, 'Empty word entered OOV mass')
                contributions.append(log_prefix+self.log_next(history,' ')); eos_terms += 1
            for letter,target,logweight in self.edges[node]:
                next_prefix = prefix+letter
                next_log = log_prefix+logweight
                if next_prefix in prefixes:
                    pending.append((next_prefix,target,next_log))
                else:
                    # Once a prefix is absent, no extension can be a lexicon word.
                    contributions.append(next_log+self.logs[target]); exit_terms += 1
        need(contributions, 'Nontrivial observation should have infinite OOV support')
        return (log_add(contributions)-self.model.log_nonempty,
                dict(lexical_words=len(lexical),prefix_nodes=len(prefixes),visited_prefix_nodes=visited,
                     first_exit_terms=exit_terms,nonlexical_prefix_eos_terms=eos_terms,
                     combination='Nonnegative first-prefix-exit/EOS terms; no full-minus-lex subtraction'))


def comparison(actual_log, archived_log):
    log_error = actual_log-archived_log
    actual = math.exp(actual_log)
    archived = math.exp(archived_log)
    need(abs(log_error)<=1e-8, f'Independent mass differs from archived estimate: log error {log_error}')
    return dict(independent_log_mass=actual_log,archived_log_mass=archived_log,
                independent_mass=actual,archived_mass=archived,
                absolute_mass_error=abs(actual-archived),relative_mass_error=abs(math.expm1(log_error)),
                absolute_log_error=abs(log_error),audit_log_agreement_tolerance=1e-8)


def audit_v1():
    integrity=Integrity()
    inputs=integrity.verify('data/omission-oracle/input-freeze.json')
    prepared=integrity.verify('data/omission-oracle/prepared-manifest.json')
    terminal=integrity.verify('results/omission-oracle/failed-v1-freeze.json')
    need(prepared['input_hashes']==inputs['hashes'], 'Input/prepared chain mismatch')
    failure_path='results/omission-oracle/computation-failure.json'
    need(failure_path in terminal['hashes'], 'Failure absent from terminal freeze')
    failure=read(failure_path)
    need(failure['status']=='COMPUTATION_INCOMPLETE_NO_GATE_OR_NATIVE' and
         failure['observation']=='ppl' and failure['run_context']['case']=='la/ordered-omit-aeiou',
         'This first audit is restricted to the stopped v1 ppl case')
    archived=failure['analysis']
    need(archived['numerical_complete'] is False and archived['search_complete'] is True,
         'Original failed numerical status changed')
    model=IndependentModel('la',integrity)
    skeleton='ppl'
    lexical=sorted(word for word in model.words if omit(word)==skeleton)
    need(lexical==archived['lexicon_candidates'], 'Archived lexical candidate set incomplete')
    graph=GraphMass(model,skeleton)
    full=graph.full_log_mass()
    oov,oov_accounting=graph.oov_log_mass(lexical)
    lex_logs={word:model.char_log(word) for word in lexical}
    partition=log_add([oov]+list(lex_logs.values()))
    need(abs(partition-full)<=1e-10, 'Independent full/OOV/lexicon partition disagrees')
    concrete=[]
    for label,items in [('lexical',[dict(word=w,log_char_probability=archived['lexicon_log_char_probabilities'][w]) for w in lexical]),
                        ('archived_top_oov',archived['top_oov'])]:
        for item in items:
            word=item['word']
            need(omit(word)==skeleton and (word in model.words)==(label=='lexical'), 'Concrete candidate classification')
            log_probability=model.char_log(word)
            need(abs(log_probability-item['log_char_probability'])<=1e-10, 'Concrete word probability disagreement')
            concrete.append(dict(kind=label,word=word,independent_log_probability=log_probability,
                                 archived_log_probability=item['log_char_probability'],
                                 absolute_log_error=abs(log_probability-item['log_char_probability'])))
    numerical=archived['numerical_certificate']
    return dict(status='V1_SINGLE_FAILURE_MASS_AUDIT_PASS_NOT_CALIBRATION',skeleton=skeleton,
                case=failure['run_context']['case'],language='la',unique_skeletons_audited=1,
                full_mass=comparison(full,archived['log_char_mass']),
                oov_mass=comparison(oov,archived['log_oov_char_mass']),
                independent_partition_log_error=abs(partition-full),lexical_candidates=lexical,
                concrete_word_probabilities=concrete,reachable_graph_states=len(graph.states),
                reachable_graph_edges=sum(map(len,graph.edges.values())),
                position_equation_diagnostics=graph.position_diagnostics,oov_accounting=oov_accounting,
                original_numerical_failure=dict(numerical_complete=False,
                    archived_conservative_log_mass_error_bound=numerical['conservative_log_mass_error_bound'],
                    unchanged_tolerance=numerical['numerical_log_tolerance'],
                    conclusion='Original run remains stopped; agreement of estimates does not retroactively validate its failed bound.'),
                scope=dict(independent_algorithm='Explicit reachable graph; per-position numpy.linalg.solve vowel core; iterative noncore fill; forward prefix-exit OOV sum.',
                    probability_input='Frozen018 tables independently normalized; complete prepared training lexicon.',
                    production_modules_imported=False,independent_audit_helper='src/audit_omission_recovery.py',
                    infinite_mass='Full and lexicon-OOV mass independently recomputed with unbounded vowel insertion.',
                    residuals='Residuals of represented scaled equations, not rigorous interval or forward-error certificates.',
                    top2_search='Only probabilities of the archived concrete words checked; global top2 search completeness is not claimed.',
                    v2_execution_or_inspection=False,manuscript_scoring=False,calibration_claim=False),
                hash_manifests=integrity.manifests,distinct_frozen_files=len(integrity.files),
                failure_sha256=sha(failure_path),input_freeze_sha256=sha('data/omission-oracle/input-freeze.json'),
                audit_source_sha256=sha('src/audit_omission_mass.py'),
                independent_helper_sha256=sha('src/audit_omission_recovery.py'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='results/audit-016-omission-mass-v1.json')
    args=parser.parse_args()
    path=ROOT/args.output
    need(not path.exists(),'Refusing to overwrite an audit; choose another --output')
    try:
        report=audit_v1()
    except Exception as error:
        report=dict(status='MASS_AUDIT_FAILED',error_type=type(error).__name__,error=str(error),
                    traceback=traceback.format_exc(),audit_source_sha256=sha('src/audit_omission_mass.py'))
        stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        failed=path.with_name(path.stem+'.failure-'+stamp+'.json')
        failed.parent.mkdir(parents=True,exist_ok=True)
        with failed.open('x',encoding='utf-8') as handle:json.dump(report,handle,ensure_ascii=False,indent=2)
        print(json.dumps(report,ensure_ascii=False,indent=2),file=sys.stderr)
        raise
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as handle:
        json.dump(report,handle,ensure_ascii=False,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({key:report[key] for key in ('status','skeleton','full_mass','oov_mass',
                     'independent_partition_log_error','reachable_graph_states','reachable_graph_edges')},indent=2))


if __name__=='__main__':
    main()
