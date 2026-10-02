"""Archive complete020 recovery tables and descriptive diagnostics; never decode."""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'results/omission-oracle-v2'
LANGUAGES = ('la', 'de', 'it')
MECHANISMS = ('ordered', 'iid', 'copy')
RULES = ('identity', 'omit-aeiou')
SPLITS = ('development', 'reserve')


def read(path):
    return json.loads((ROOT/path).read_text())


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def main():
    frozen = read('results/omission-oracle-v2/result-freeze.json')
    for path, expected in frozen['hashes'].items():
        assert sha(path) == expected, path
    summary = read('results/omission-oracle-v2/summary.json')
    assert summary['status'] == 'COMPLETE_ORACLE_CALIBRATION_NO_NATIVE_SCORES'
    rows, diagnostics = [], {}
    fields = ['known_words', 'exact_words', 'reference_characters',
              'char_edit_accuracy', 'exact_word_accuracy', 'lossy_words',
              'exact_lossy_word_accuracy', 'lossy_char_edit_accuracy',
              'nontrivial_words', 'nontrivial_exact_word_accuracy',
              'nontrivial_char_edit_accuracy', 'eligible_complete_rows',
              'exact_complete_rows', 'exact_complete_row_accuracy', 'excluded_rows',
              'eligible_row_word_coverage', 'nontrivial_confident_words',
              'nontrivial_confident_correct', 'nontrivial_confident_accuracy',
              'nontrivial_confident_coverage', 'omitted_symbol_edit_accuracy',
              'baseline_char_edit_accuracy', 'baseline_exact_word_accuracy',
              'ordered_recovery_thresholds_pass', 'high_confidence_error_flag',
              'identity_exact']
    for language in LANGUAGES:
        certs = read(f'results/omission-oracle-v2/{language}-completion-certificates.json')
        omission = {k: v for k, v in certs.items() if k.startswith('omit-aeiou|')}
        assert all(v['numerical_complete'] and v['search_complete'] for v in omission.values())
        notes = dict(unique_omission_skeletons=len(omission),
            maximum_log_mass_allowance=max(v['numerical_certificate']['conservative_log_mass_error_bound'] for v in omission.values()),
            maximum_search_expanded_prefixes=max(v['search_counters']['expanded_prefixes'] for v in omission.values()),
            oov_top2_numerical_ties=sum(v['search_counters'].get('top2_numerical_tie', False) for v in omission.values()),
            all_complete=True, cases={})
        for mechanism in MECHANISMS:
            for rule in RULES:
                case = f'{language}/{mechanism}-{rule}'
                path = f'results/omission-oracle-v2/scores/{case}.json'
                score = read(path)
                for split in SPLITS:
                    rows.append(dict(language=language, mechanism=mechanism, rule=rule,
                                     split=split, **{key: score['tests'][split][key] for key in fields}))
                if rule == 'identity':
                    continue
                pred = read(f'results/omission-oracle-v2/predictions/{case}.json')
                truth = read(f'data/omission-oracle/truth/{case}.json')
                errors = []
                for p, t in zip(pred['rows'], truth['rows']):
                    assert p['locus'] == t['locus']
                    for index, (word, original) in enumerate(zip(p['tokens'], t['tokens'])):
                        if original['units'] is None:
                            assert word is None
                            continue
                        text = ''.join(original['units'])
                        if word['prediction'] != text and word['posterior'] >= .9:
                            errors.append(dict(locus=p['locus'], split=p['split'], slot=index,
                                truth=text, prediction=word['prediction'], observation=word['observation'],
                                posterior=word['posterior']))
                errors.sort(key=lambda x: (-x['posterior'], x['locus'], x['slot']))
                allowance_sums = [sum(omission[pred['rows'][r]['tokens'][t]['completion_certificate_key']]
                    ['numerical_certificate']['conservative_log_mass_error_bound']
                    for r, t in segment['positions']) for segment in pred['segments']]
                notes['cases'][mechanism] = dict(
                    confidence_bins={s: score['tests'][s]['confidence_bins'] for s in SPLITS},
                    nontrivial_confidence_bins={s: score['tests'][s]['nontrivial_confidence_bins'] for s in SPLITS},
                    row_length_strata={s: score['tests'][s]['row_length_strata'] for s in SPLITS},
                    high_confidence_error_words=len(errors),
                    descriptive_top_five_errors=errors[:5],
                    error_example_selection='Posterior descending, then locus/slot; descriptive after scoring, never fed into decoding.',
                    maximum_segment_sum_of_log_mass_allowances=max(allowance_sums, default=0),
                    numerical_scope='Sum of recorded completion allowances only; not a complete HMM posterior interval certificate.')
        diagnostics[language] = notes
    assert len(rows) == 36
    target = ROOT/'results/omission-oracle-summary.csv'
    with target.open('x', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['language', 'mechanism', 'rule', 'split']+fields)
        writer.writeheader(); writer.writerows(rows)
    with (ROOT/'results/omission-oracle-diagnostics.json').open('x') as out:
        json.dump(dict(status='DESCRIPTIVE_SUMMARY_OF_FROZEN020_V2_OUTPUTS',
            result_freeze_sha256=sha('results/omission-oracle-v2/result-freeze.json'),
            summary=summary, languages=diagnostics,
            limitation='Known-key controls, not manuscript readings; complete computations do not establish recovery.'),
            out, ensure_ascii=False, indent=2, allow_nan=False)
        out.write('\n')
    print(json.dumps(dict(status='COMPLETE', rows=len(rows),
        omission_skeletons=sum(x['unique_omission_skeletons'] for x in diagnostics.values()),
        ordered_gate_pass=summary['ordered_oracle_gate_pass'])))


if __name__ == '__main__':
    main()
