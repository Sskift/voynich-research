# Experiment 002 — context as a check on the m hypothesis

Exploratory follow-up after experiment 001, 2026-10-02. The earlier held-out
vocabulary results have already been inspected; this is NOT a new pristine
confirmation set. Context statistics have not yet been inspected.

Experiment 001 found robust final-m enrichment but ambiguous lexical restoration:
r, l, and iin all achieve substantial held-out vocabulary coverage. Run a fixed
backward bigram model: target = current word; predictor = previous word. Training
uses non-final targets only, from the same training leaves as experiment 001.
No pair crosses a locus boundary or an unreadable token. Test is line-final m
with a known preceding word, from held-out leaves.

For each replacement, compute log2 P(target), log2 P(target | previous), and
context lift (the difference). Use additive 0.5 unigram smoothing with an unknown
bucket; bigram posterior predictive (count(prev,target) + tau * P(target)) /
(count(prev) + tau), with tau=10, 50, 200. Use actual glyph readings unchanged
for the preceding word. Candidate outputs are treated as unknown if absent from
the training vocabulary. Only descriptive comparisons, no success p-value from
selecting the best of several candidates.

Repeat for final r, l, n, y as calibration controls; applying a rule to these is
not assumed valid. Report m type-level examples and Currier-stratified scores.
If context rankings vary with smoothing, tokenization, or transcription, keep
the original m and retain competing hypotheses. Even stable predictive gains
would require semantic anchors to establish a decipherment.
