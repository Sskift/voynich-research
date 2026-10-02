# Experiment016: necessary frequency bounds for fixed length1/2 units

Selected after015 calibrated and fitted the actual manuscript, before computing
this bound on reference corpora. No016 scores yet.015's top candidate expands
training text only1.103-fold and leaves the most frequent glyph as a monograph.
That observation motivates a solver-independent diagnostic; it is not a proof.
Keep all015 fits, failures and freezes. Do not tune the bound to the manuscript.

## Derivation to audit before use

For a complete cipher with fixed word boundaries, let f_c be each cipher symbol's
token frequency, and λ be decoded letters / cipher units. Every code emits one or
two letters, so λ=1+sum_{c:|K(c)|=2}f_c and1<=λ<=2.

For a given true plaintext, define a as its maximum single-letter count / total
letters, and b as its maximum WITHIN-WORD adjacent-letter-pair count / total
letters. Pair counts use all overlapping starts but the denominator is LETTERS,
not pair opportunities. Each occurrence of a digraph-emitting code contributes
one occurrence of that same adjacent pair, even if other occurrences straddle
code boundaries. Thus f_c<=λa for monograph codes and f_c<=λb for digraph codes.
Homophones only increase required plaintext counts, so these inequalities remain
necessary. Also b<=a, hence λ>=max_c f_c/a is necessary regardless of code length.

Given an upper bound t on λ, all f_c>t*b must be monographs. Therefore replace t
by min(t, 2−sum_{f_c>t*b}f_c), beginning at t=2. Iterate to a fixed point with an
explicit floating-point tolerance or exact rational arithmetic. This is a safe
upper bound, not a sufficient feasibility test. Reject an exact candidate corpus
only if the necessary lower bound exceeds that upper bound. Independently check
small alphabets by enumerating all1/2 length assignments and compare feasible λ.

## Reference scope and calibration before manuscript result

Use existing archived Latin sources without changing their normalization:
classical Caesar1–7, historical Macer/Aquinas. Retain unknown characters as opaque
and derive CONSERVATIVE upper bounds for a,b rather than silently discard them.
The exact handling must be written and frozen before source counts. Single-letter
upper bound can assign every opaque position to the most frequent known letter;
pair upper bound can add every within-word pair touching an opaque position to
the most frequent fully known pair. Retain resulting looseness and source defects.

Primary empirical envelope: maximum a and maximum b over all existing >=5000-unit
word-end blocks of these references; include full texts separately and report both.
Maxima may come from different blocks/works: this deliberately loose envelope is
NOT a natural language theorem. No empirical p-value from overlapping blocks or
same-work blocks. Do not use average frequency or exact majorization as a universal
language criterion;012 already demonstrated extensive cross-work false rejection.

Before actual use, regenerate all015 control ciphers and verify the true plaintext
counts/true λ satisfy the inequalities. Add exhaustive small-codebook checks with
homophones and non-overlapping encoded digraphs. Then score held-out Cicero/Celsus
compressed controls against the source-only envelope; report false rejections
without changing envelope. If this empirical transfer fails, the bound remains
an exact-corpus invariant only; no manuscript-language exclusion is warranted.

Apply a frozen, calibrated envelope to the same79H nativeGC training frequencies;
16H and16P remain development sensitivity checks. Do not consume f49v. Report
lower/upper λ, every forced-monograph iteration, dominant-symbol witnesses, and
unknown-position treatment. A pass provides no reading. A rejection only limits
compatibility with these source-frequency envelopes; unusual language/genre,
different unit sizes, contextual codes or different spaces remain outside scope.

This is a bounded diagnostic for whether more015 search is informative, not a
replacement for the overall goal of a connected, independently verified reading.
