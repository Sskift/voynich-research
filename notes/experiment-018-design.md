# Experiment018: infer plaintext spaces as ordinary cipher outputs

Protocol before new source-plus-space frequency counts or key scores.017rejected
all72prespecified1/2-letter, fixed-space conditions. This motivates testing another
specific family; it does not justify repairing its old decoded words.

## Representation and candidate languages

Primary input is nativeGC, with case and @nnn; retained. The archived original
transcription explanation (data/raw/transcr.html) describes its intended complete
glyph units; EVA is a component notation. This fixes GC before any new gate or
decoding, not by selecting whichever score looks better. It is still not a proved
linguistic unit. Use the same frozen013loci:79Htrain,16H/16Pdevelopment. No f49v.

Every known glyph maps to exactly ONE output character: a letter OR an ordinary
SPACE. Multiple glyphs may emit the same character. No NULLs, multiletter outputs,
contextual plaintext values, reordered text, compressed spaces or dictionary repairs.
Observed visual word gaps are ignored by the model, with their positions retained
in metadata. Unknown whole words and unseen test glyphs break the context; never
join known runs across them. Physical locus boundaries also break context without
forcing a plaintext space or sentence ending. IVTFFfile continuations are not new
physical loci. Preserve every predicted space, including repeated or edge spaces.

Candidates are Latin, historical German and Italian. The old014German gate used
fixed observed spaces and is not reused for language choice here.

- Latin LM: archived Caesar1–7; alphabet abcdefghiklmnopqrstuxyz plus SPACE.
  Gate references Caesar/Macer/Aquinas as016; control work CiceroCatilinarians1–4.
- German LM and gate references:014v2Brant base projection only; alphabet
  abcdefghijklmnopqrstuvwxyzæß plus SPACE. Control work014v2Tepl base projection.
  The global projection stays frozen, including loss of historical combining marks.
- Italian LM and gate references:Dante primary012normalization; alphabeta–z plus
  SPACE. Control workBoccaccio primary012biography. No posthocaccent/letterchanges.

Control works are EXCLUDED from the gate reference pool. Gate and solver checks
are cross-work within each language, but these archived works have all been seen
in previous experiments. Three encodings of one work are not independent works.

## Cheap source-only gate, calibrated separately from decoding

Keep the existing word-end >=5000LETTER-unit block boundaries from016(Latin) or
012(Brant/Dante). For Brant apply014v2baseletter projection using the SAME original
word indices, not new block boundaries. For Nprojectedletters and Wwords in a
block insert exactly W−1SPACE characters. M=N+W−1. If C is maximum known letter
count and U is unknown letter positions, use max(C+U,W−1)/M. Unknown positions
cannot add to the SPACEcount. Keep whole-work profiles as a named sensitivity.

Primary language ceiling is the maximum of these block character-frequency bounds.
Compare nativeGCtraining's maximum glyph count/N to it using exact fractions.
Every one-character code must meet this necessary condition, including homophones.
This gate is deliberately weak: it does not certify a feasible space assignment,
joint frequency vector, valid syntax or a language. Do not add exact majorization.

Generate the known ciphers below before using this screen for native interpretation.
Check their unmasked full train/test streams against their OWN exact plaintext
character frequency, then report whether the source-only reference ceiling admits
them. If a language fails empirical gate transfer, retain the failure and do not
make a dependable language exclusion; do not widen its ceiling after seeing it.

Fit EVERY language whose native primary gate admits AND whose own new decoder
calibration passes. Do not choose a language by comparing raw LM scores. A failed
gate can be reported only as incompatibility with those numerical assumptions.

## Probability model and budget

For each source record, concatenate its existing words with single spaces, with
no leading/trailing spaces. Latinrecords are archivedparagraphlines; Germanrecords
and Italianverses remain their saved record boundaries. Train counts of all ngrams
length1..4 within records. Unigram+.1; recursive Dirichletbackoff mass10. Save all
four conditional probability arrays. SPACE is an ordinary member of each array.

For each known cipher segment, predict every output character using order
min(4,position+1); initial positions use the lower-order models. Every cipher unit
is predicted once. Add the channel likelihood sum Nc logNc−sum Nt logNt, i.e.
trainingMLE homophone emission. If a generated character has no cipher code its
probability is a failure outcome, not a silent null. Divide by all input units,
including those decoded as spaces. No division by variable non-space letter count.

16restarts,3cycles each,25000proposals/cycle,temperatures N*.02→N*.00005;
at most100coordinatepasses. Moves80%singlevaluechange(halfweightedbycipherfrequency),
15%swap two glyph values,5%globallyswap two output-character names. SPACE participates
in every move. Initial greedy frequency allocation uses source unigram INCLUDING
SPACE. Restart0 keeps it; others randomly rename all output characters, and even
restarts independently replace30%glyphvalues uniformly. Seed20261102pluslanguage
index*10000, indicesLatin0/German1/Italian2, with separate frozen offsets for
encoding, fitting, actualsearch, nulls. Select training likelihood only. Retain keys.
Pythoninitialization uses numpydefault_rng(solver_seed+100000), before C++search.

## New known-cipher controls with hidden plaintext spaces

For each control work take the first complete-word span reaching24000characters
(letters plus spaces) as training, then the next reaching20000 as test, or all
remaining words if shorter. No word overlap. Encode the intervening TRUEspaces
as unknown cipher symbols just like letters. Cut physical rows at60cipherpositions,
including inside words. False visible groups have uniform random lengths2..8
within each row and carry no plaintext boundary meaning.
The final visual group is truncated to remaining row length and can have length1.
If the last physical row contains just one character, the tailcondition uses its
first-position rule (first takes precedence over last).

Three conditions per language:
1. bijection: one code per output character, SPACEunknown too.
2. uniform127:127codes allocated one per output character plus extras proportional
   to source unigram counts, largestremainders with alphabet-order ties. SPACE
   therefore gets multiple codes. Uniform choice within each homophone group.
3. tail127-position-mask: same allocation; interiorpositions choose primary.97,
   otherwise uniformly extra codes. At first/last physical-row position choose
   respectively first extra/last extra with probability.97, primary otherwise.
   Groups with one code use it. Mask whole false-visible groups independently.02.
   This retains globally fixed outputs while making glyph choice position-dependent.

Use numpydefault_rng with seedbase+conditionindex for randomlabels/emissions/groups/
masks in that order. Keep true labels and plaintext ONLY in calibration files;
solver input has cipher segments, aggregate counts, model, initials and budget.
Map only training-seen labels to solver IDs; unseen test labels break context.
Record masked/unseen coverage separately. Invariant test: changing false visual
groups while keeping symbols/masks fixed leaves parsed segments and scores identical.

Every condition must achieve >=.95 predictable-character accuracy AND >=.95
non-space letter accuracy AND SPACEprecision and recall both>=.95. Report keytype
accuracy, coverage and complete-word recovery too; rare mappings can remain weak.
All four recovery criteria apply to predictable positions; masked/unseen positions
are excluded from SPACErecall as well as letteraccuracy, with their coverage
reported separately. These conditional accuracies are not full-stream recovery.
Words qualify only if their true bracketing spaces and all letters are predictable
within the same physical row. A correct word requires letters AND both boundaries
correct, with no interior predicted space. Never select restart by these metrics.
Failures kept; any solver/objective revision versioned before more scores.

## Actual fits, diagnostics and integrity

After successful per-language calibration freeze code/model/input/config/controls
before actualGCfit. Every admitted language gets the same budget. Report all train,
16H and16Pdevelopment character likelihoods, raw predicted spaces, sourceword hits,
key disagreement and complete outputs markedNOT_A_TRANSLATION. Never infer success
from isolated words. Do not use future f49v to choose between these fits.

For each fitted language perform9equal-budget fits on independently shuffled known
glyphs within each physical row, preserving unknown positions and row lengths.
Evaluate corresponding shuffled16Pdevelopment. Add999global output-character-name
permutations, INCLUDING SPACE. Neither null identifies a language. Describe scores
within their own model only; cross-language raw scores have different entropies.

Independently reconstruct input segments, LMarrays, allcontrolrecoveries, allfit
objectives and nulls. Existing freezes remain unchanged. A solver, gate or audit
passing does not fulfill the overall goal of a reproducible connected reading.
