# Experiment015: bounded one-or-two-letter outputs, before score design and calibration

Selected after completed013/014. Not implemented or scored yet. Fixed one-letter
Latin substitution on GC and historical German substitution on EVA failed to give
a connected reading despite successful synthetic recovery.013 showed that known
representation changes alone can move the coarse language screen; none is phonetic
proof. Do not repeat completed simple models or repair their decoded words.

## Candidate family

Return to nativeGC, with case and escaped glyphs preserved. Each observed symbol
gets one globally fixed nonempty plaintext string of length1 or2. Several symbols
may emit the same string. Spaces remain fixed. No context-dependent value, per-word
exception, null output, reordering or hand-written crib. This bounds flexibility
while allowing the mean plaintext word length to exceed nativeGC word length.
Do not automatically split EVA 'ii' into repeated plaintext letters under this model.

The next concrete work is a correctly specified variable-length objective and a
known-cipher encoder/decoder benchmark, followed by an ACTUAL manuscript attack if
calibration supports this family. Mean word length alone does not support any key.
Begin in Latin using the already archived Caesar training and separate Cicero/Celsus
control works, whose one-letter controls were recoverable. Later medieval Latin
transfer must remain separate; do not silently replace a source after manuscript fit.
The language choice is a bounded starting point, not identification or an exclusion
of German/Italian. Keep all prior failed experiments and source limitations.

## Resolve score-length bias before fitting

Unnormalized sum log P(text) usually favors shorter decoded strings, while naive
per-letter averaging can favor spurious repetitive expansions. Derive and record
the probability model or explicitly named optimization objective BEFORE fitting
manuscript keys. In particular, explore a word-language probability conditioned on
output length with its normalization computed by finite-state forward sums, plus
an explicit length/segmentation prior. If a heuristic is used instead, label it as
such and test its selection bias on known competing outputs. No default shortcut
that pretends variable-length scores are directly comparable.

The initial derivation is a development step, not a frozen numerical solver. Once
chosen, save the exact objective, priors, alphabet, search budget and source splits
before calibration scores; preserve any calibration failures and revisions. Avoid
hand-setting each glyph's expansion from apparent manuscript words.

## Required calibration and actual use

- Create known plaintexts encoded with globally consistent one/two-letter code
  units, including both unit lengths and skewed/homophonic emission. Save actual
  segmentation and true values; do not pass them to the solver. Keep separate
  plaintext training/test spans and source-separated language-model training.
- Include a pure one-letter control so the new family cannot obtain apparent gains
  by gratuitous expansion. Include length-prior/normalization sanity checks whose
  alternatives differ in output length. Report normalized edit/letter accuracy,
  word accuracy, key-string accuracy, output length and coverage, with a criterion
  fixed before solver results. Recovery of only the short-output case is inadequate.
- Freeze successful solver/configuration then apply it to actual GC input. Record
  every restart and choose only by the training objective. Use matched-budget
  shuffled controls and full decoded output, explicitly NOT_A_TRANSLATION unless
  an independently supported continuous interpretation is established.

## Validation scope

All79H training pages,16H development pages,16P bodies and now f1r have been globally
score-inspected. They are development material, including for a new cipher family.
Do not call f1r unused merely because this future model is different. Before any new
confirmatory scoring, inventory remaining metadata/loci and freeze page scope;
preserve Currier/hand/section differences. If no same-stratum untouched pages exist,
label actual output exploratory, and any other-stratum prediction a distribution
shift rather than an equivalent replication. Existing photographs/transcriptions
are not themselves newly acquired evidence.

Goal remains the full decipherment, requiring reproducible connected readings and
successful predictions beyond fitted text. Positive controls, statistical structure,
new code and green audits are progress, not completion.
