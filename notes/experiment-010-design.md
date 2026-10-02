# Experiment 010: calibrated, restricted global Latin substitution

Frozen before any manuscript decoded-output score. Goal remains undeciphered.

## Scope and representation

Native GC2a v101, Currier A, hand 1, illustration H only. Eligible loci are
@P0,+P0,*P0 except those followed by =/~ paragraph continuation. Page assignment:
SHA256(ASCII folio) as integer modulo 5 == 0 gives held-out pages; all others train.
This reserves decoding scores, not previously unobserved manuscript acquisition.
One global bijective key maps the 23 most frequent unambiguous TRAINING symbols
(frequency descending, lexical tie break) to abcdefghiklmnopqrstuxyz (j->i,v->u).
Preserve native case and @ddd; units. This is an explicitly restricted model:
less common symbols are masked in place, NOT merged or claimed to be nulls.
Words containing ? or [...] are wholly masked because their extent may be unclear.
Dot and comma both become space; braces merely delimit ligatures and are removed.
No letter sequence across line edges; one boundary space at each line end.
Unknown tokens invalidate any overlapping score window, never join neighbours.
All exclusions and positional masks are retained. No manuscript word deletion,
invented pronunciations, wordwise keys, or semantic crib values.

## Language model and calibration

Train a smoothed character four-gram model on Caesar, De bello Gallico books1–7,
Perseus Latin edition; extract text body only, exclude notes/headings/foreign text.
Lowercase, fold j/i and v/u, normalize diacritics, punctuation becomes space.
Objective: sum of log joint four-gram probabilities, additive smoothing0.1,
alphabet23 letters plus fixed space. Four-grams are a search heuristic, not a
test of grammatical meaning or a normalized generative model of whole sentences.
Cipher solver: C++ permutation simulated annealing plus exhaustive swap ascent.
Solver parameters may be calibrated ONLY on synthetic control texts, frozen
with source/code hashes before first manuscript training search.

Calibration corpus disjoint by author/work: Cicero In Catilinam (Perseus) and
Celsus De medicina book1 (LacusCurtius Latin, notes/header/footer excluded).
At least three random substitutions, including a medical-domain control and
a control with the manuscript training unknown-letter fraction randomly masked.
Train keys on the first up-to24000 characters of each selected control; reserve
the remaining text or an explicitly separate speech for evaluation. Chunk into
eight-word lines; score no cross-line windows. Report token accuracy (excluding
masks/spaces), supported-letter key accuracy and fullword recovery, per restart.
Require >=95% unmasked test-letter recovery for the selected key in each control
before interpreting the actual attack; failure calls for solver repair, not a
claim that manuscript is nonlinguistic. Never tune against manuscript held-out.

## Manuscript attack and comparisons

Use the frozen search configuration on training. Select best training objective
only, serialize/hash key before held-out score/decoding. Output complete lines
with unknown marks, plus explicit status NOT_A_TRANSLATION. Compare held-out
scores to999 random bijections and known-control correct-key scores.
Also fit9 independent null corpora: within each locus, permute all known glyphs
among known slots, preserving spaces/masks/line lengths and line glyph counts.
Use identical search effort and independently shuffled held-out loci per null.
The null changes word shape; it only measures low-order structure beyond such
shuffling and is NOT evidence for Latin when rejected. A Monte Carlo comparison
with9 nulls has minimum attainable p0.1; treat as descriptive screening.
Save ALL restart keys and scores; inspect whether multiple incompatible keys
produce similar training fits. Report model-coverage and held-out score limits.
No extra unit model, language, space regime or mapping is tuned this round.

## Interpretation and audit

Synthetic success validates this restricted solver on these sources; it does not
validate glyph assumptions, medieval orthography, an underlying Latin language,
or historical cipher design. Manuscript failure only limits this exact baseline.
An ngram gain, isolated Latin word or high significance versus random keys is
not a reading. Goal success still requires independently justified connected
plaintext and prediction of unused manuscript material. Independently recompute
key bijectivity, masks, training/held-out separation, objective, synthetic accuracy,
and frozen hashes. Keep failed calibration runs and design amendments if any.
