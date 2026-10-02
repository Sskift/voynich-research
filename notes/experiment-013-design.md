# Experiment013: fixed transcription components before a new decoding model

Selected after experiment012, before computing this experiment's matched-corpus
unit profiles or attempting an EVA-component cipher. Objective still undeciphered.

## Bounded model selection

Experiment012 found nativeGC training words shorter and its most frequent unit
more concentrated than six reference works. This does not exclude their languages.
The next explicitly selected candidate model is a fixed many-to-one, one-component
to-one-plaintext-letter mapping on an existing EVA transcription, with original
word boundaries retained. EVA components are an encoding convention, NOT known
phonemes or necessarily individual physical glyphs. Test its admissibility before
training any new language model or searching keys. Do not hand-split native units
to make a candidate plaintext, or tune new merges to the reference statistics.

## First stage: source-controlled unit comparison

- Reuse only experiment011's79H training loci for selection and the already-seen
  16H/16P bodies descriptively. No f1r or CurrierB global decoding scores.
- Compare GC2a-n with ZL3b-n and IT2a-n at the same folio/locus. RF is dependent
  combined material, not a new replicate. Preserve source missingness/disagreements.
- Count EVA case-preserving encoded atoms, high-ASCII escape as one atom, no
  semantic merges and no lowercasing rule invented after profile inspection.
  Parse IVTFF syntax, uncertainty and optional spacing explicitly; do not count
  markup or escape code digits as character positions. Establish parser rules
  from archived format/source headers before scores. If the format defines an
  uppercase rendering as multiple components, state that fact and show a separate
  declared sensitivity; do not silently equate it to lowercase.
- Primary matched-word comparison: same eligible locus in all three files, same
  number of dot/comma-separated word positions, and no uncertainty/alternative in
  any of its three words. Report how much this agreement restriction removes.
  No fuzzy alignment or lexical correction. Unequal word-count loci are excluded
  from paired word statistics but retained in per-source full summaries.
- Report word-token triplets, unit length ratios, alphabet counts, max frequency,
  source-specific 5000-unit blocks and same word-length/frequency diagnostics as012.
  Paired GC counts must be recomputed on the SAME accepted word sample; do not
  compare a filtered EVA subset to the entire GC population without labeling it.
- Report recurring GC-word/EVA-word relationships descriptively. They are not
  letter values, translations, or evidence for a language. Different transcriptions
  concern one manuscript and are not independent manuscript validation samples.
- Use the frozen012 reference corpora and both apostrophe conventions, with no
  new language/source/chunk selection. Freeze matched selection before scores.

## Continuation gate

Apply the same coarse screen (some reference block can admit the maximum frequency,
and mean word length lies in that language's reference-block range) to the PRIMARY
EVA representation, declared as ZL with IT sensitivity. Do not choose whichever
transcription happens to favor a candidate. The gate uses all legible ZL words at
the frozen79H training loci; the strictly matched triplet subset is a declared
sensitivity, not a replacement chosen after scoring. A pass only authorizes calibration of
this bounded decoding model, not a linguistic conclusion. If any primary candidate
passes, prepare a source-separated language model and known-cipher recovery controls;
freeze their configuration before actual manuscript attack, and APPLY it after
calibration succeeds. Declare validation pages and any distribution shift before
new decoding-score inspection. Prior79H/16H/16P pages are development data now.

If no primary candidate passes, preserve the full negative result and assess a
fixed variable-length code separately (one symbol emits one or two letters), with
an explicit length model and synthetic recovery controls. Do not add per-word
exceptions, guessed plant names or manual grammatical repairs. Do not claim that
the coarse screen mathematically excludes all alphabetic languages.

Independent audit must rebuild source alignment, count all exclusions, rederive
atoms and profile statistics, verify frozen012 references, and retain010/011 freezes.
No new plaintext output has been produced for this experiment at design creation.
