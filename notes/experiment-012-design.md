# Experiment012: source-controlled language and unit constraints

Freeze before computing new-language frequency or word-length comparison scores.
Earlier011 is completed progress, not decipherment. No goal completion is claimed.

## Declared source set

Three bounded candidates: medieval Latin, historical German, historical Italian.
This is not exhaustive language identification. Two different works per candidate:
- Latin: Macer Floridus, De viribus herbarum (Latin column of Remacle's historical
  edition transcription); Thomas Aquinas, Summa theologiae Ia question1 (Corpus
  Thomisticum's Leonine1888 text). Expand to further fixed adjacent questions only
  if q1 has fewer than10000 usable letters, documenting the reason before scores.
- German: Brant Narrenschiff1494, Narragonien digital's edited reading; Johannes
  von Tepl, Ackermann, Knieschek1877 edition on Gutenberg47465, MAIN text only.
- Italian: Dante Commedia, Gutenberg1012; Boccaccio's Vita di Dante within vol1
  of Guerri's edition, Gutenberg21424, MAIN prose only, exclude quoted poem blocks.

Brant and Locher1497 diplomatic transcriptions are supplementary GLYPH-reference
conditions, not additional independent works or canonical alphabetic Latin.
Locher's downloaded reg branches are EMPTY; do not invent an expanded Latin text.
Locher translates/adapts Brant, so no cross-language independence claim for these.
Download/provenance/access failures and editorial rules retained. Do not use
automatic catalogue summaries, translated English, modern editorial apparatus,
web navigation, or unverified OCR as historical-language samples.

## Extraction and units

Inspect markup before selecting exact body boundaries. Save extracted passages,
source locators, extraction rules and hashes before computing comparative scores.
Brant primary: mainText in numbered B chapters, reg and corr/lem branches only;
exclude paratext P chapters, notes, marginal references, apparatus alternatives.
This avoids the Latin prose paratext. Diplomatic companion uses orig/sic branches.
Locher: mainText in B chapters, orig branches, no notes/marginal glosses.
Keep raw transcription and normalized analysis tokens as separate fields.

Primary reading-unit normalization: lowercase; NFC grapheme clusters (base plus
combining marks). Italian accent marks folded to base; German marks preserved.
Latin j/i and v/u folded; ae/oe ligatures expanded to two letters. Standard
typographical ligatures expanded for reading-text conditions, long-s->s and
round-r->r where present. Unknown MUFI private glyphs retained as opaque units,
never deleted/guessed. Report opaque mass and a conservative maximum-frequency
upper bound that could assign all opaque units to the most frequent clear letter.
Diplomatic companions retain ligatures, abbreviation signs and combining clusters
as encoded grapheme units; no phonemic claims. Hyphenated line continuations are
rejoined only where explicitly marked by line-end hyphen or XML break convention.
Words follow spaces/punctuation. Apostrophes join adjacent segments in PRIMARY;
one prespecified SPLIT sensitivity for all sources, not selected after scores.
No word deletion for unusual glyphs. Nonletter isolated punctuation/numbers and
printed line/page numbers excluded, with source-specific rules recorded.

NativeGC comparison is previously scored79H A/hand1 training pages, same127
case-preserving units, high-ASCII escape as one token, dot/comma split. Also show
already-inspected16H development and16P bodies descriptively, not new validation.
Do not access global decoding outputs for unscored f1r/CurrierB pages this round.

## Fixed diagnostics

One ciphertext symbol->one plaintext letter, possibly many-to-one, implies:
1. maximal output-letter frequency >= maximal input-symbol frequency;
2. sorted output-frequency cumulative sums majorize input cumulative sums;
3. word length is unchanged if spaces/word segmentation are preserved.
These are exact for an actual paired plaintext/ciphertext. DIFFERENT works only
provide empirical references; a mismatch is not a proof excluding a whole language.

Partition each work into consecutive nonoverlapping blocks of at least5000 letter
units, stopping at the first word boundary at/above5000; discard the final shorter
remainder for block summaries but retain it in full-work statistics. Repeat the
same block procedure on native text. Report every block, full-work frequencies,
max clear/opaque frequency, ranked cumulative curve and word-length distribution.
For each source block versus native H full distribution report max positive
cumulative deficit D=max_k(sum_top_k(q_cipher)-sum_top_k(p_reference),0), plus
word-length KS distance and mean length. Padding by zero keeps total mass1.
No ad hoc corpus selection or normalization to minimize D or KS.

Calibration: generate random bijections and fixed many-to-one homophonic ciphers
from each source's first two complete blocks (up to four homophones per plaintext
unit, independent selections). Verify exact majorization, max bound and word-length
invariance against their OWN plaintext. Then measure cross-work reference variation
for same-language blocks, recording how a naive exact D==0 criterion can reject
genuine same-language text. No p-values from overlapping samples or tiny language
counts. These diagnostics are constraints/screens, not a posterior language label.

## Decision and audit

Declare empirical support only when at least one source block admits the native
top-frequency necessary condition AND native mean word length lies within the
combined within-language block range (including opaque uncertainty as reported).
This is a coarse, preregistered screen, not a hypothesis test or proof of encoding.
If a language meets these minimum constraints, follow with calibrated actual
decoding under that explicitly defined alphabet; do not call the screen a decipherment.
If none does, do not run a knowingly unsupported one-letter model just to generate
more pseudo-words: quantify what word-length/unit assumptions would need to change,
using diplomatic/reading pairs where available to assess realistic contractions.
Choose the next bounded decoding model from these measured limitations.

Independently recompute source-selection integrity, unit counts, block assignments,
all diagnostic values and synthetic invariant controls. Audit cannot authenticate
all source spellings or establish the manuscript's language/meaning. Freeze source
and design hashes; preserve previous experiment010/011 freezes unchanged.
