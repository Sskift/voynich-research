# Experiment 001 — layout before phonetics

Recorded before computing corpus statistics, 2026-10-02.

The objective remains a verifiable decipherment of the manuscript. This experiment
is a constraint-building step, not a translation. EVA letters encode shapes, not
known sounds. Existing work already describes layout effects; replication is not
a claim of novelty.

Sources: ZL3b and IT2a from René Zandbergen's collection. RF1b-e is a derived
reference, not an independent third reading. Preserve raw files and SHA-256.

Questions:

1. How strongly do p/f and gallows-initial words depend on paragraph position?
2. Does final m depend on line-final position within the same physical line?
3. Could m be a positional spelling of a common suffix? Compare predefined
   replacements: empty, r, l, n, y, in, iin, iiin. No pronunciation is assigned.
4. Does exact word repetition remain above a shuffle preserving each line's
   vocabulary? If yes, independent drawing of words within a line is inadequate.

Primary processing: split certain and uncertain spaces; discard ambiguous tokens
from numerical tests without closing gaps or redefining their positions. Lowercase
EVA and remove ligature braces, retaining original token and locus. Exclude rare
characters, illegible characters, alternative readings, and other non-letters.
Also rerun ZL with uncertain spaces joined. No text from comments enters counts.

Restrict layout tests to P0 lines with locator @, +, or *; these mark ordinary
left-justified text rather than label lists, circular text, or right-hand line
fragments. Preserve paragraph start/end explicitly from <%> and <$> and verify
balanced marks. No inference from file newlines alone.

Position tests: paired contrasts within each line; paragraph contrasts within
pages with both first and subsequent paragraph lines. Randomization nulls retain
line/page composition. A low p-value rejects only the stated exchangeability null,
not every language/cipher/hoax family. Confidence intervals resample pages, not
individual words. Use seed 20261002 and 1999 randomizations.

Suffix experiment: train vocabulary using eligible NON-final words on training
folios only. Folio split uses SHA-256 of the numeric leaf ID: hash integer modulo
3 equals 0 for held-out leaves (recto and verso remain together). Report all suffix
candidates' held-out coverage and smoothed unigram score, including leave-m-alone.
This is a dictionary-compatibility screen, not a valid proof of allography or
meaning. Common short outputs can win for unrelated reasons. A promising candidate
must later predict contexts and illustrated referents, and survive matched controls.

No successful-decipherment claim without a fixed decoding algorithm and mapping,
readable independently verified plaintext on held-out pages, consistency across
scribes/sections, explicit treatment of exceptions, and independently grounded
image/text or historical anchors. An attractive isolated word does not suffice.
