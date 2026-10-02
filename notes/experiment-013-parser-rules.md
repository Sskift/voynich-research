# Parser rules fixed before matched profiles

Original archived IVTFF_format.txt sections6.6/6.7 define period as an apparent
space, comma as uncertain space, <-> and <~> as drawing interruptions implying
a space. Primary counting splits all four, as specified in experiment013. One
declared sensitivity deletes commas (no split); other spaces remain. The primary
gate is not selected from this sensitivity. Retain all counts/exclusions.

Curly braces mark connected characters; count their explicitly encoded components
without counting braces. Square brackets denote uncertain alternatives, '?' one
unreadable character and '???' an unknown-length gap; exclude each affected whole
word from full known-word profiles and paired slots. Do not collapse an unreadable
word to one known glyph. High-ASCII @128;..@255; is one atom; malformed codes raise
an error rather than silently contribute digits. Source text and word positions
remain saved, including unknown words. Inline comments/tags are metadata; drawing
interruptions are replaced by boundaries before stripping other inline metadata.

Original transcr.html describes EVA as a shape encoding, not semantic units. It
explicitly contrasts sequences 'in'/'iin' with single units in older systems.
Capitalization marks connections, not an extra symbol replacing multiple letters.
Keep case and ligature components unchanged. Syntax inspection at the frozen
training/development/reserve locus set found no uppercase EVA atoms; therefore
no capitalization sensitivity is needed. An apostrophe within ligatures such as
{c'y} is a shape component, not a word apostrophe: retain it as a literal atom.
Do not apply the historical-plaintext apostrophe rules to manuscript transcription.

Eligible locus must exist in the relevant source, have @P0/+P0/*P0 type, and not be
followed by a =P*/~P* continuation in that source. Scope is the prior1388GC locus
IDs, not newly selected pages. Record missing/ineligible source loci individually.
For paired slots, all three source loci must be eligible and have equal word counts
under the same spacing rule. Every retained slot must be fully legible in all three.
Equal counts are an alignment convention, NOT a guarantee of matching physical word
boundaries: compensating split/join disagreements may remain. Describe word triplets
as positional correspondences, not paleographically verified character mappings.

GC counts on the complete frozen locus set should reproduce012's23721training
units and6586known words; compare012's other splits with its unreadable-word rule.
Reference corpora/normalization/blocks remain frozen012. No new profile frequencies
or language-screen results inspected before saving these parser rules.
