# Zodiac label order and repeated sequence constraints

2026-10-02. Previous goal turn classified PROGRESS after inspecting authoritative
STATE, rule outputs, and audit-003. This experiment is specified after reading
layout comments and checking transcription coverage, before any order statistic.
No phonetic value, calendar date, or degree number is assigned.

Inventory all 299 ZL Lz labels, matching IT by page/locus. A primary label may
contain several words, but every token must be known and identical between ZL
and IT (split uncertain spaces). Keep a fixed missing/disagreement mask; never
close gaps. Coverage check: 186 agreed labels, 277 fully known ZL labels. This
coverage check is not an order test. Retain . between words in similarity scores.

Build explicit clockwise groups from IVTFF descriptions and confirm layout from
original scans/catalogue. Do not sort by approximate clock strings, which contain
errors. Pisces has 19 outer +10 inner human labels and one central star label;
Gemini has 5 top-arc +15 outer +9 inner labels, with an unlabelled outer figure
between loci19 and20. Insert a missing slot there. No artificial cycle connects
separate rings or the two 15-label panels of a zodiac sign. Top arcs do not wrap.

Discovery diagrams: all zodiac pages on folios70–72. Confirmation diagrams:
f73r/f73v (Scorpio/Sagittarius), reserved from score inspection until discovery
outputs have been saved. The text is in the same downloaded corpus; this is a
prospective computational split, not blind acquisition or independent scribes.
All diagrams are assigned the same scribal hand in the corpus metadata.

Primary question: do diagrams with the SAME ring-size layout share an ordered
label sequence? Compare six pairs among the four 10+5-label panels, plus Leo vs
Virgo (18+12). Confirmation compares Scorpio vs Sagittarius (16+10, ignoring
the separate four-item top arcs). Within each corresponding ring, allow every
cyclic shift and reversal, but no arbitrary reorder or variable skip. Similarity
= 1 - Levenshtein distance / maximum string length, in EVA including preserved
word separators. Each phase needs at least3 jointly accepted labels. Report
phase support; exclude ring pairs with no eligible phase. Aggregate = mean best
similarity over eligible corresponding-ring pairs, equal ring-pair weights.
It is an order-coherence screen, not an identification of the ordered items.

Null: independently shuffle accepted strings within each ring, keeping unknown
positions fixed. In every randomization re-optimize shift/reversal. Reuse each
shuffled ring for every comparison involving it, preserving dependence between
pair scores. 9,999 permutations, seed20261007. As a planned layout sensitivity,
shuffle only within ring AND EVA letter-length bins <=4,5–6,7–8,>=9, preserving
spatial length clustering. State when few positions remain exchangeable.
A same-order exact-equality score is secondary; it cannot identify meaning.

Second question: are immediate clockwise neighbours more similar than this same
within-ring shuffle? Use complete circular groups only, with wraparound where
physically present. Unknown/missing slots break edges. Pool all valid adjacent
edges, weighting each edge equally. Rings with fewer than3 accepted labels do
not contribute. Arcs and the central Pisces star are excluded. A positive result
could reflect naming structure, copying, or layout, not necessarily ordinals.

The primary order score and neighbour score form a two-test family; report Holm
adjustments within each stage/null/processing run. Also report every local ring
and pair score descriptively, without promoting the strongest one post hoc.
Use the same primary inclusion loci for IT split and ZL joined-spacing checks;
do not select a new best semantic interpretation from their differences.

A result earns further investigation only if it survives the layout sensitivity
and appears in the preassigned confirmation diagrams. Even that is not a
translation: it must predict independently grounded referents or readable text.

Implementation note before scoring: under the agreed split-label criterion,
retained IT split strings are identical to ZL split by construction. Verify that
identity instead of redundantly simulating the same data and presenting it as
independent replication. ZL join remains a separate spacing sensitivity.
