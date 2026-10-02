# Experiment016 frozen operational rules, before source-pair counts

Companion to the already archived016design. Previous turn was PROGRESS:015 was
actually calibrated/fitted/audited. Its nonsensical output is not a plaintext crib.
No016 corpus digraph frequencies or envelope results have yet been computed.

## Units, sources and blocks

References: Caesar1–7 from data/latin/caesar-gallic-N.txt concatenated in numeric
book order; Macer/Aquinas primary `records[*].words` in data/language-screen.
Caesar words become lists of characters. The other two already contain lists of
grapheme strings: retain these atomic units, including any combining marks. Use
their existing normalization and primary apostrophe treatment without change.
No pair crosses a word boundary. Denominators count all letter-unit positions,
not number of available pairs. Whole-source summaries retained separately.

Blocks accumulate complete words until >=5000 units, then begin a new block.
Final shorter tails remain in full-source summaries but not the block envelope.
This is identical to012 for Macer/Aquinas; Caesar blocks are newly constructed
by that same rule. No pairs cross words or blocks. Source boundaries separate
Macer/Aquinas/Caesar. No source or cutoff selection after results.

Opaque unit rule: any code point with Unicode categoryCo, any digit, or~,&,⁊
makes that atomic unit opaque, as in012. Each opaque position stands for ONE
unknown letter-unit for this screen; expansions of opaque editorial glyphs are
outside this normalization assumption. Let N be all positions, U opaque positions,
C the maximum fully clear single-unit count, P the maximum fully clear pair count,
and V the count of within-word pair positions touching any opaque unit. Then use
a=(C+U)/N and b=min((P+V)/N,a). The clamp is valid because every true pair count
is <= one of its constituent letter counts. Save all counts and ties. Do not
drop opaque words or turn them into guessed text.

Primary envelope takes maximum a and maximum b separately over ALL reference
blocks. Full-source envelope is a prespecified sensitivity. Report each maximum's
witness source/block. These independent maxima deliberately loosen the screen.
Neither is asserted to cover every possible Latin text, genre or historical hand.

## Exact computation and checks

Use fractions.Fraction throughout the inequality iteration. Starting upper t=2,
force every code with f_c>t*b to length1, replace t by min(t,2−their total mass),
repeat until unchanged. Lower=max(1,max(f_c)/a). Exact equality is admitted.
Archive every iteration, forced codes and rational endpoints; reject iff lower>t.
The relaxation may admit an impossible key; never interpret admission as recovery.

Exhaustive length test:1..4 code types, each count1..3; a=i/12 for i1..12,
b=j/12 for j0..i. Enumerate every length assignment. Every truly feasible length
must lie in the reported interval; every bound rejection must have no feasible
assignment. Keep false-admissions as evidence that this is necessary, not sufficient.

Exact synthetic cipher test: enumerate all code sequences of length1..4 for
four fixed codebooks [a,b], [a,b,aa,ab,ba,bb], the latter plus a second a-code,
and [a,b,ab]. Decode each sequence to one word, measure its exact plaintext a,b,
and assert its true expansion and each code satisfy the derived necessary bounds.

Regenerate all015 controls from original plain_train/plain_test and archived
source-defined codebook rules/seeds, including the labels unseen in training.
Those true labels are available to a known-cipher invariant test, so do not apply
015's decoder-coverage mask. Both full train/test spans must pass using their OWN
plaintext frequencies and true lengths. Save differences in coverage explicitly.

## Transfer gate and actual application

Before native frequency application, test the source-only block envelope against
all three full TEST spans. Primary transfer PASS requires true expansion inside
the interval AND every known monograph/digraph code obeying its own a/b inequality
in all three. Also report interval-admission alone, since it is a weaker check.
The first two tests share Cicero; Celsus is a second work, not a third.

Prespecified size sensitivity: split those same test words into >=5000 plaintext
letter-unit blocks, repeat transfer checks and report every failure. These blocks
are not independent replications. A block failure limits transfer to that scale,
even if the full-span primary gate passes. Never widen envelope to repair failure.
Repeat full-test transfer for full-source envelope as sensitivity, not reselection.

Then compute the necessary interval for the same79H nativeGC training set and
16H/16P development sensitivities under both envelopes. Count ORIGINAL nativeGC
atoms, including training-unseen but legible symbols. Exclude entire words with
literal? (the old source parser's whole-word unreadability mask). Retain case and
@nnn; codes as one unit. Record omitted words. Training must remain23721 units.
No015decoder key is used. Preserve f49v/f58r/f58v as unscored.

If transfer fails, report manuscript intervals only as compatibility with the
specified numerical envelopes; no dependable Latin-language exclusion follows.
Even successful finite transfer is empirical and cannot exclude all Latin.
No statistical p-values from deterministic bounds or shared text blocks.

Freeze protocol/code/input checksums before counting reference pairs. Freeze
reference envelope before control transfer, then freeze control results before
native application. Independently reconstruct counts, bounds and all decisions.
