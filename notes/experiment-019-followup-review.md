# After019: require continuous recovery even with a known key

Status: proposal only, NOT a frozen020 protocol. No new source counts, generated
ciphers, manuscript scores, or f49v text inspection were performed for this idea.
This is a proposed computational hypothesis, not a historical or linguistic claim.

019 did not identify a language or copying mechanism. Its negative prediction
gains also occurred on known ordinary text, while cipher channels obscured the
local-copy statistic. Repeating a weak classifier is not a route to a reading.
The next candidate should face a direct continuous-plaintext recovery test before
any manuscript search.

One specific candidate is a deterministic omission channel: within each plaintext
word remove the normalized symbols a/e/i/o/u, except that a word consisting only
of those symbols is retained unchanged. Then apply a fixed homophonic encoding,
initially retaining true word boundaries. This does not allow each manuscript
glyph an arbitrary longer output. Omitted symbols must instead be inferred from
the whole word and its context. The rule is only an invented, falsifiable cipher
channel; it is not asserted to be a historical abjad or shorthand practice. In
particular, i/u in the archived Latin normalization are not purely phonetic vowels.

Start with an oracle-key control: reveal the correct substitution key on the
existing three source pairs and test whether the recovery model can reconstruct
continuous held-out text after this loss of information. Report character edits,
whole words, exact multiword rows, coverage and ambiguity for every required
source pair. Retain alternative candidates and all failures. Separate ambiguity
introduced by the omission channel from limitations of the chosen source prior;
failure of a prior is not a proof that the channel is intrinsically undecodable.

Only if a prespecified continuous-recovery gate passes on every required control
should an unknown-key end-to-end solver be considered. Its own positive controls
must then pass before application to manuscript rows. Apply the same omission,
encoding and complete search to the019 iid/copy controls at the same budget, so
plausible-looking reconstructed sentences on those controls remain visible.
Correctly recovering the known unordered plaintext is a success, not a failed
negative control. The relevant failure is falsely repairing it into a different,
fluent, high-confidence sentence or treating recovery confidence as semantic proof.
Renaming an already-frozen key is only an additional, weaker negative control.
Do not repair words by hand or present a single attractive path as established text.

A future020 protocol must freeze the normalization, preserved boundaries, handling
of unknown slots and all-omitted words, model/training rules, search budget, exact
acceptance thresholds, output alternatives and test partitions before new counts.
No threshold is selected here. Existing manuscript pages remain development;
even successful control recovery would neither identify this mechanism in the
manuscript nor establish a translation.019 supplies no positive evidence for
this particular omission rule. It supplies a reason to test recoverability first.
