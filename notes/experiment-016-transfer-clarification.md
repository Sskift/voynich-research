# Independent-review clarification: what the transfer gate establishes

Read with experiment016-rules.md and report012. The registered transfer criterion
checks whether each KNOWN CIPHER CODE's contribution and the true overall expansion
obey the source-envelope inequalities. It does not assert that every overlapping
plaintext digraph's frequency is below the envelope.

A concrete archived example is Celsus/homophone127 TEST: its own maximum plaintext
digraph frequency is156/6397 =0.02438643, above the full-source sensitivity envelope
6315/275084 =0.02295662. Its full-source transfer still passes. That is possible
because the digraph's occurrences include starts across code-unit boundaries and
different homophones divide contributions among multiple cipher labels.

This does not contradict the registered gate or change any score. It shows why
'the known ciphers pass' must not become 'the envelope covers every plaintext
frequency in these works', much less 'the envelope covers all Latin'. Primary
block envelope and full-source sensitivity remain separate and unmodified.

Other retained conditions from the independent review:

- Each opaque reference position is assumed to denote exactly one unknown
  letter-unit. Expansion, nulls, variable missing spans or hidden word boundaries
  fall outside the conservative count calculation.
- Excluding unreadable manuscript words gives a theorem about the retained-word
  subset under its candidate plaintext frequencies. It is not a worst-case bound
  over every possible completion of the full manuscript.
- A two-letter code must emit within a true plaintext word for the within-word
  digraph ceiling to apply. Reinterpreting visible spaces alone does not remove
  the constraint if every code still emits wholly within true words.
- The upper-endpoint length witness satisfies both scalar inequalities only
  when lower<=upper, as clarified in experiment-016-report-correction.md.

The independent mathematical review found no algorithmic error changing016's
conditional incompatibility results. No phonetic or semantic values follow.
