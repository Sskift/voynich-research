# Clarification of the inherited “79H” name

The frozen019 protocol inherits the name “79H” from earlier experiments. That
name describes the candidate page cohort, not the number of pages with retained
rows. A post-preparation bookkeeping check found:

| Split | Candidate pages from GC metadata and fixed split rule | Pages with selected rows | Rows | Known word slots | Unknown word slots |
|---|---:|---:|---:|---:|---:|
| train | 79 | 78 | 980 | 6,586 | 0 |
| development | 16 | 16 | 185 | 1,242 | 0 |
| reserve (section transfer) | inherited P cohort | 16 | 223 | 2,449 | 31 |

The difference is **f49v**: its header places it in the candidate training cohort,
but the old frozen selection contains no rows from that page. No019 selected row
comes from it, and no019 score uses it. The full archived GC source remains a
frozen dependency. This check compared metadata headers, the fixed page
hash rule in `src/substitution.py`, and existing selected rows; it did not score
or inspect that page's text to change the selection.

The actual training rows are identical in `data/substitution-manuscript.json`,
`data/homophonic/manuscript.json`, `data/transcription-units/selection.json` and
the new019 template. The previously reported980rows/23,721known atoms and all
numerical results remain unchanged. Descriptions that called these “79 used
training pages” should read **78 used pages from a79candidate-page cohort**.

This is a scope clarification, not a changed split or retrospective protocol
change. All old frozen artifacts remain unchanged. Actual page coverage, rather
than the inherited cohort name, is used in the new report. f49v remains unscored.
