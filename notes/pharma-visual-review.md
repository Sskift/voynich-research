# Visual review of labels and proposed plant correspondences

2026-10-02. Exploratory inspection by the assistant, with transcription and
proposed correspondence visible. Not blinded, not a botanist's identification,
not a new catalogue of species. Shape descriptions below are observations of
these particular drawings, not claims about real plants. The available mirror
scans are modest resolution; edges, folds, and overpainting limit interpretation.

## Photograph coverage

All 12 downloaded pharmaceutical photographs were viewed. The explicit mapping
in `data/pharma-image-coverage.json` replaces the original, erroneous assumption
that a mirror filename is an IVTFF page ID. The archive preserves original pixels
and filenames. Notable corrections:

- `088v.jpg` includes f88v, f89r1, f89r2 from left to right.
- `089r.jpg` is the foreground f89v2, not f89r1/r2.
- `101r.jpg` contains the left half of logical f101v (old f101v2).
- `101v.jpg` contains its right half at left, then f102r1 and f102r2.
- `102r.jpg` is the foreground f102v2; `102v.jpg` is f102v1.

The [quire 15](https://www.voynich.nu/q15/index.html) and
[quire 19](https://www.voynich.nu/q19/index.html) catalogues provide the foldout
layout. Their linked thumbnails for f89r1 and f101v2 were independently viewed
against the mirror photographs. Label inventory image references are now fixed;
text-only statistics did not depend on photograph filenames.

## All six repeated agreed label types

Selection preceded image review. Each locus is retained even when object
association is uncertain. These descriptions do not establish a one-to-one
mapping between word types and depicted species.

| EVA | Loci | Visible comparison and limitation |
|---|---|---|
| okol | f88r.15; f101v.10 | First is beside an elongated green, toothed form in the middle group, with another label near its roots; second is beside a broad curled green leaf in the left panel's second row. Both may describe a part; no distinctive identical object established. |
| otaly | f88r.6; f99v.6 | First accompanies pale, long ribbon-like roots with two green leaves; second a brown root crown with several curved roots and a cut/curved stalk. The drawings are not copies; this cannot exclude different preparations or parts of one plant. |
| otoky | f88r.5; f99v.9 | First is beside a thick brown root mass with lateral strands; the second is at the folded right edge of the top row. Its object is insufficiently visible to judge a match or mismatch. |
| otaldy | f88r.12; f101v.2 | First is catalogued as a container label, between the middle jar and the adjacent root group; second is above the left panel's second plant, with tiers of downcurving leaves and a small root. Container-vs-plant classification itself prevents assuming two identical naming contexts. |
| otoldy | f89r1.4; f99r.14; f99v.2 | These accompany respectively a fine-rooted top-row plant, a broad unpainted ribbon-root form at the lower right of the top group, and a narrow red central root with many short side branches. No visually distinctive identical drawing established. |
| oky | f99r.4; f99r.24 | Same-page repetition: a narrow upright root in the top group and an extended horizontal root with long pointed projections in the middle group. A common category remains possible. |

## Five labelled external picture correspondences

The full selection contains nine source-proposed identities; four have no
individual plant label and remain in `data/herbal-image-pairs.json` as exclusions.
Only the five labelled pairs below were reviewed against their separate herbal
scans in this round. The other four target scans were downloaded, not yet visually
assessed. All five labelled pharmaceutical objects were visible in the reviewed
photographs; no botanical names were assigned.

| Fragment and target | Observed shared drawing feature | Limitation |
|---|---|---|
| 54 on f89v2 → f48v | Broad leaves with deep, irregular marginal lobes | Fragment has two leaves, target many and flowers; a leaf motif alone is not an identity test. |
| 61 on f89v2 → f48r | Rounded green forms with long projections at one end | Whole plant has additional unpainted flower heads; still only a proposed relation. |
| 212 on f102r2 → f18v | Long horizontal tan rhizome with knobby upright projections and a curling end | Fragment is mainly a root; the page adds leaves/flower. Strong visual lead, no confirmed species. |
| 213 on f102r2 → f23r | Horizontal segmented root/rhizome with fine dangling roots | Small foldout drawing lacks the whole plant's leaves and flowers. |
| 240 on f102v1 → f19r | Truncated root crown and long curved leaf with downward finger-like divisions | Fragment has one prominent leaf, whole plant two and several additional structures. |

The fifth item at f89v2.28 has ZL `daseky` versus IT `dareky`. It is excluded
from the agreed-word primary test, and both readings are shown in sensitivity
output. Source proposals and visible likeness constrain a test; neither supplies
a plaintext translation. The source author's original selection method is not
fully established here, so these are not independent blind evidence.
