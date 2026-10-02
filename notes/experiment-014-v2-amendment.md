# Version2 amendment, known-cipher development only

Version1 completed all three control searches and failed the declared95% threshold:
94.2964%,94.3060%,94.2789%. No manuscript German keys were searched. Its complete
artifacts are retained with a separate failed-v1 freeze. Do not treat this calibration
failure as evidence against the manuscript's language.

In the bijective control,896/1126 wrong predicted letters were u->y; z->t contributed
209. The historical Brant reading stores many u/o/a characters with superscript
letters as private-use glyphs. Masking them removed contexts containing these base
letters. This motivates a new, explicitly different NORMALIZED-letter model.

Official MUFI v4 search evidence and the downloaded academic TITUS character tables
agree on four definitions: E42C=a with e above; E644=o with e above; E72B=u with e
above; E72D=u with o above. TITUS tables at0061.htm,006f.htm,0075.htm are archived.
MUFI's PDF failed direct download and its character page returned a browser check;
neither failed page is treated as successfully archived character content.

For BOTH Brant and Tepl in v2, apply one uniform base-letter projection: convert
those four verified private codes to a/o/u/u and remove Unicode combining marks
after NFD. Preserve æ and ß as single units (no new ae,ss,sz expansion), and preserve
all other letters/spelling, word positions and order. Output alphabet is a-z,æ,ß.
All resulting source units must fall in this28-unit alphabet; no residual masks
may be silently deleted. Save projected words as new files;012's originals and
original constraint results remain unchanged. This is deliberately lossy graphemic
normalization, not historical orthographic recovery or a Voynich phonetic claim.

New control conditions:28-code bijection,64-code uniform,63-code97% primary with2%
mask; otherwise same budget/objective95% criterion. Seed20261024. Tepl text and its
60/40 split were already used to diagnose v1, so v2 calibration is DEVELOPED on
these known controls, not fresh independent validation. No source words, grammar,
or dictionary entries from Tepl are added to the LM. Source-normalization is global,
fixed before v2 scores and follows character definitions rather than per-word edits.

Keep f1r score reserve and original ZL units exactly as v1. After normalized
calibration passes, freeze v2 then execute the actual attack and same controls.
If it fails, preserve failure and identify whether objective/model or search is
responsible on known ciphers before proceeding. No cutoff reduction is permitted.
