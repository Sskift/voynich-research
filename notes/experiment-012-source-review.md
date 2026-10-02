# Source review before comparison scores

This records markup and transcription-quality decisions made before experiment012
frequency/length comparisons. The language candidates and diagnostics are unchanged.

- Macer: lang=LA is also applied to French. Select LEFT column of each two-cell
  table row, exclude bold headings. This gives2269 Latin verse records, first
  Herbarum quasdam... and last una diagridii..., not the translation. The digital
  transcription credits Thierry Vebr; source index lists Choulant1832 and
  Baudet/Panckoucke1845, not Nisard. No claim of error-free or newly proofread text.
  Visible defects `f1uxum` and `sciasique~medetur` remain in the saved text: embedded
  digit/tilde between letters is one OPAQUE unit and does not silently disappear or
  split a word. This follows the design's rule to preserve unusual glyphs; no
  conjectural repair. Occasional middle dots count as punctuation. Report this
  source limitation and opaque upper bound. This is an edited website transcription,
  not new unverified OCR generated from a scan for this experiment.
- Aquinas:66 titled paragraphs of question1 alone exceed10000 units. Exclude all
  span.ref section locators; the Leonine1888 edition is editorially normalized.
- Tepl: div#textdiv main paragraphs after first chapter;87 critical apparatus
  paragraphs in blockquote and90 navigation/empty paragraphs removed. All34
  chapters remain. Remove bracketed Z. line markers and footnote callouts. Retain
  the edition's bracketed additions to historical text; this is a critical reading,
  not a diplomatic medieval manuscript. Explicit printed line-end hyphens rejoin.
- Boccaccio: first full Vita only, p ids00020..00174; exclude two stand-alone Latin
  poem blocks00084,00142;95 prose paragraphs remain. Later compendia, Guerri's
  editorial warning and commentary excluded. Source title-page edition is Bari,
  Laterza1918; that is NOT the composition date. Short inline authorial Latin
  quotations remain. No modern AI-generated catalogue summary used.
- Dante:14233 verses in100 cantos, exact first/last verse boundaries; headings,
  special-character table and PG boilerplate omitted. Embedded authorial Latin
  and other language passages are part of this work, not separately selected data.
  PG source edition unidentified; no claim of original fourteenth-century spelling.
- Brant:112 B chapters, reading and diplomatic branches of the same302 mainText
  blocks. Some app/lem elements are nested in rdgGrp: recursively select lem once.
  P paratext excluded. Reading still contains four private-use glyph types;
  diplomatic is a glyph-level comparison, not another independent language sample.
- Locher:110 B chapters,1482 mainText blocks. All reg branches empty; select orig
  only. Ampersand and Tironian et, if present, count as opaque abbreviation/logogram
  units rather than punctuation. They are never silently expanded to Latin letters.
- NFC grapheme unit = one base plus following combining marks. Unknown private-use
  characters and literal anomalous embedded digits/tilde or &/⁊ are opaque. All opaque
  mass can conservatively be assigned to the most frequent clear unit in the upper
  bound, which does not assume that expansion preserves length.
- Native training and development have no '?' words. The descriptive pharmaceutical
  split has31 whole unreadable-word masks ('?' can stand for more than one glyph).
  Exclude those31 words from its length/frequency denominator; do not pretend each
  unreadable word is one glyph. Preserve otherwise legible previously unseen glyphs
  in descriptive native counts. Original79H training constraint is unaffected.

No word-frequency or comparative length results were inspected in making these
source-boundary/encoding decisions. Source-freeze.json is regenerated only now,
before scoring, to incorporate the extraction fixes and this review.
