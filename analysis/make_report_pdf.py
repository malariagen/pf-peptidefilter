"""Render the collaborator report to PDF (print layout, light theme, paginated).

Usage:  python make_report_pdf.py analysis/report.html analysis/output/R2-sensitivity-report.pdf

WeasyPrint is NOT in the repository's requirements.txt and must not be added there: that file
is byte-compared against live ``main`` by verify_provenance.sh, and it is the app's dependency
set, not the analysis's. The renderer is pinned in ``analysis/requirements-report.txt`` instead.

Pagination is version-dependent, so if page numbers are ever cited, pin the version from
requirements-report.txt and confirm the rendered file with pdfinfo.
"""
import os
import re
import sys

from weasyprint import HTML, CSS

SRC = sys.argv[1] if len(sys.argv) > 1 else "report.html"
DST = sys.argv[2] if len(sys.argv) > 2 else "report.pdf"

html = open(SRC).read()


def tag_wide_tables(doc, min_cols=8):
    """Mark any .tablebox whose header has many columns; those get landscape pages."""
    out, pos, n = [], 0, 0
    for m in re.finditer(r'<div class="tablebox">', doc):
        end = doc.find("</div>", m.end())
        block = doc[m.end():end]
        cols = len(re.findall(r"<th", block.split("</tr>")[0]))
        out.append(doc[pos:m.start()])
        if cols >= min_cols:
            out.append('<div class="tablebox wide">')
            n += 1
        else:
            out.append(m.group(0))
        pos = m.end()
    out.append(doc[pos:])
    print(f"  {n} wide table(s) set landscape")
    return "".join(out)


html = tag_wide_tables(html, min_cols=99)  # nothing is wide enough to need it now

# the artifact host supplies the document skeleton; supply our own for print
title = re.search(r"<title>(.*?)</title>", html)
title = title.group(1) if title else "Report"

PRINT_CSS = """
@page {
  size: A4; margin: 17mm 15mm 16mm 15mm;
  @bottom-center { content: counter(page) " / " counter(pages);
                   font-family: sans-serif; font-size: 8pt; color: #7c8896; }
}
@page wide { size: A4 landscape; margin: 14mm; }
/* the widest tables get their own landscape page; no forced breaks around them --
   those made WeasyPrint overlay the following block */
/* a named landscape page does not itself widen the containing block, so the
   wide tables need the landscape content width (297mm - 2x14mm) set explicitly,
   and must escape the max-width:100% that constrains their flex siblings */
section > .tablebox.wide { page: wide; width: 265mm; max-width: none; }
html { background: #fff; }
body { background: #fff !important; font-size: 9.6pt; line-height: 1.55; }
/* ---------------------------------------------------------------------------
   Flex/grid containers do not fragment across pages in WeasyPrint: whatever
   does not fit on the current page is dropped silently rather than carried to
   the next one. That is what removed five rows from the 17-row completeness
   table in section 02 and the whole of section 05 from the first PDF. Paged
   media therefore gets block layout, with the `gap` values restated as
   margins so the spacing is unchanged.
   Verified by re-render: all 17 rows and all six sections now appear.
   --------------------------------------------------------------------------- */
.wrap, header, section, figure, .callout, .decision, .kv, ul, ol {
  display: block !important;
}
.wrap > * { margin: 0 0 26pt; }
.wrap > *:last-child { margin-bottom: 0; }
section > *, header > *, figure > *, .callout > *, .decision > * { margin: 0 0 12pt; }
section > *:last-child, header > *:last-child, figure > *:last-child,
.callout > *:last-child, .decision > *:last-child { margin-bottom: 0; }
.shead { display: block !important; }
.shead .snum { margin-right: 10pt; }
.shead h2 { display: inline; }
.kv dt { font-weight: 600; color: #495566; margin: 10pt 0 2pt; }
.kv dt:first-child { margin-top: 0; }
.kv dd { margin: 0 0 0 14pt; }
li { margin-bottom: 6pt; }
.meta span { margin-right: 20pt; }
/* flex items default to min-width:auto, so one wide table stretches the whole
   column and everything else (the figure included) then overflows with it */
.wrap { max-width: 100%; width: 100%; padding: 0; gap: 26pt; }
.wrap > *, section > *, header > * { min-width: 0; max-width: 100%; }
.tablebox { overflow: hidden; }
.tablebox, figure, .draft, .callout, .decision { max-width: 100%; }
/* sticky headers repeat oddly in paged media */
thead th { position: static !important; }
.tablebox { overflow: visible !important; max-width: 100% !important; }
/* print has no horizontal scroll: wrap cells so wide tables fit the page */
/* auto layout + universally wrappable cells: columns size to content and the
   table's minimum width stays under the page width, so nothing clips or overlaps */
/* fixed layout is what actually forces width:100% to be honoured; it only
   overlapped before because overflow-wrap:anywhere is a no-op in WeasyPrint --
   with word-break:break-word the cells wrap instead */
table { font-size: 7.2pt !important; width: 100% !important;
        max-width: 100% !important; table-layout: auto !important; }
.tablebox.wide table { font-size: 7pt; }
/* WeasyPrint places stylesheets passed here BEFORE the document's own <style>,
   so rules that must beat the screen CSS need !important */
th, td { padding: 3.5pt 5pt !important; white-space: normal !important;
         overflow-wrap: break-word; word-break: break-word; }
/* only data cells hold their line; headers must be free to wrap or they set
   a minimum column width the page cannot honour */
td.num, td.mono { white-space: nowrap !important; }
th { white-space: normal !important; overflow-wrap: break-word; word-break: break-word; }
th.num { text-align: right; }
h1 { font-size: 22pt; }
h2 { font-size: 14pt; }
h3 { font-size: 10.5pt; }
section { break-inside: auto; }
.shead { break-after: avoid; }
h3 { break-after: avoid; }
/* tall blocks must be allowed to split, or they get pushed whole and leave a blank page */
.draft, .decision, .tablebox, table, section { break-inside: auto; }
figure, .callout { break-inside: avoid; }
tr { break-inside: avoid; }
thead { display: table-header-group; }
tfoot { display: table-footer-group; }
.tablebox { overflow: visible; }
.note-under { break-before: avoid; }
figure img { max-width: 100%; }
footer { break-before: avoid; }
a { text-decoration: none; }
"""

doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
       f'<title>{title}</title></head><body>{html}</body></html>')

HTML(string=doc, base_url=os.path.dirname(os.path.abspath(SRC)) or ".") \
    .write_pdf(DST, stylesheets=[CSS(string=PRINT_CSS)])
print(f"wrote {DST}")
