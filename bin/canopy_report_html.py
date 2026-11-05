#!/usr/bin/env python3
"""Render a Canopy Markdown report as one self-contained HTML file next to it.

  canopy_report_html.py <report>.md [...]

Used for <sample>_assembly_summary.md and pipeline_tools.md. Handles the Markdown those files use:
headings, tables, quotes (a quote starting with [!IMPORTANT] becomes a highlighted callout), fenced
code, rules, bulleted and numbered lists, images on a line of their own, `code`, bold, [links](url),
<https://autolinks>, and paragraphs that are wholly italic. Underscores inside a line stay literal,
because file and rule names use them. Images are embedded, so the HTML needs no other file; a
missing image, an unbalanced bold or code marker, or a malformed table stops the conversion rather
than producing a wrong report.
"""
import base64
import html
import mimetypes
import os
import re
import sys

CSS = """
:root { --bg:#fff; --fg:#1f2933; --muted:#5f6c7b; --line:#d9dfe5; --head:#f1f4f7; --note:#f3f7f4;
        --accent:#2f6f4f; --code:#eef1f4; --callout:#fff6e0; --callout-line:#c58a00; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#14181c; --fg:#e4e8ec; --muted:#9aa6b2; --line:#333b44; --head:#1e242a; --note:#1a231e;
          --accent:#6fbf94; --code:#252c33; --callout:#2a2414; --callout-line:#e0a82e; }
}
* { box-sizing:border-box; }
body { margin:0 auto; max-width:60rem; padding:2rem 1rem 3rem; background:var(--bg); color:var(--fg);
       font:16px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }
h1 { font-size:1.8rem; margin:0 0 .25rem; }
h2 { font-size:1.3rem; margin:2.25rem 0 .75rem; padding-bottom:.3rem; border-bottom:1px solid var(--line); }
h3 { font-size:1.05rem; margin:1.5rem 0 .5rem; }
nav { margin:1.25rem 0 0; padding:.75rem 1rem; background:var(--head); border-radius:6px; font-size:.92rem; }
nav ol { margin:0; padding:0; list-style:none; }
nav li { margin:.15rem 0; }
a { color:var(--accent); }
.gen, footer { color:var(--muted); font-size:.88rem; }
footer { margin-top:2.5rem; }
code { font-family:"Courier New",Courier,monospace; font-size:.92em; background:var(--code);
       padding:.05em .3em; border-radius:3px; overflow-wrap:anywhere; }
pre code { padding:0; background:none; overflow-wrap:normal; }
.table { overflow-x:auto; margin:.75rem 0 1rem; }
table { border-collapse:collapse; font-size:.9rem; }
th, td { padding:.35rem .6rem; border:1px solid var(--line); text-align:left; vertical-align:top; }
th { background:var(--head); }
td.num, th.num { text-align:right; font-variant-numeric:tabular-nums; white-space:nowrap; }
blockquote { margin:.75rem 0; padding:.6rem .9rem; background:var(--note); border-left:3px solid var(--accent);
             border-radius:0 6px 6px 0; }
blockquote.callout { background:var(--callout); border-left:5px solid var(--callout-line); }
pre { margin:.75rem 0; padding:.7rem .9rem; overflow-x:auto; background:var(--code); border-radius:6px;
      font-size:.85rem; tab-size:8; }
ul, ol { margin:.6rem 0; padding-left:1.5rem; }
li { margin:.3rem 0; }
figure { margin:1rem 0; }
figure img { display:block; max-width:100%; height:auto; border:1px solid var(--line); border-radius:4px;
             background:#fff; }
figcaption { color:var(--muted); font-size:.88rem; margin-top:.35rem; }
hr { margin:2.25rem 0 0; border:0; border-top:1px solid var(--line); }
p { margin:.6rem 0; }
@media print {
  body { max-width:none; padding:0; font-size:10.5pt; }
  nav { display:none; }
  .table { overflow:visible; }
  h2, h3 { break-after:avoid; }
  table, blockquote, pre, figure { break-inside:avoid; }
}
"""

NUMERIC = re.compile(r"^[-+]?[\d,]*\.?\d+\s*(%|bp|kb|Mb|Gb|x)?$")
LIST_ITEM = re.compile(r"(-|\d+\.) +(.*)")
IMAGE = re.compile(r"!\[([^\]]*)\]\(([^)\s]+)\)")
# A line that starts a block other than a paragraph.
BLOCK_START = re.compile(r"(#{1,3} |```|\||>|-{3,}\s*$|- |\d+\. |!\[)")


def inline(text, source):
    if text.count("**") % 2:
        raise SystemExit(f"ERROR: {source}: unbalanced bold marker in: {text[:80]}")
    if text.count("`") % 2:
        raise SystemExit(f"ERROR: {source}: unbalanced code marker in: {text[:80]}")
    # Code spans first, so nothing inside them is read as bold or a link.
    codes = []

    def keep(m):
        codes.append(f"<code>{html.escape(m.group(1), quote=False)}</code>")
        return f"\x00{len(codes) - 1}\x00"

    text = html.escape(re.sub(r"`([^`]+)`", keep, text), quote=False)
    text = re.sub(r"\[([^\]]+)\]\((https?://[^)\s]+)\)", r'<a href="\2">\1</a>', text)
    text = re.sub(r"&lt;(https?://[^\s&]+)&gt;", r'<a href="\1">\1</a>', text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    return re.sub(r"\x00(\d+)\x00", lambda m: codes[int(m.group(1))], text)


def paragraph(text, source):
    if len(text) > 2 and text.startswith("_") and text.endswith("_"):
        return f'<p class="gen"><em>{inline(text[1:-1], source)}</em></p>'
    return f"<p>{inline(text, source)}</p>"


def table_row(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def table(rows, source):
    width = len(rows[0])
    for row in rows[1:]:
        if len(row) != width:
            raise SystemExit(f"ERROR: {source}: table row has {len(row)} cells, expected {width}: {row[0]}")
    # A column is right-aligned only when every value in it is a number, so a column never mixes
    # left- and right-aligned values.
    numeric = [i > 0 and len(rows) > 1 and all(NUMERIC.match(r[i]) for r in rows[1:]) for i in range(width)]
    out = ['<div class="table"><table><thead><tr>']
    for i, cell in enumerate(rows[0]):
        cls = ' class="num"' if numeric[i] else ""
        out.append(f'<th scope="col"{cls}>{inline(cell, source)}</th>')
    out.append("</tr></thead><tbody>")
    for row in rows[1:]:
        out.append("<tr>")
        for i, cell in enumerate(row):
            tag = "th" if i == 0 else "td"
            scope = ' scope="row"' if i == 0 else ""
            cls = ' class="num"' if numeric[i] else ""
            out.append(f"<{tag}{scope}{cls}>{inline(cell, source)}</{tag}>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def figure(alt, path, base, source):
    """An image embedded as a data URI, linked to the file for full resolution."""
    if re.match(r"[a-z]+:", path) or os.path.isabs(path):
        raise SystemExit(f"ERROR: {source}: image path must be relative to the report: {path}")
    full = os.path.normpath(os.path.join(base, path))
    if not os.path.isfile(full):
        raise SystemExit(f"ERROR: {source}: image not found: {full}")
    mime = mimetypes.guess_type(full)[0]
    if not mime or not mime.startswith("image/"):
        raise SystemExit(f"ERROR: {source}: not an image: {full}")
    with open(full, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    alt = html.escape(alt)
    return (f'<figure><a href="{html.escape(path)}" title="Open the full-resolution image">'
            f'<img src="data:{mime};base64,{data}" alt="{alt}"></a>'
            f"<figcaption>{alt}</figcaption></figure>")


def list_block(lines, i, source):
    """A bulleted or numbered list; items may wrap onto indented lines or be separated by blank lines."""
    first = LIST_ITEM.match(lines[i].strip())
    ordered = first.group(1) != "-"
    items = []
    while i < len(lines):
        m = LIST_ITEM.match(lines[i].strip())
        if m and (m.group(1) != "-") == ordered and not lines[i][:1].isspace():
            items.append([m.group(2).strip()])
            i += 1
        elif items and lines[i][:1].isspace() and lines[i].strip():
            items[-1].append(lines[i].strip())
            i += 1
        elif not lines[i].strip():
            # A blank line continues the list only when another item of the same kind follows.
            j = i
            while j < len(lines) and not lines[j].strip():
                j += 1
            nxt = LIST_ITEM.match(lines[j].strip()) if j < len(lines) else None
            if nxt and (nxt.group(1) != "-") == ordered and not lines[j][:1].isspace():
                i = j
            else:
                break
        else:
            break
    tag = "ol" if ordered else "ul"
    start = int(first.group(1)[:-1]) if ordered else 1
    attr = f' start="{start}"' if ordered and start != 1 else ""
    body = "".join(f"<li>{inline(' '.join(item), source)}</li>" for item in items)
    return f"<{tag}{attr}>{body}</{tag}>", i


def render(md, source="report", base="."):
    """Return (title, table of contents, body) for the report text."""
    lines = md.split("\n")
    body, toc = [], []
    title = None
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            i += 1
        elif stripped.startswith("```"):
            i += 1
            code = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code.append(lines[i].rstrip())
                i += 1
            if i == len(lines):
                raise SystemExit(f"ERROR: {source}: unclosed code block")
            i += 1
            body.append(f'<pre><code>{html.escape(chr(10).join(code), quote=False)}</code></pre>')
        elif re.match(r"#{1,3} ", line):
            level = len(line.split(" ", 1)[0])
            text = line[level + 1:].strip()
            if level == 1 and title is None:
                title = text
                body.append(f"<h1>{inline(text, source)}</h1>")
            elif level == 2:
                anchor = f"s{len(toc) + 1}"
                toc.append((anchor, text))
                body.append(f'<h2 id="{anchor}">{inline(text, source)}</h2>')
            else:
                body.append(f"<h{level}>{inline(text, source)}</h{level}>")
            i += 1
        elif re.fullmatch(r"-{3,}", stripped):
            body.append("<hr>")
            i += 1
        elif IMAGE.fullmatch(stripped):
            m = IMAGE.fullmatch(stripped)
            body.append(figure(m.group(1), m.group(2), base, source))
            i += 1
        elif LIST_ITEM.match(stripped) and not line[:1].isspace():
            block, i = list_block(lines, i, source)
            body.append(block)
        elif stripped.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(table_row(lines[i]))
                i += 1
            if len(rows) < 2 or not all(re.fullmatch(r":?-+:?", c) for c in rows[1]):
                raise SystemExit(f"ERROR: {source}: malformed table near: {stripped[:80]}")
            body.append(table([rows[0]] + rows[2:], source))
        elif stripped.startswith(">"):
            quote = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                quote.append(lines[i].strip()[1:].strip())
                i += 1
            cls = ""
            if quote and quote[0] == "[!IMPORTANT]":
                cls = ' class="callout"'
                quote = quote[1:]
            body.append(f"<blockquote{cls}>{inline(' '.join(quote), source)}</blockquote>")
        else:
            text = []
            while i < len(lines) and lines[i].strip() and not BLOCK_START.match(lines[i].strip()):
                text.append(lines[i].strip())
                i += 1
            if not text:
                raise SystemExit(f"ERROR: {source}: cannot parse line {i + 1}: {stripped[:80]}")
            body.append(paragraph(" ".join(text), source))
    if title is None:
        raise SystemExit(f"ERROR: {source}: no title heading")
    return title, toc, "\n".join(body)


def convert(path):
    with open(path, encoding="utf-8") as fh:
        md = fh.read()
    name = os.path.basename(path)
    title, toc, body = render(md, name, os.path.dirname(os.path.abspath(path)))
    links = "".join(f'<li><a href="#{a}">{inline(t, name)}</a></li>' for a, t in toc)
    # The contents list goes after the title and the first paragraph (the generated-on line).
    split = body.index("</p>") + len("</p>")
    page = (
        '<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{html.escape(title, quote=False)}</title><style>{CSS}</style></head><body>\n"
        f'<main>\n{body[:split]}\n<nav aria-label="Sections"><ol>{links}</ol></nav>\n{body[split:]}\n</main>\n'
        f"<footer>Rendered from <code>{html.escape(name, quote=False)}</code>.</footer>\n</body></html>\n"
    )
    out = os.path.splitext(path)[0] + ".html"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(page)
    return out


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for p in sys.argv[1:]:
        print(convert(p))
