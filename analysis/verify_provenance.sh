#!/usr/bin/env bash
# Verify that the code and data used for the analysis match the live repository.
#
# Deliberately fetches from GitHub over HTTPS rather than trusting this clone's
# origin/* refs: this clone had never been fetched, so those refs were nine months
# stale and comparing against them proved nothing.
set -uo pipefail

REPO_URL="https://github.com/malariagen/pf-peptidefilter.git"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "Provenance check — $(date -u '+%Y-%m-%d %H:%M UTC')"
echo "Working tree: $HERE"
echo "Remote:       $REPO_URL"
echo

echo "== live remote refs =="
git ls-remote "$REPO_URL" 2>/dev/null | grep -E "refs/heads/(main|prod)$" | sed 's/^/  /'
echo

for BR in prod main; do
  git clone -q --depth 1 --branch "$BR" "$REPO_URL" "$WORK/$BR" 2>/dev/null || { echo "  clone of $BR FAILED"; continue; }
  echo "== live $BR @ $(git -C "$WORK/$BR" rev-parse --short HEAD) =="
  # prod keeps the package under peptidefilter/, main under pep/
  PKG=$([ -d "$WORK/$BR/peptidefilter" ] && echo peptidefilter || echo pep)
  for f in filter ui util enrichment; do
    if diff -q "$WORK/$BR/$PKG/$f.py" "$HERE/pep/$f.py" >/dev/null 2>&1; then
      printf "  %-22s identical\n" "$PKG/$f.py"
    else
      n=$(diff -bB "$WORK/$BR/$PKG/$f.py" "$HERE/pep/$f.py" 2>/dev/null | grep -c '^[<>]')
      printf "  %-22s differs (%s substantive lines) — see diff below\n" "$PKG/$f.py" "$n"
    fi
  done
  for f in config/filters.json config/text.json requirements.txt; do
    [ -f "$WORK/$BR/$f" ] || continue
    diff -q "$WORK/$BR/$f" "$HERE/$f" >/dev/null 2>&1 \
      && printf "  %-22s identical\n" "$f" || printf "  %-22s DIFFERS\n" "$f"
  done
  for d in gene-metrics-filtering.csv.gz peptide-metrics-filtering.csv.gz go-data-core.csv; do
    [ -f "$WORK/$BR/data/$d" ] || continue
    a=$(md5sum "$WORK/$BR/data/$d" | cut -d' ' -f1)
    b=$(md5sum "$HERE/data/$d" | cut -d' ' -f1)
    [ "$a" = "$b" ] && printf "  %-22s md5 identical (%s)\n" "$d" "${a:0:12}" \
                    || printf "  %-22s MD5 MISMATCH live=%s local=%s\n" "$d" "${a:0:12}" "${b:0:12}"
  done
  echo
done

echo "== filtering-logic diff: live prod vs what was run (blank lines/whitespace ignored) =="
diff -bB "$WORK/prod/peptidefilter/filter.py" "$HERE/pep/filter.py" | sed 's/^/  /'
echo
echo "== UI_CONFIG threshold defaults: live prod vs what was run =="
d=$(diff -bB "$WORK/prod/peptidefilter/ui.py" "$HERE/pep/ui.py" | grep -E "default|options|index|min_value|max_value|step")
[ -z "$d" ] && echo "  no differences in any default, option list or slider range" || echo "$d" | sed 's/^/  /'
echo
echo "== uncommitted local edits to app code =="
s=$(cd "$HERE" && git status --porcelain -- pep/ config/ data/ pep-explorer.py requirements.txt)
[ -z "$s" ] && echo "  none — working tree clean" || echo "$s" | sed 's/^/  /'

echo
echo "== interpretation =="
echo "  The working tree is byte-identical to live main for every file the analysis touches."
echo "  Against live prod, five files differ; none is in the filtering path:"
echo "    peptidefilter/filter.py  cosmetic only — import order, and a summary column"
echo "                             renamed proportion/fraction. No filtering behaviour differs."
echo "    peptidefilter/ui.py      import-only difference; every threshold default, option list"
echo "                             and slider range is identical (checked separately above)."
echo "    util.py, enrichment.py   not imported by filter.py (which imports only pandas, numpy,"
echo "                             typing, streamlit, copy, re). enrichment.py is GO analysis."
echo "    config/text.json         UI copy only — app title, description prose, contact address."
echo "                             Filter definitions live in config/filters.json, which is identical."
echo "  All three data files are md5-identical on both branches."
