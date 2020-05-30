#!/bin/sh

base="${1%.*}"

rmd_preview(){
  if grep -q "xaringan\|bookdown" "$1"; then 
    use chrome "$base.html"
  elif grep -q "pdf_document\|pandoc-crossref" "$1"; then
    use skim "$base.pdf"
  fi
}

case "$1" in
    *.html) use chrome "$1";;
    *.Rmd) rmd_preview "$1";;
    *.pdf) use skim "$1";;
esac
