#!/usr/bin/tcsh

echo "==================================="
echo " Generating XML-REFERENCE document"
echo "==================================="

set DOCUMENT=xmlReference

pdflatex $DOCUMENT
bibtex   $DOCUMENT
pdflatex $DOCUMENT
pdflatex $DOCUMENT
pdflatex $DOCUMENT
thumbpdf $DOCUMENT
pdflatex $DOCUMENT

echo "==================================="
echo " See xmlReferemce.pdf for details"
echo "==================================="
