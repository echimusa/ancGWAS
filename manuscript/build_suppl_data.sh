#!/bin/bash
/usr/bin/latex supplementary_data
/usr/bin/bibtex supplementary_data
/usr/bin/latex supplementary_data
/usr/bin/latex supplementary_data
/usr/bin/pdflatex supplementary_data

rm supplementary_data.aux
rm supplementary_data.bbl
rm supplementary_data.blg
rm supplementary_data.dvi
rm supplementary_data.log

mv supplementary_data.pdf output_plos_format

