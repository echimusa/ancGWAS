#!/bin/bash
/usr/bin/latex september_ancGWAS
/usr/bin/bibtex september_ancGWAS
/usr/bin/latex september_ancGWAS
/usr/bin/latex september_ancGWAS
/usr/bin/pdflatex september_ancGWAS

rm september_ancGWAS.aux
rm september_ancGWAS.bbl
rm september_ancGWAS.blg
rm september_ancGWAS.dvi
rm september_ancGWAS.log
#rm -r output_plos_format
mkdir output_plos_format

mv september_ancGWAS.pdf output_plos_format
