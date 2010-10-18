#!/bin/bash

rst2html.py src/fileformat.txt > fileformat.html
rst2newlatex.py src/fileformat.txt > fileformat.tex
pdflatex fileformat.tex

/bin/rm -f fileformat.aux fileformat.log fileformat.out fileformat.tex

rst2html.py src/usersguide.txt > usersguide.html
rst2newlatex.py src/usersguide.txt > usersguide.tex
pdflatex usersguide.tex

/bin/rm -f usersguide.aux usersguide.log usersguide.out usersguide.tex
