FILE=A-MapReduce-like-Model-for-Madness

pdf: ${FILE}.pdf
dvi: ${FILE}.dvi

${FILE}.pdf: ${FILE}.tex
	pdflatex $<

${FILE}.dvi: ${FILE}.tex
	latex $<

spell: ${FILE}.tex detex.sed
	sed -f detex.sed ${FILE}.tex | spell | sort | uniq | less

xdvi: ${FILE}.dvi
	xdvi $<

okular: ${FILE}.pdf
	okular $<

evince: ${FILE}.pdf
	evince $<

emacs: ${FILE}.tex
	emacs -nw $<

vi: ${FILE}.tex
	vi $<

renoir: ${FILE}.pdf
	lpr -Prenoir-duplex $<

clean:
	rm -f ${FILE}.aux ${FILE}.log ${FILE}.tex~ .${FILE}.swp ${FILE}.dvi \
	  ${FILE}.ps

dist: clean ${FILE}.tex ${FILE}.pdf
	dir=`basename $$PWD`; cd ..; tar zcvf $$dir.tar.gz ./$$dir
	dir=`basename $$PWD`; ls -l ../$$dir.tar.gz
