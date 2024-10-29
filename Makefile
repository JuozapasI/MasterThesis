all: thesis.pdf

%.pdf: %.tex
	latexmk -pdf -silent -deps-out=.depend $*

clean:
	latexmk -c
	rm *.bbl
	rm *.bib.bak
