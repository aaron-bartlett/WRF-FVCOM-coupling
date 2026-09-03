TEX = pdflatex -shell-escape -interaction=nonstopmode -file-line-error
BIB = bibtex

.PHONY: all view clean

FILE=nuwrf_userguide

all: $(FILE).pdf

view:
	evince $(FILE).pdf

$(FILE).pdf: $(FILE).tex
	$(TEX) $(FILE)
	$(BIB) $(FILE)
	$(TEX) $(FILE)
	$(TEX) $(FILE)
	$(TEX) $(FILE)
	#$(BIB) $(FILE)
	#$(TEX) $(FILE)
	#$(TEX) $(FILE)

clean:
	rm -f $(FILE).pdf
	rm -f *.aux
	rm -f *.bbl
	rm -f *.blg
	rm -f *.dvi
	rm -f *.log
	rm -f *.nav
	rm -f *.out
	rm -f *.snm
	rm -f *.toc
	rm -f *.vrb
