# partition-methods.pdf: the figures are computed by figures.py
DOC = partition-methods

all: $(DOC).pdf

fig-counts.tex: figures.py
	python3 figures.py

$(DOC).pdf: $(DOC).tex definitions.tex algorithms.tex fig-counts.tex
	pdflatex -interaction=nonstopmode -halt-on-error $(DOC).tex > /dev/null
	pdflatex -interaction=nonstopmode -halt-on-error $(DOC).tex > /dev/null

clean:
	rm -f $(DOC).pdf $(DOC).aux $(DOC).log $(DOC).out $(DOC).toc fig-*.tex

.PHONY: all clean
