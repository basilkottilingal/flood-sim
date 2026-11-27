# hmesh.pdf: libhmesh as it is, with the plans to come
DOC = hmesh

all: $(DOC).pdf

$(DOC).pdf: $(DOC).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(DOC).tex > /dev/null
	pdflatex -interaction=nonstopmode -halt-on-error $(DOC).tex > /dev/null

clean:
	rm -f $(DOC).pdf $(DOC).aux $(DOC).log $(DOC).out $(DOC).toc

.PHONY: all clean
