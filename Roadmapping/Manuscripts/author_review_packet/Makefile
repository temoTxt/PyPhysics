PAPER = paper

.PHONY: all clean

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex references.bib sections/*.tex equations/*.tex figures/cross_comparison_table.tex
	pdflatex -interaction=nonstopmode $(PAPER) >/dev/null
	bibtex $(PAPER) || true
	pdflatex -interaction=nonstopmode $(PAPER) >/dev/null
	pdflatex -interaction=nonstopmode $(PAPER) >/dev/null

references.bib:
	uv run python ../../History/Bibliography/build_bibtex.py --output references.bib || \
	  echo "% build_bibtex.py unavailable; placeholder references.bib" > references.bib

clean:
	rm -f *.aux *.log *.bbl *.blg *.out *.toc *.fls *.fdb_latexmk
