# Makefile for Sphinx documentation

SPHINXBUILD   = sphinx-build

.PHONY: clean html 

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  clean     to delete all generated files in _build"

clean:
	-rm -rf _build/*

# sphinx options:
#   -E forces it to rebuild the docs from scratch every time.
#   -d keeps the intermediate doctree files out the html directory.
html:
	$(SPHINXBUILD) -b html -E -d _build/doctrees . _build/doc
	cd _build && zip -rq snappy_doc.zip doc
	@echo
	@echo "HTML build done. Result was _build/doc and _build/snappy_doc.zip."
