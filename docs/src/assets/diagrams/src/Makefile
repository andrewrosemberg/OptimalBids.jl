# Builds the market diagram from market.tex: PDF with pdflatex, SVG with
# pdftocairo, then a white rectangle under the drawing so it reads on dark pages too.
#
#   make -C docs/src/assets/diagrams/src
#
# Needs a TeX Live with pdflatex and poppler's pdftocairo.
OUT := ..
TMP := build
NAME := optimalbids--market

all: $(OUT)/$(NAME).svg $(OUT)/$(NAME)-compact.svg

$(TMP)/%.pdf: market.tex
	mkdir -p $(TMP)
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=$(TMP) -jobname=$* market.tex >/dev/null

$(OUT)/%.svg: $(TMP)/%.pdf
	pdftocairo -svg $< $(TMP)/$*.raw.svg
	perl -0pe 's/(<svg [^>]*>)/$$1\n<rect width="100%" height="100%" fill="#ffffff"\/>/' $(TMP)/$*.raw.svg > $@

clean:
	rm -rf $(TMP)

.PHONY: all clean

.SECONDARY:
