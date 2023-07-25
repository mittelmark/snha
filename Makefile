VERSION := 0.2.0
PKG     := $(shell basename `pwd`)
build: man/snha.Rd man/mgraph.Rd
	R CMD build --no-build-vignettes .

check: build
	R CMD check --no-build-vignettes $(PKG)_$(VERSION).tar.gz

man/%.Rd: R/%.R
	Rscript bin/rman.R $<

install: check
	R CMD INSTALL $(PKG)_$(VERSION).tar.gz

install-ubuntu: 
	echo "install.packages(c('knitr','rmarkdown');" | Rscript -e -
