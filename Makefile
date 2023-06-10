VERSION := 0.2.0
PKG     := $(shell basename `pwd`)
build: man/snha.Rd man/mgraph.Rd
	R CMD build .

check: build
	R CMD check $(PKG)_$(VERSION).tar.gz

man/%.Rd: R/%.R
	Rscript bin/rman.R $<
