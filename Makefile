VERSION := 0.1.1
PKG     := $(shell basename `pwd`)
build: man/snha.Rd
	R CMD build .

check: build
	R CMD check $(PKG)_$(VERSION).tar.gz

man/%.Rd: R/%.R
	Rscript bin/rman.R $<
