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
	echo "install.packages(c('knitr','rmarkdown'),repos='http://cran.us.r-project.org',lib=getwd());tinytex::install_tinytex();" | Rscript -
	echo "install.packages('https://github.com/rstudio/tinytex/archive/refs/tags/v0.45.tar.gz',repos=NULL)" | Rscript -
	export R_LIBS=. && R CMD build --no-build-vignettes .
	export R_LIBS=. && R CMD check --no-build-vignettes $(PKG)_$(VERSION).tar.gz
	-mkdir build
	cp $(PKG)_$(VERSION).tar.gz build/$(PKG)_$(VERSION).tar.gz
	
