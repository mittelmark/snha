VERSION := $(shell grep -E '^Version:' DESCRIPTION | sed 's/Version: //')
PKG     := $(shell basename `pwd`)
build: man/snha.Rd man/mgraph.Rd
	R CMD build .

check: build
	R CMD check $(PKG)_$(VERSION).tar.gz
	echo "@page { size: 12in 9in; }  body { font-size: 80% ; }" > small.css
	weasyprint snha.Rcheck/snha/doc/tutorial.html -s small.css snha-tutorial.pdf
	### --no-build-vignette


man/%.Rd: R/%.R
	Rscript bin/rman.R $<

install: check
	R CMD INSTALL $(PKG)_$(VERSION).tar.gz

install-ubuntu: 
	echo "Starting build!"
	sudo apt-get install texlive-base
	sudo apt-get install texlive-latex-recommended
	sudo apt-get install texlive-latex-extra
	sudo apt-get install texlive-pictures
	sudo apt-get install pandoc
	mkdir rlibs
	echo "install.packages(c('knitr','rmarkdown'),repos='http://cran.us.r-project.org',lib=file.path(getwd(),'rlibs'));" | Rscript -
	export R_LIBS=./rlibs && R CMD build --no-build-vignettes .
	export R_LIBS=./rlibs && R CMD check --no-build-vignettes $(PKG)_$(VERSION).tar.gz
	#mkdir build
	#cp $(PKG)_$(VERSION).tar.gz build/$(PKG)_$(VERSION).tar.gz
	
