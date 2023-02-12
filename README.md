# Rpkg

A generic R package template.

Just use the green "Use this template" link on top to create your own R-package in Github.

You need to edit thereafter at least two files where your place the right
package name and change the other relevant information like author, license
etc.

* DESCRIPTION
* tests/test-add.R

If you do not want to create a package vignette just remove the entry in the
description file with the vignette information, so the lines starting with
`Suggests:` and `VignetteBuilder:`. Your package has then zero dependencies from other R packages.


## Build and Install test

To build, check and install this package you need only an installed R. If you
create your repository you have to change the relevant parts entering, package
name, author etc in the description file and the library call in the file
tests/test-add.R. Then you can start a first check if the package is
installable.


The commands are the following:

```
R CMD build .                   # creates a tar.gz file
R CMD check pkgname_0.1.tar.gz  # replace pkgname with your name
```

## Documenting your code

If this is all ok, you can start entering your own code and update your documentation and tests. The template provides as well a simple R script which copies the documentation out of the R-files into the man-directory. If you do not like this approach of combining the documentation with the code in one file you can as well just skip this file and create directly your own Rd-files. The documentation out of the R files is extracted like this:

```
Rscript bin/rman.R R/add.R
```

This approach, embedding the documentation in the R file and then extracting
it into the man-folder is very similar to that of the documentation tool
[roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html), but
here you do not need the roxygen2 package with all the dependencies, just this simple extraction script
`bin/rman.R`. 

## Makefile

There is as well a Makefile file which provides the necessary commands to
build and install the package and which does as well a check if the R file is
newer than the Rd file and only then does the required update. If you use make
you can then simple write:

```
make check VERSION=0.1
```

if your version number is 0.1.

If you do not use make just execute the required commands manually in the terminal.

```
Rscript bin/rman.R R/add.R
R CMD build .
R CMD check pkgname_version.tar.gz
R CMD INSTALL pkgname_version.tar.gz
```

## Files

The following files are included in this package template and will be required for the package installation:

* `DESCRIPTION` - the file containing the essential information about your package
* `NAMESPACE` - the file containing the export statements, usually all functions with lower case methods are exported
* `LICENSE` - license files can be replaced if you prefere other licenses, such as GPL or BSD
* `R/add.R` - example R code file with embedded documentation
* `man/add.Rd` - example R documentation file 

The following files are usually not part of your package:

* `Makefile` - a file for make containing essential commands to build and check the package, should be not part of the installed package, but might be helpful during development
* `bin/rman.R` - simple R script to extract documentation from the R files, can be used for your own package instead of using roxygen2 documentation, if you edit the files in the `man` folder directly you do not need this file
* `README.md` - this Readme your are reading, should be not part of your package, replace it wiht your own Readme file


## Author and Copyright

Author: Detlef Groth, University of Potsdam, Germany

License: MIT License see the file [LICENSE](LICENSE) for details.

## Bug reporting

In case of bugs and suggestions, use the [https://github.com/mittelmark/Rpkg/issues](issues) link on top.
