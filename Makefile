SHELL   := /bin/bash
PACKAGE := $(shell perl -aF: -ne 'print, exit if s/^Package:\s+//' DESCRIPTION)
VERSION := $(shell perl -aF: -ne 'print, exit if s/^Version:\s+//' DESCRIPTION)
BUILD   := $(PACKAGE)_$(VERSION).tar.gz

.PHONY: doc build install test vignette $(BUILD)

check: $(BUILD)
	export _R_CHECK_FORCE_SUGGESTS_=false && R CMD check --as-cran $<

check-no-vignette: $(BUILD)
	export _R_CHECK_FORCE_SUGGESTS_=false && R CMD check --as-cran --no-build-vignettes --ignore-vignettes --no-manual $<


check-cran: $(BUILD)
	R --interactive --no-save --args $< <<<'rhub::check_for_cran(commandArgs(T)[1])'
	Rscript -e 'rhub::check("$(BUILD)", platform = c("solaris-x86-patched"))'

compile:
	find src/ -type f -exec chmod 644 {} \;
	Rscript -e "library(Rcpp); compileAttributes('.');"
	Rscript -e "devtools::load_all(); roxygen2::roxygenise('.');"
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;

build:
	# autoconf
	# chmod 755 cleanup
	# chmod 755 configure
	find src/ -type f -exec chmod 644 {} \;
	chmod 644 ChangeLog DESCRIPTION Makefile NAMESPACE README.md
	# ./configure
	# ./cleanup
	Rscript -e "library(Rcpp); compileAttributes('.');"
	Rscript -e "devtools::load_all(); roxygen2::roxygenise('.');"
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;
	R CMD build .

install:
	# autoconf
	# chmod 755 cleanup
	# chmod 755 configure
	find src/ -type f -exec chmod 644 {} \;
	chmod 644 ChangeLog DESCRIPTION Makefile NAMESPACE README.md
	# ./configure
	# ./cleanup
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;
	Rscript -e "library(Rcpp); compileAttributes('.');"
	Rscript -e "devtools::load_all(); roxygen2::roxygenise('.');"
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;
	R CMD build . # --no-build-vignettes
	R CMD INSTALL $(BUILD)
	
install-fast:
	find src/ -type f -exec chmod 644 {} \;
	chmod 644 ChangeLog DESCRIPTION Makefile NAMESPACE README.md
	find . -iname "*.a" -exec rm {} \;
	find . -iname "*.o" -exec rm {} \;
	find . -iname "*.so" -exec rm {} \;
	R CMD INSTALL .

vignette:
	Rscript -e "rmarkdown::render(input='vignettes/vignette.rmd', output_format='html_vignette')"
	GITHUB_README=Yes Rscript -e "rmarkdown::render(input='vignettes/vignette.rmd', output_file='../README.md', output_format=rmarkdown::github_document(html_preview=FALSE))"; unset GITHUB_README

test:
	IS_LOCAL=Yes Rscript tests/test_pairwise.R && unset IS_LOCAL
	IS_LOCAL=Yes Rscript tests/test_RadixTree.R && unset IS_LOCAL
	IS_LOCAL=Yes Rscript tests/test_RadixForest.R && unset IS_LOCAL

local-bench:
	Rscript inst/extra_tests/benchmark.R


R_INCLUDE=$(shell R CMD config --cppflags)
Rcpp_INCLUDE=$(shell Rscript -e 'cat(system.file("include", package = "Rcpp"))')
BH_INCLUDE=$(shell Rscript -e 'cat(system.file("include", package = "BH"))')
RcppParallel_INCLUDE=$(shell Rscript -e 'cat(system.file("include", package = "RcppParallel"))')

clang-tidy:
	clang-tidy src/CharCounter.cpp -header-filter=inst/include/.* -checks=-*,clang-analyzer-*,clang-analyzer-cplusplus* -extra-arg=-std=c++11 -- $(R_INCLUDE) -Iinst/include -I$(Rcpp_INCLUDE) -I$(BH_INCLUDE) -I$(RcppParallel_INCLUDE)
	clang-tidy src/RadixForest.cpp -header-filter=inst/include/.* -checks=-*,clang-analyzer-*,clang-analyzer-cplusplus* -extra-arg=-std=c++11 -- $(R_INCLUDE) -Iinst/include -I$(Rcpp_INCLUDE) -I$(BH_INCLUDE) -I$(RcppParallel_INCLUDE)
	clang-tidy src/RadixTree.cpp -header-filter=inst/include/.* -checks=-*,clang-analyzer-*,clang-analyzer-cplusplus* -extra-arg=-std=c++11 -- $(R_INCLUDE) -Iinst/include -I$(Rcpp_INCLUDE) -I$(BH_INCLUDE) -I$(RcppParallel_INCLUDE)
	clang-tidy src/pairwise.cpp -header-filter=inst/include/.* -checks=-*,clang-analyzer-*,clang-analyzer-cplusplus* -extra-arg=-std=c++11 -- $(R_INCLUDE) -Iinst/include -I$(Rcpp_INCLUDE) -I$(BH_INCLUDE) -I$(RcppParallel_INCLUDE)

