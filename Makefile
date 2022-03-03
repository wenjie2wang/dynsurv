objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell grep "^Version" DESCRIPTION | awk '{print $$NF}')
pkg := $(shell grep "^Package" DESCRIPTION | awk '{print $$NF}')
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log


.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

.PHONY: pkgdown
pkgdown:
	Rscript -e "library(methods); pkgdown::build_site();"

$(tar): $(objects)
	@$(MAKE) updateTimestamp
	@Rscript -e "library(methods); devtools::document();";
	R CMD build .

$(checkLog): $(tar) $(tinytest)
	R CMD check $(tar)

.PHONY: check-as-cran
check-as-cran: $(tar)
	R CMD check --as-cran $(tar)

## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateTimestamp
updateTimestamp:
	@bash misc/update_timestamp.sh

## make tags
.PHONY: TAGS
TAGS:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"
	gtags

.PHONY: clean
clean:
	@$(RM) -rf *~ */*~ src/*.o src/*.so *.Rhistroy *.tar.gz *.Rcheck/ .\#*
