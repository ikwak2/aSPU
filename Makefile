doc:
	R -e 'library(devtools);document(roclets=c("namespace", "rd"))'

build:
	R -e 'library(devtools);build(“../aSPU”)’

vig:
	cd vignettes; R -e 'library(knitr);knit2html("aspu.Rmd")'

update:
	R -e 'library(devtools);install_github("ikwak2/aSPU”)