R package, aSPU
=================

Il-Youp Kwak <ikwak@umn.edu>

R/aSPU is an R package for Genetic association testing methods such as aSPU, aSPUw, aSPUpath, aSPUs, aSPUsPath, GEEaSPU, MTaSPUs, GATES, GATE-Simes, HYST etc.

### Summary table 

| Function Name | Data Type  | Description                                    |
|---------------|------------|------------------------------------------------|
| aSPU, aSPUw   | Individual | Single trait; gene-based                       |
| aSPUs         | Summary    | Single trait; gene-based                       |
| aSPUpath      | Individual | Single trait; pathway-based                    |
| aSPUsPath     | Summary    | Single trait; pathway-based                    |
| GEEaSPU       | Individual | Multiple traits; single SNP based              |
| MTaSPUs       | Summary    | Multiple traits; single SNP based              |
|---------------|------------|------------------------------------------------|

*Data type indicate the structure of data set. "Individual" for individual level data. "Summary" for summary statistics data (such as Z scores or p-values of each SNP) 


### installation
From `CRAN` :
```S
install.packages("aSPU")
```

Or, with `devtools`:
```S
library(devtools)
install_github("ikwak2/aSPU")
```
