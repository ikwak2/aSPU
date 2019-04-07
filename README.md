## R package, aSPU

[![CRAN](http://www.r-pkg.org/badges/version/aSPU)](https://CRAN.R-project.org/package=aSPU)


Il-Youp Kwak <ilyoup.kwak@gmail.com>

R/aSPU is an R package for Genetic association testing methods such as aSPU, aSPUw, aSPUpath, aSPUs, aSPUsPath, GEEaSPU, MTaSPUs, GATES, GATE-Simes, HYST etc.

### Summary table 

| Function Name         | Data Type  | Description                                    |
|-----------------------|------------|------------------------------------------------|
| aSPU, aSPUw, aSPUr, aSPUd    | Individual | Single trait; gene-based                       |
| aSPUs                 | Summary    | Single trait; gene-based                       |
| aSPUpath              | Individual | Single trait; pathway-based                    |
| aSPUsPath             | Summary    | Single trait; pathway-based                    |
| GEEaSPU               | Individual | Multiple traits; single SNP based              |
| MTaSPUs               | Summary    | Multiple traits; single SNP based              |
| MTaSPUsSet            | Summary    | Multiple traits; gene-based                    |
| MTaSPUsSetPath        | Summary    | Multiple traits; pathway-based                 |


*Data type indicate the structure of data set. "Individual" for individual level data. "Summary" for summary statistics data (such as Z scores or p-values of each SNP) 

* aSPU is the function for original aSPU test (Pan et al. 2014), aSPUw is weghted version of it (Kim et al. 2014), aSPUr is robust version of aSPU test (Wei et al. 2016), aSPUd is aSPU test using asymptotic distribution of SPU statistics (Gong et al. 2016). The original version of aSPUd is two sample mean comparison available in R [highmean](https://CRAN.R-project.org/package=highmean) package.

### Tutorials 

Some [tutorials](https://github.com/ikwak2/aSPU_tutorials) for aSPUs, aSPUsPath and MTaSPUsSet are available. 
 - [Mapping Snp to Gene](https://ikwak2.github.io/tutorials/mappingSnpToGene2.html)
 - [Getting correlation estimate among SNPs from reference panel](https://ikwak2.github.io/tutorials/CorrFromRef.html)
 - [R and Perl codes to perform aSPUs and MTaSPUsSet](https://ikwak2.github.io/tutorials/ForMTgenes.html)
 - [Pathway data manipulation](https://ikwak2.github.io/tutorials/pathwayMani.html)
 - [Speed comparison using R, awk and Perl](https://ikwak2.github.io/tutorials/SpeedComp.html)
 - [Vignette for aSPUs and aSPUsPath](https://ikwak2.github.io/tutorials/aSPUstat.html)
 - [Vignette for MTaSPUsSet](https://ikwak2.github.io/tutorials/MTaSPUsSet.html)

### Citations

For 'aSPU'
```
Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
A powerful and adaptive association test for rare variants,
Genetics, 197(4), 1081-95
```

For 'aSPUw'
```
Junghi Kim, Jeffrey R Wozniak, Bryon A Mueller, Xiaotong Shen and Wei Pan (2014)
Comparison of statistical tests for group differences in brain functional networks,
NeuroImage, 1;101:681-694
```

For 'aSPUr'
```
Peng Wei, Ying Cao, Yiwei Zhang, Zhiyuan Xu, Il-Youp Kwak, Eric Boerwinkle, Wei Pan (2016)
On Robust Association Testing for Quantitative Traits and Rare Variants, 
G3, 6(12) 3941-3950. 
```

For 'aSPUd'
```
Gongjun Xu, Lifeng Lin, Peng Wei and Wei Pan (2016) 
An adaptive two-sample test for high-dimensional means, 
Biometrika (2016) 103 (3): 609-624.
```

For 'aSPUpath'
```
Wei Pan, Il-Youp Kwak and Peng Wei (2015)
A Powerful and Pathway-Based Adaptive Test for Genetic Association With Common or Rare Variants,
The American Journal of Human Genetics 97, 86-98
```

For 'aSPUs' and 'aSPUsPath'
```
Il-Youp Kwak, Wei Pan (2015)
Adaptive Gene- and Pathway-Trait Association Testing with GWAS Summary Statistics,
Bioinformatics, 32(8), 1178-1184
```

For 'GEEaSPU'
```
Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014)
Testing for association with multiple traits in generalized estimation equations, with application to neuroimaging data,
Neuroimage. 96, 309-325
```

For 'MTaSPUs'
```
Junghi Kim, Yun Bai and Wei Pan (2015)
An Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics,
Genetic Epidemiology, 8:651-663
```

For 'MTaSPUsSet' and 'MTaSPUsSetPath'
```
Il-Youp Kwak, Wei Pan (2017)
Gene- and pathway-based association tests for multiple traits with GWAS summary statistics, Bioinformatics. 33(1), 64-71
```


For 'GATES'
```
Miao-Xin Li, Hong-Sheng Gui, Johnny S.H. Kwan and Pak C. Sham (2011)
GATES: A Rapid and Powerful Gene-Based Association Test Using Extended Simes Procedure,
The American Journal of Human Genetics 88, 283-293
```

For 'GATES-Simes'
```
Hongsheng Gui, Miaoxin Li, Pak C Sham and Stacey S Cherny (2011)
Comparisons of seven algorithms for pathway analysis using the WTCCC Crohn's Disease
BMC Research Notes, 4:386
```

For 'HYST'
```
Miao-Xin Li, Johnny S.H. Kwan and Pak C. Sham (2012)
HYST: A Hybrid Set-Based Test for Genome-wide Association Studies, with Application to Protein-Protein Interaction-Based Association Analysis
The American Journal of Human Genetics, 91, 478-488.
```



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


### License

The R/aSPU package is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>
