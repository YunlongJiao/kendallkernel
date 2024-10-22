# kendallkernel

This repo contains `R` code for reproducing results and more experiments from

> Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." _IEEE_ Transactions on Pattern Analysis and Machine Intelligence(TPAMI), 2017. [DOI: 10.1109/TPAMI.2017.2719680](https://doi.org/10.1109/TPAMI.2017.2719680)

# Cluster full rankings (APA)

```r
setwd("[somepath]/kendallkernel_demo/APA/")
install.packages(c("rmarkdown", "combinat", "caret", "mvtnorm", "reshape2", "ggplot2", "gplots", "flexclust", "Rankcluster", "devtools"))
devtools::install_github("YunlongJiao/kernrank")
rmarkdown::render("analysis.Rmd")
```

# Cluster multivariate partial rankings (eurovision)

```r
setwd("[somepath]/kendallkernel_demo/eurovision/")
install.packages(c("rmarkdown", "combinat", "reshape2", "ggplot2", "gplots", "kernlab", "Rankcluster", "cluster", "devtools"))
devtools::install_github("YunlongJiao/kernrank")
rmarkdown::render("analysis.Rmd")
```

# Classify gene expression data (geneexpr)

```r
setwd("[somepath]/kendallkernel_demo/geneexpr/")
install.packages(c("rmarkdown", "kernlab", "pcaPP", "caret"))
system("R CMD SHLIB src/utiles.c", intern = TRUE)
system("R CMD SHLIB src/tsp.c", intern = TRUE)
rmarkdown::render("analysis.Rmd")
```

# Notes

- Some implementations depend on `parallel::mclapply` which may be inconvenient for MS Windows users. See [here](http://www.r-bloggers.com/implementing-mclapply-on-windows-a-primer-on-embarrassingly-parallel-computation-on-multicore-systems-with-r/) for an easy solution, or manually revert all `mclapply` to `lapply`.
