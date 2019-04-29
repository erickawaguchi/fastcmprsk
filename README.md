Introduction
============

fastcmprsk is an R package for performing Fine-Gray regression via a forward-backward scan algorithm.

Features
========
 - Scalable Fine-Gray estimation procedure for large-scale competing risks data.
 - Currently supports unpenalized and penalized (LASSO, ridge, SCAD, MCP) regression.
 - Can perform CIF estimation with interval/band estimation via bootstrap.

Examples
========

```r
library(fastcmprsk)
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2, 200, replace = TRUE)
cov <- matrix(runif(1000), nrow = 200)
dimnames(cov)[[2]] <- c('x1','x2','x3','x4','x5')
fit1 <- fastCrr(ftime, fstatus, cov, variance = FALSE)
summary(fit1)
```
 
Implementation
============
fastcmprsk in an R package with most functionality implemented in C. The package uses cyclic coordinate descent to optimize the likelihood function.


Installation
============
To install the latest development version, install from GitHub. 

```r
install.packages("devtools")
devtools::install_github(“erickawaguchi/fastcmprsk”)
```

System Requirements
===================
Requires R (version 3.5.0 or higher). 

 
User Documentation
==================
* Package manual: Currently unavailable. 

License
=======
fastcmprsk is licensed under GPL-3.  

Development
===========
fastcmprsk is being developed in R Studio.

