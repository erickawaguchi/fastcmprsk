Introduction
============
fastcmprsk is an R package for performing Fine-Gray regression via a forward-backward scan algorithm.

Official CRAN release (v1.1.1) is available [here](https://cran.r-project.org/web/packages/fastcmprsk/index.html).

NOTE TO USERS: We plan to make monthly/quarterly updates to the package!

What’s New in Version 1.21.1?
========

1. Allows for Fine-Gray regression w.o presence of right censoring.

Features
========
 - Scalable Fine-Gray parameter estimation procedure for large-scale competing risks data.
 - Currently supports unpenalized and penalized (LASSO, ridge, SCAD, MCP, elastic-net) regression.
 - Can perform CIF estimation with interval/band estimation via bootstrap.
 
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
* Please cite [Kawaguchi et al. (2019)](https://arxiv.org/abs/1905.07438).

License
=======
fastcmprsk is licensed under GPL-3.  

Development
===========
fastcmprsk is being developed in R Studio.
