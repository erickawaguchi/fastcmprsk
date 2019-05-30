Introduction
============

fastcmprsk is an R package for performing Fine-Gray regression via a forward-backward scan algorithm.

Features
========
 - Scalable Fine-Gray estimation procedure for large-scale competing risks data.
 - Currently supports unpenalized and penalized (LASSO, ridge, SCAD, MCP) regression.
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

License
=======
fastcmprsk is licensed under GPL-3.  

Development
===========
fastcmprsk is being developed in R Studio.

