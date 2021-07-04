# ACAT (Aggregated Cauchy Association Test)
This package implemented ACAT, a generic method for combing p-values, and ACATV, a variant set test for rare variant analysis using ACAT.
## Installation
```
library(devtools)
devtools::install_github("yaowuliu/ACAT")
```
## Usage
Please see the [ACAT user manual](https://github.com/yaowuliu/ACAT/blob/master/doc/ACAT_manual.pdf).

## FAQ
**Q: How to deal with p-values that are exactly 1 in ACAT？If some of the p-values combined by ACAT are exactly 1, ACAT will give a combined p-value of 1.**

**A:** When a test statistic follows a continuous ditribution, the corresponding p-value should follow a uniform distribution between 0 and 1 under the null. Hence, for continuous distributions, one should never get a p-value that is exactly 1. In practice, we may have p-values being 1 becasue the calibrated p-value is an approximation of the "true" p-value. For example, if simulation-based methods (e.g., permutation) is used, one could have p-values equal to 1 since the number of simulations is always finite. 

Therefore, if there are p-values equal to 1, we need to first find out the reason and then try to "correct" the calibrated p-values. For example, if it is due to permutation, we can replace 1 by 1-1/N, where N is the number of permutations. Another simple way is to replace 1 by 1-1/d, where d is the number of p-values combined by ACAT.   
