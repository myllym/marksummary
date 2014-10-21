# marksummary

An R library for calculating summary functions for marked point patterns.
Currently it implements the (mark-weighted) K_f functions for different mark
test functions f. It also provides the ability to perform simulations for
the random labelling hypothesis.

The main function in `marksummary` is summ_func_random_labelling.
The recommended use of `marksummary` is through the R library `spptest`
(functions: random_labelling).

## Installation

Install the library to R with the following two R commands:

```R
library(devtools)
install_github('myllym/marksummary')
```

If you do not have the library ´devtools´ installed, install it first by

```R
install.packages("devtools")
```

