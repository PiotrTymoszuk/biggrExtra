# Tests during the development

# packages -----

  library(tidyverse)
  library(rlang)
  library(stringi)

  library(biggrExtra)

# Testing of internal C++ functions --------

  biggrExtra:::bca(1:10)
  biggrExtra:::perci(1:10)

  biggrExtra:::Median(1:10)
  biggrExtra:::vecDraws(1:10, size = 5, n_iter = 10)
