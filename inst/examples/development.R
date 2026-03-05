# Tests during the development

# packages -----

  library(tidyverse)
  library(stringi)

  library(biggrExtra)

# Testing of internal C++ functions --------

  biggrExtra:::bca(1:10)
  biggrExtra:::perci(1:10)

  biggrExtra:::Median(1:10)
  biggrExtra:::vecDraws(1:10, size = 5, n_iter = 10)

# Testing of the extraction/parsing function -------

  extract_genes(Recon2D)
  extract_genes(Recon2_2D)


