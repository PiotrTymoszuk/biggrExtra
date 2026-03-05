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

# Testing of the distribution drawing function -------

  tst_x <- set_names(tcga_data$estimate,
                     tcga_data$entrez_id)

  tst_err <- set_names(tcga_data$se,
                       tcga_data$entrez_id)

  draw_matrix <- biggrExtra:::draw_norm(x = tst_x,
                                        err = tst_err,
                                        scale = "log2",
                                        n = 100)

  draw_matrix[1:10, 1:10]
  dim(draw_matrix)




