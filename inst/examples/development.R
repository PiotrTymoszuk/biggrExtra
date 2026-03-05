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

# Expression estimates and databases used in the tests -------

  tst_x <- set_names(tcga_data$estimate,
                     tcga_data$entrez_id)

  tst_err <- set_names(tcga_data$se,
                       tcga_data$entrez_id)

  recon2_db <- extract_genes(Recon2D)
  recon2_2_db <- extract_genes(Recon2_2D)

# Testing of the distribution drawing function -------

  draw_matrix <- biggrExtra:::draw_norm(x = tst_x,
                                        err = tst_err,
                                        scale = "log2",
                                        n = 100)

  draw_matrix[1:10, 1:10]
  dim(draw_matrix)

# Testing of internal evaluation tools: differential gene expression ----------

  ## error derived from normal distribution

  biggrExtra:::get_regulation(x = tst_x,
                              err = tst_err,
                              scale = "log2",
                              database = recon2_db,
                              err_method = "norm")

  biggrExtra:::get_regulation(x = tst_x,
                              err = tst_err,
                              scale = "log2",
                              database = recon2_2_db,
                              err_method = "norm")

  ## errors derived from MC simulations

  biggrExtra:::get_regulation(x = tst_x,
                              err = tst_err,
                              scale = "log2",
                              database = recon2_db,
                              err_method = "mc")

  biggrExtra:::get_regulation(x = tst_x,
                              err = tst_err,
                              scale = "log2",
                              database = recon2_2_db,
                              err_method = "mc",
                              return_mc = TRUE)

# Testing of internal evaluation tools: gene expression data --------






