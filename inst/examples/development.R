# Tests during the development

# packages -----

  library(tidyverse)
  library(stringi)
  library(trafo)

  library(microViz) ## gene expression data and plots
  library(fastTest) ## statistical hypothesis tests

  library(hgnc)

  library(biggrExtra)

# Testing of internal C++ functions --------

  biggrExtra:::bca(1:10)
  biggrExtra:::perci(1:10)

  biggrExtra:::Median(1:10)
  biggrExtra:::vecDraws(1:10, size = 5, n_iter = 10)

# Reaction annotation databases used in the test ---------

  recon2_db <- extract_genes(Recon2D)
  recon2_2_db <- extract_genes(Recon2_2D)

# Prostate cancer differential gene expression estimates in the tests -------

  ## the DGE estimates and their errors will be used by tools
  ## computing reaction activity with differential gene expression data

  tst_x <- set_names(tcga_data$estimate,
                     tcga_data$entrez_id)

  tst_err <- set_names(tcga_data$se,
                       tcga_data$entrez_id)

# Breast cancer gene expression data used in the tests -------

  ## the breat cancer expression data will be used in testing
  ## of tools which compute reaction activity scores

  brca_data <- brca %>%
    column_to_rownames("sample_id") %>%
    select(-patient_id, -histology, -timepoint, -metastasis, -er_status)

  ## annotation of gene symbols with Entrez IDs

  hgnc_annotation <- import_hgnc_dataset()

  brca_annotation <-
    tibble(gene_symbol = names(brca_data),
           entrez_id = exchange(names(brca_data),
                                hgnc_annotation,
                                key = "symbol",
                                value = "entrez_id")) %>%
    filter(complete.cases(.))

  ## selection of genes, re-naming with EntrezIDs,
  ## log2 transformation of the expression data

  brca_data <- set_names(brca_data[, brca_annotation$gene_symbol],
                         brca_annotation$entrez_id)

  brca_data <-
    log2(as.matrix(brca_data) + 1) %>%
    as.data.frame

  ## Z-scores used for caclulation of the reaction activity scores

  brca_z_scores <- zScores(brca_data)

# Testing of the distribution drawing function -------

  draw_matrix <- biggrExtra:::draw_norm(x = tst_x,
                                        err = tst_err,
                                        scale = "log2",
                                        n = 100)

  draw_matrix[1:10, 1:10]
  dim(draw_matrix)

# Testing of evaluation tools: differential gene expression ----------

  ## error derived from normal distribution

  tst_bare_estimates <-
    get_regulation(x = tst_x,
                   scale = "log2",
                   database = recon2_db)

  get_regulation(x = tst_x,
                 err = tst_err,
                 scale = "log2",
                 database = recon2_db,
                 err_method = "norm")

  tst_norm_estimates <-
    get_regulation(x = tst_x,
                   err = tst_err,
                   scale = "log2",
                   database = recon2_2_db,
                   err_method = "norm")

  ## errors derived from MC simulations

  get_regulation(x = tst_x,
                 err = tst_err,
                 scale = "log2",
                 database = recon2_db,
                 err_method = "mc")

  tst_mc_estimates <-
    get_regulation(x = tst_x,
                   err = tst_err,
                   scale = "log2",
                   database = recon2_2_db,
                   err_method = "mc",
                   return_mc = TRUE)

# Testing of evaluation tools: gene expression data --------

  brca_activity_scores <-
    get_activity(x = brca_z_scores,
                 database = recon2_2_db,
                 as_data_frame = TRUE,
                 .parallel = FALSE)

  reaction_variables <- names(brca_activity_scores)[-1]

  ## appending the activity scores with clinical meta-data

  brca_activity_scores <-
    inner_join(brca[, c("sample_id",
                        "patient_id",
                        "timepoint",
                        "metastasis",
                        "histology",
                        "er_status")],
               brca_activity_scores,
               by = "sample_id")

  ## comparison of the activity scores between ER-negative and -positive cancers

  brca_activity_er_test <-
    f_wilcox_test(brca_activity_scores[, reaction_variables],
                  f = brca_activity_scores[["er_status"]],
                  exact = FALSE,
                  adj_method = "BH",
                  as_data_frame = TRUE) %>%
    as_tibble

# Diagnostic plots ----------

  tst_norm_estimates %>%
    plot_errors(fun = log2,
                line_color = "orangered3",
                fill = "gray80")

  tst_mc_estimates %>%
    plot_errors(subsystems = "Fatty acid oxidation",
                fun = sqrt,
                line_color = "orangered3",
                line_alpha = 0.5,
                fill = "gray80")

  tst_norm_estimates %>% plot_mc

  tst_mc_estimates %>%
    plot_mc(subsystems = "Oxidative phosphorylation",
            fun = sqrt)

  tst_mc_estimates %>%
    plot_mc(subsystems = "Oxidative phosphorylation",
            fun = sqrt,
            type = "violin",
            point_alpha = 0.15,
            plot_title = "Oxidative phosphorylation")

# END --------
