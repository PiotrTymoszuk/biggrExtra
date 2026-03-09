# Calculation of reaction activity scores with log2 gene expression data
# in a breast cancer data set (delivered with `microViz` package)

# packages --------

  library(tidyverse)

  library(microViz) ## data and plots
  library(fastTest) ## statistical hypothesis tests
  library(hgnc) ## gene symbol - Entrez ID conversion

  library(biggrExtra)

# analysis data ------

  brca_data <- brca %>%
    column_to_rownames("sample_id")

  ## annotation of gene symbols with Entrez IDs

  hgnc_annotation <- import_hgnc_dataset() %>%
    select(symbol, entrez_id) %>%
    filter(complete.cases(.))

  hgnc_entrez_vector <- set_names(hgnc_annotation$entrez_id,
                                  hgnc_annotation$symbol)

  ## selection of genes, re-naming with EntrezIDs,
  ## log2 transformation of the expression data

  brca_data <- brca_data %>%
    select(any_of(hgnc_annotation$symbol))

  names(brca_data) <- hgnc_entrez_vector[names(brca_data)]

  brca_data <-
    log2(as.matrix(brca_data) + 1) %>%
    as.data.frame

  ## Z-scores used for calculation of the reaction activity scores

  brca_z_scores <- zScores(brca_data)

# Calculation of the activity scores ---------

  ## reaction annotation database created from
  ## the `Recon2_2D` data set from `biggrExtra` package

  recon2_2_db <- as_reactDB(Recon2_2D)

  brca_activity_scores <-
    get_activity(x = brca_z_scores,
                 database = recon2_2_db,
                 as_data_frame = TRUE)

  reaction_variables <- names(brca_activity_scores)[-1] ### to be used later in tests

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

# comparison of activity scores between major histologies ----------

  ## differences in activity scores between major histologies:
  ## Mann-Whitney test

  brca_test_data <- brca_activity_scores %>%
    filter(histology %in% c("IDC", "ILC"))

  brca_histology_test <-
    f_wilcox_test(brca_test_data[, reaction_variables],
                  f = brca_test_data$histology,
                  as_data_frame = TRUE,
                  safely = TRUE,
                  adj_method = "BH") %>%
    as_tibble

  ## significant differences in activity scores:
  ## pFDR < 0.05 and large (abs(r) > 0.5)

  brca_histology_significant <- brca_histology_test %>%
    filter(p_adjusted < 0.05,
           abs(biserial_r) >= 0.5) %>%
    .$variable

  ## visualization of the activity score in histological
  ## subtypes in a heat map

  brca_histo_heat_map <- brca_test_data %>%
    heat_map(variables = brca_histology_significant,
             split_fct = "histology",
             normalize = FALSE,
             midpoint = 0,
             hide_x_axis_text = TRUE,
             limits = c(-5, 5),
             oob = scales::squish,
             x_lab = "cancer sample",
             y_lab = "metabolic reaction",
             fill_lab = "activity\nscore",
             plot_title = "Predicted differences in metabolism, IDC and ILC") +
    theme(strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          axis.text.y = element_blank())

# END -------
