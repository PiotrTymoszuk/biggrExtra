# Example of calculation of fold-differences in activity of Recon2.2 metabolic
# reactions between collagen subsets of prostate cancer
# (TCGA data set, DOI: 10.1016/j.euo.2024.05.014).

# packages --------

  library(tidyverse)
  library(biggrExtra) ### modeling of metabolism

# analysis data ------

  ## the `tcga_data` data set, delivered with `biggrExtra` package
  ## differential regulation: estimates of log2 fold-changes
  ## with their errors, named with Entrez IDs of the genes

  dge_estimates <- set_names(tcga_data$estimate,
                             tcga_data$entrez_id)

  dge_errors <- set_names(tcga_data$se,
                          tcga_data$entrez_id)

  dge_estimates[1:10]
  dge_errors[1:10]

  ## reaction annotation database created from
  ## the `Recon2_2D` data set from `biggrExtra` package

  recon2_2_db <- as_reactDB(Recon2_2D)

# Calculation of fold changes of reaction activity --------

  ## predictions of fold-changes in reaction activity,
  ## errors derived from Monte Carlo simulations

  reaction_regulation <-
    get_regulation(x = dge_estimates, ### vector of DGE estimates
                   err = dge_errors, ### vector of DGE errors
                   database = recon2_2_db, ### reaction annotation database
                   scale = "log2", ### scale of the DGE estimates
                   err_method = "mc", ### error and p value calculation method
                   return_mc = TRUE, ### reaction activity for each iteration saved
                   n_iter = 2000, ### number of iterations
                   .parallel = TRUE)

  ## calculation of log2 estimates of differential expression
  ## identification of differentially regulated reactions
  ## (pFDR < 0.05, at least 1.1-fold change in activity)

  reaction_regulation <- reaction_regulation %>%
    transform_estimates(fun = log2, prefix = "log2_")

  reaction_regulation <- reaction_regulation %>%
    identify_regulated(p_type = "adjusted",
                       p_cutoff = 0.05,
                       fold_cutoff = 1.1)

  ## visual diagnostic:
  ## plot of log2 errors and log2 estimates of
  ## regulation of reaction activity

  est_error_plot <- reaction_regulation %>%
    plot(plot_type = "errors",
         fun = log2,
         plot_title = "Estimates and errors of reaction regulation",
         x_lab = "log2 fold-regulation estimate",
         y_lab = "log2 error")

# Enrichment analyses for metabolic subsystems ---------

  ## identification of significantly over-represented
  ## metabolic subsystems among activated and inhibited reactions
  ## by random sampling test

  subsystem_enrichment <- reaction_regulation %>%
    suba(type = "random",
         n_iter = 2000,
         .parallel = TRUE)

  ## significantly enriched subsystems:
  ## OR >= 2 and pFDR < 0.05

  significant_subsystems <- subsystem_enrichment %>%
    filter(or >= 2, p_adjusted < 0.05)

  ## visualization of the FDR-adjusted p values and OR
  ## for the significantly enriched subsystems

  significant_subs_or_plot <- significant_subsystems %>%
    ggplot(aes(x = or,
               y = reorder(subsystem, or),
               size = -log10(p_adjusted),
               fill = reaction_set)) +
    facet_grid(reaction_set ~ .,
               scales = "free",
               space = "free") +
    geom_point(shape = 21) +
    scale_size_area(labels = function(x) signif(10^-(x), 2),
                    max_size = 4.5,
                    name = "pFDR") +
    scale_fill_manual(values = c("activated" = "firebrick",
                                 "inhibited" = "steelblue"),
                      name = "reaction set") +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    labs(title = "RECON metabolic subsystems, collagen hi vs low",
         x = "enrichment, OR")

  ## counts of differentially regulated reaction
  ## in the significantly enriched subsystems

  reaction_regulation %>%
    select(subsystems = significant_subsystems$subsystem) %>%
    count_regulated

  ## visualization of percentages of activated and inhibited reactions
  ## in the significantly enriched subsystems

  significant_subs_percent_plot <- reaction_regulation %>%
    select(subsystems = significant_subsystems$subsystem) %>%
    plot(plot_type = "numbers",
         scale = "percent",
         plot_title = "RECON metabolic subsystems, collagen hi vs low")

# END ---------
