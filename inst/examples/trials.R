# some tools -----

  library(tidyverse)
  library(rlang)
  library(stringi)
  library(BiGGR)
  library(biggrExtra)

# Finding and escaping numbers in a string ------

  biggrExtra:::escape_numbers('(75.1) %OR% (275.1) %AND% (1738.5) %AND% (275.1) %AND% (2653.1)')

# Extracting gene IDs -----

  data("Recon2")

  gene_info <- extract_genes(Recon2D)

  gene_info %>%
    filter(react_id == 'R_ATPS4m') %>%
    .$entrez_id

# Patching of the Recon2D errors -------

  #Recon2D <- Recon2

  #atp_rule <- read_file('recon1_rule.txt') %>%
    #stri_replace_all(fixed = '_AT', replacement = '.')

  #atp_rule <- paste0('<p>GENE_ASSOCIATION: ',
    #                 atp_rule,
    #                 '</p>')

  #Recon2D@model@reactions$R_ATPS4m@notes <-
    #Recon2D@model@reactions$R_ATPS4m@notes %>%
    #stri_replace(regex = '<p>GENE_ASSOCIATION:\\s{1}.*</p>',
       #          replacement = atp_rule)

  #save(Recon2D, file = './data/Recon2D.RData')

# Expression regulation estimates and their errors------

  de_vec <- rlang::set_names(tcga_data$estimate,
                             tcga_data$entrez_id)

  de_vec <- 2^de_vec

  err_vec <- rlang::set_names(tcga_data$se,
                              tcga_data$entrez_id)

  de_native <- rlang::set_names(tcga_data$estimate,
                                tcga_data$entrez_id)

  err_native <- rlang::set_names(tcga_data$se,
                                 tcga_data$entrez_id)

# Generation of a model --------

  reco2_model <- build_geneSBML(x = de_native,
                                #err = err_vec,
                                scale = 'log2',
                                database = Recon2D,
                                x_default = 1,
                                save_memory = FALSE)

  plot(reco2_model, type = 'top', n_top = 20)
  plot(reco2_model, type = 'regulation')
  plot(reco2_model, type = 'errors')
  plot(reco2_model, type = 'volcano', n_top = 20) +
    theme(plot.tag.position = 'bottom')

  plot(reco2_model, type = 'top_forest', n_top = 20)

# Annotation ------

  reco2_model$reg$react_id %>%
    annotate_bigg()

  reco2_top <- plot(reco2_model,
                    type = 'top') +
    scale_y_discrete(labels = annotate_bigg(reco2_model$reg$react_id,
                                            annotation_db = reco2_model))

  annotate_bigg('R_A4GNT1g', annotation_db = Recon2D)

  annotate_bigg('R_A4GNT1g', annotation_db = reco2_model)

# Extraction and checking for reactions and metabolites -------

  calculate(reco2_model, c("R_HEX1", "R_PGI", "R_PFK", "R_FBA", "R_TPI",
                           "R_GAPD", "R_PGK", "R_PGM", "R_ENO", "R_PYK",
                           "R_G6PDH2r", "R_PGL", "R_GND", "R_RPE", "R_RPI",
                           "R_TKT1", "R_LDH_D", "R_LDH_Lm", 'R_r0358',
                           'R_r0364', 'R_r0363', 'R_r0355'))

  check_geneSBML(geneSBML = reco2_model,
                 react_id = c("R_HEX1", "R_PGI", "R_PFK", "R_FBA", "R_TPI",
                              "R_GAPD", "R_PGK", "R_PGM", "R_ENO", "R_PYK",
                              "R_G6PDH2r", "R_PGL", "R_GND", "R_RPE", "R_RPI",
                              "R_TKT1", "R_LDH_D", "R_PDHm", "R_r0358", 'R_r0363', 'R_r0354'),
                 metab_id = c("M_glc_D_c", "M_g6p_c", "M_f6p_c",
                              "M_fdp_c", "M_dhap_c", "M_g3p_c",
                              "M_13dpg_c", "M_3pg_c", "M_2pg_c",
                              "M_pep_c", "M_pyr_c", "M_pyr_m",
                              "M_6pgl_c", "M_6pgc_c", "M_ru5p_D_c",
                              "M_xu5p_D_c", "M_r5p_c", "M_g3p_c",
                              "M_s7p_c", "M_lac_D_c", "M_lac_D_m", 'M_gam6p_c'))

  components(reco2_model, 'reactions')[stri_detect(components(reco2_model,
                                                                'reactions'),
                                                     fixed = 'PDH')]

  components(reco2_model, 'genes') %>% length

# LIM modeling -----

  ## Baustelle!!!

 # max_proc <- "R_ATPS4m + R_NDPK1m - R_HEX1 - R_PFK - R_PGK + R_PYK"

  #check_geneSBML(c('R_ATPS4m', 'R_NDPK1m', 'R_HEX1', 'R_LACD', 'R_LDH_D'),
             #    geneSBML = reco2_mc_model)

#  ext_metabolites <- c("M_glc_D_e", "M_lac_L_e", "M_ala_L_e",
#                       "M_gln_L_e", "M_h2o_e", "M_co2_e",
 #                      "M_o2_e", "M_h_e", "M_o2s_m",
  #                     "M_adp_c", "M_atp_c", "M_pi_c",
   #                    "M_h_c", "M_nadp_c", "M_nadph_c",
    #                   "M_na1_c", "M_na1_e", "M_gln_L_c",
     #                  "M_nh4_c", "M_pyr_e")

#  check_geneSBML(metab_id = ext_metabolites,
 #                geneSBML = reco2_mc_model)

  #fit(reco2_mc_model,
   #   maximize = max_proc,
    #  externals = ext_metabolites,
     # file.name = 'test.lim')

  #test_rates <- getRates('test.lim')

# Model fitting with Monte Carlo errors ------

  reco2_mc_model <- build_geneSBML(x = de_native,
                                   err = err_native,
                                   scale = 'log2',
                                   database = Recon2D,
                                   x_default = 1,
                                   err_method = 'mc',
                                   n_iter = 1000,
                                   mc_estimate = 'median',
                                   burn_in = 100,
                                   ci_method = 'bca',
                                   .parallel = TRUE)

  plot(reco2_mc_model, type = 'regulation', signif_type = 'fdr')
  plot(reco2_mc_model, type = 'top', label_names = TRUE)

  plot(reco2_mc_model, type = 'errors')
  plot(reco2_mc_model, type = 'volcano') +
    theme(plot.tag.position = 'bottom')

  plot(reco2_mc_model,
       type = 'top_forest',
       n_top = 10,
       label_names = FALSE,
       signif_type = 'fdr')

  reco2_mc_model$reg$lower_ci %>%
    range(na.rm = TRUE)

  reco2_mc_model$reg %>%
    filter(react_id == 'R_ATPS4m')

  reco2_mc_model$gene_map %>%
    filter(react_id == 'R_ATPS4m') %>%
    .$entrez_id

# Visualization ------

  test1 <-
    visualize(reco2_mc_model,
              rate_sep = ', ',
              suffixes = 'none',
              suffix_sep = '\u00B1',
              signif_type = 'fdr',
              relevant.species = c("M_glc_D_c", "M_g6p_c", "M_f6p_c",
                                   "M_fdp_c", "M_dhap_c", "M_g3p_c",
                                   "M_13dpg_c", "M_3pg_c", "M_2pg_c",
                                   "M_pep_c", "M_pyr_c", "M_pyr_m",
                                   "M_6pgl_c", "M_6pgc_c", "M_ru5p_D_c",
                                   "M_xu5p_D_c", "M_r5p_c", "M_g3p_c",
                                   "M_s7p_c", "M_lac_D_c", "M_lac_D_m",
                                   "M_fru_c", 'M_gam_c',
                                   'M_man6p_c', 'M_man_c', 'M_itp_c'),
              relevant.reactions = c("R_HEX1", "R_PGI", "R_PFK",
                                     "R_FBA", "R_TPI",
                                     "R_GAPD", "R_PGK", "R_PGM",
                                     "R_ENO", "R_PYK",
                                     "R_G6PDH2r", "R_PGL", "R_GND",
                                     "R_RPE", "R_RPI",
                                     "R_TKT1", "R_LDH_D", "R_LDH_Lm",
                                     "R_PDHm", "R_r0358", "R_r0364",
                                     'R_r0361', 'R_r0363', 'R_r0357'),
              layoutType = "dot",
              plt.margins = c(10, 50, 50, 10))


  plot(test1)

  test2 <-
    visualize(reco2_mc_model,
              rate_sep = '\n',
              suffixes = 'p_fdr',
              suffix_sep = ', p = ',
              signif_type = 'fdr',
              relevant.species = c("M_cit_m", "M_icit_m" , "M_akg_m",
                                   "M_succoa_m", "M_succ_m", "M_fum_m",
                                   "M_mal_L_m", "M_oaa_m"),
              relevant.reactions = c("R_CSm", "R_ACONTm", "R_ICDHxm",
                                     "R_AKGDm", "R_SUCOAS1m", "R_SUCD1m",
                                     "R_FUMm", "R_MDHm", "R_ICDHyrm", "R_ME1m",
                                     "R_ME2m", "R_ASPTAm","R_AKGMALtm", "R_GLUDym",
                                     "R_ABTArm", "R_SSALxm","R_CITtam"),
              layoutType = "circo",
              plt.margins = c(150, 235, 150, 230))

  plot(test2)

  steroid_rcts <- biggrExtra::reactions %>%
    filter(stri_detect(reaction_string,
                       regex = 'andrstndn_r|tststerone_r|estradiol_r')) %>%
    .$bigg_id %>%
    paste0('R_', .)

  steroid_metabs <- annotate_bigg(steroid_rcts,
                                  value = 'reaction_string',
                                  annotation_db = biggrExtra::reactions) %>%
    map(stri_split_regex,
        pattern = '(\\s{1}<->\\s{1})|(\\s{1}\\+\\s{1})') %>%
    unlist %>%
    unname %>%
    unique

  steroid_metabs <- steroid_metabs[!steroid_metabs %in% c('h_r', 'nadph_r',
                                                          'nadp_r', 'nad_r',
                                                          'nadh_r', '0.5 o2_r',
                                                          'h2o_r', '2.0 o2_r',
                                                          'o2_r', 'h_c',
                                                          '2.0 h2o_r',
                                                          'udp_r', 'ac_r',
                                                          'for_r',
                                                          'udpglcur_r',
                                                          'acald_r')]

  steroid_metabs <- paste0('M_', steroid_metabs)

  test3 <-
    visualize(reco2_mc_model,
              rate_sep = '\n',
              suffixes = 'error',
              suffix_sep = ' \u00B1 ',
              signif_type = 'fdr',
              fontsize = 0.65,
              relevant.species = steroid_metabs,
              relevant.reactions = steroid_rcts,
              layoutType = "dot",
              lwd.min = 2,
              plt.margins = c(20, 50, 20, 150),
              colors = c(activated = 'coral4',
                         inhibited = 'darkolivegreen',
                         ns = 'gray50'))

  plot(test3)

  plot(reco2_mc_model,
       type = 'forest',
       relevant.reactions = steroid_rcts,
       label_name = TRUE)

  signif_steroid_react <- components(reco2_mc_model, 'regulation') %>%
    filter(react_id %in% steroid_rcts,
           p_adjusted < 0.05)

  plot(reco2_mc_model,
       type = 'mc',
       relevant.reactions = signif_steroid_react$react_id,
       line_alpha = 0.25) +
    facet_wrap(nrow = 3, facets = 'react_id') +
    theme(legend.position = 'none')


# Gene, reaction and metabolite mapping -------

  gene_annot <- react_to_gene(c("R_HEX1", "R_PGI", "R_PFK",
                                "R_FBA", "R_TPI",
                                "R_GAPD", "R_PGK", "R_PGM",
                                "R_ENO", "R_PYK",
                                "R_G6PDH2r", "R_PGL", "R_GND",
                                "R_RPE", "R_RPI",
                                "R_TKT1", "R_LDH_D", "R_LDH_Lm",
                                "R_PDHm", "R_r0358", "R_r0364",
                                'R_r0361', 'R_r0363', 'R_r0357'),
                              reco2_mc_model)

  react_annot <- gene_to_react(c('1738', '8050'), reco2_mc_model)

  metab_annot <- react_to_metab(c("R_HEX1", "R_PGI", "R_PFK",
                                  "R_FBA", "R_TPI",
                                  "R_GAPD", "R_PGK", "R_PGM",
                                  "R_ENO", "R_PYK",
                                  "R_G6PDH2r", "R_PGL", "R_GND",
                                  "R_RPE", "R_RPI",
                                  "R_TKT1", "R_LDH_D", "R_LDH_Lm",
                                  "R_PDHm", "R_r0358", "R_r0364",
                                  'R_r0361', 'R_r0363', 'R_r0357'),
                                exc_regex = '(^h_)|(^pi_)')

# New annotation functions -------

  react_to_metab(c('R_ATPS4m', 'R_PFK'),
                 detailed = TRUE)

  extract_subsystems(Recon2D, as_list = FALSE)

# Enhanced output of the components and plot method ------

  components(reco2_mc_model, type = 'gene_map')

  components(reco2_mc_model, type = 'regulation')

  components(reco2_mc_model, type = 'subsystems')

  plot(reco2_mc_model,
       type = 'bar',
       relevant.reactions = c("R_HEX1", "R_PGI", "R_PFK",
                              "R_FBA", "R_TPI",
                              "R_GAPD", "R_PGK", "R_PGM",
                              "R_ENO", "R_PYK",
                              "R_G6PDH2r", "R_PGL", "R_GND",
                              "R_RPE", "R_RPI",
                              "R_TKT1", "R_LDH_D", "R_LDH_Lm",
                              "R_PDHm", "R_r0358", "R_r0364",
                              'R_r0361', 'R_r0363', 'R_r0357'),
       label_names = TRUE,
       order_reactions = TRUE)

  plot(reco2_mc_model,
       type = 'volcano',
       relevant.reactions = c("R_HEX1", "R_PGI", "R_PFK",
                              "R_FBA", "R_TPI",
                              "R_GAPD", "R_PGK", "R_PGM",
                              "R_ENO", "R_PYK",
                              "R_G6PDH2r", "R_PGL", "R_GND",
                              "R_RPE", "R_RPI",
                              "R_TKT1", "R_LDH_D", "R_LDH_Lm",
                              "R_PDHm", "R_r0358", "R_r0364",
                              'R_r0361', 'R_r0363', 'R_r0357'),
       label_names = TRUE,
       order_reactions = TRUE)


# Reaction counts and enrichment analysis for the subsystems -------

  test_fisher_data <- reco2_mc_model %>%
    suba(signif_type = 'fdr',
         method = 'fisher',
         .parallel = FALSE)

  test_fisher_data %>%
    filter(status == 'activated') %>%
    ggplot(aes(x = frac_total,
               y = -log10(p_adjusted))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dashed')

  test_fisher_data %>%
    filter(p_adjusted < 0.05)

  test_simulation_data <- reco2_mc_model %>%
    suba(signif_type = 'fdr',
         method = 'simulation',
         n_iter = 1000,
         .parallel = FALSE)

  test_simulation_data %>%
    filter(status == 'inhibited') %>%
    ggplot(aes(x = OR,
               y = -log10(p_adjusted))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dashed')

  test_simulation_data %>%
    filter(p_adjusted < 0.05,
           status == 'inhibited')

# Memory saver subclass -------

  reco2_saver_model <- build_geneSBML(x = de_native,
                                   err = err_native,
                                   scale = 'log2',
                                   database = Recon2D,
                                   x_default = 1,
                                   err_method = 'mc',
                                   n_iter = 1000,
                                   mc_estimate = 'median',
                                   burn_in = 100,
                                   ci_method = 'bca',
                                   save_memory = TRUE,
                                   .parallel = TRUE)

  is_geneSBML(reco2_saver_model)

