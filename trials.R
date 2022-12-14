# some tools -----

  library(tidyverse)
  library(rlang)
  library(stringi)
  library(biggrExtra)

# Finding and escaping numbers in a string ------

  biggrExtra:::escape_numbers('(75.1) %OR% (275.1) %AND% (1738.5) %AND% (275.1) %AND% (2653.1)')

# Extracting gene IDs -----

  data("Recon2")

  gene_info <- biggrExtra:::extract_genes(Recon2)

# Expression regulation estimates and their errors------

  de_vec <- rlang::set_names(tcga_data$estimate,
                             tcga_data$entrez_id)

  de_vec <- 2^de_vec

  err_vec <- rlang::set_names(tcga_data$se,
                              tcga_data$entrez_id)

# Reaction regulation estimates ------

  reaction_vec <- biggrExtra:::get_regulation(de_vec,
                                              err_vec,
                                              Recon2,
                                              or_fun = 'min',
                                              and_fun = 'max')

# Generation of a model --------

  reco2_model <- build_geneSBML(x = de_vec,
                                err = err_vec,
                                database = Recon2,
                                x_default = 1)

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

  annotate_bigg('R_A4GNT1g')

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

  #max_proc <- "R_NDPK1m - R_HEX1 - R_PFK - R_PGK + R_PYK"

  #check_geneSBML(c('R_ATPS4m', 'R_NDPK1m', 'R_HEX1', 'R_LACD', 'R_LDH_D'), reco2_model)

  #ext_metabolites <- c('M_g6p_r')

  #fit(reco2_model,
     # maximize = max_proc,
      #externals = ext_metabolites,
      #file.name = 'test.lim')

  #test_reg <- reco2_model$reg %>%
    #top_n(10, fold_reg)

  #createLIMFromSBML(reco2_model$model,
                   # maximize = max_proc,
                    #equations = as.list(test_reg[c(1, 2)]),
                    #externals = ext_metabolites,
                   # file.name = 'test.lim')

  #test_rates <- getRates('test.lim')














# drawing from a normal distribution ------

  test_draw <- biggrExtra:::draw_norm(x = de_vec, err = err_vec, n = 10)

# Model fitting with Monte Carlo errors ------

  reco2_mc_model <- build_geneSBML(x = de_vec,
                                   err = err_vec,
                                   database = Recon2,
                                   x_default = 1,
                                   err_method = 'mc',
                                   n_iter = 1000,
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

# Visualization ------

  test1 <-
    visualize(reco2_mc_model,
              rate_sep = ', ',
              suffixes = 'error',
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

  steroid_rcts <- reactions %>%
    filter(stri_detect(reaction_string,
                       regex = 'andrstndn_r|tststerone_r|estradiol_r')) %>%
    .$bigg_id %>%
    paste0('R_', .)

  steroid_metabs <- annotate_bigg(steroid_rcts,
                                  value = 'reaction_string',
                                  annotation_db = reactions) %>%
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


# END -----

