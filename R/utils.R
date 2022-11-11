# Functions used in extraction of the gene and reaction identifiers
# and calculation of regulation estimates and their errors.

# Escape numeric in a string -----

#' Escape EntrezID in a string.
#'
#' @description Finds EntrezID in a string and escapes them with '``.
#' @return a sting.
#' @param str a string.

  escape_numbers <- function(str) {

    stopifnot(is.character(str))

    numbs <- unlist(stringi::stri_extract_all(str, regex = '\\d+\\.\\d+'))

    for(i in numbs) {

      rpl <- stringi::stri_replace(i, regex = '\\.\\d+', replacement = '')

      str <-
        stringi::stri_replace_first_fixed(str,
                                          pattern = i,
                                          replacement = paste0('`', rpl, '`'))

    }

    str

  }

# Evaluation of the gene mapping rules ----

#' Evaluate gene mapping rules.
#'
#' @description Evaluates a list of gene mapping rules stored as expression.
#' Such rules are retrieved from the SBML database e.g.
#' with \code{\link{extract_genes}}.
#' @return a numeric vector with the reaction regulation estimates.
#' Named with reaction IDs.
#' @param x a vector with gene regulation estimates named
#' with Entrez ID identifiers.
#' @param rule_exp a list of gene mapping rules.
#' @param parent_env an environments containing evaluation functions and values
#' for genes missing from the x gene regulation estimate vector.

  eval_rule <- function(x, rule_exp, parent_env) {

    stopifnot(!is.null(names(rule_exp)))
    stopifnot(rlang::is_environment(parent_env))
    stopifnot(!is.null(names(x)))
    stopifnot(is.numeric(x))

    ## entry control is realized by the upstream functions.

    eval_env <- rlang::new_environment(data = as.list(x),
                                       parent = parent_env)


    purrr::map_dbl(rule_exp,
                   eval,
                   envir = eval_env)


  }

# Simulation of gene regulation estimate distribution -------

#' Draw values from the gene regulation estimate distribution.
#'
#' @description Draws values from the gene regulation estimate
#' normal distribution given the standard deviation and expected value.
#' @param x a numeric vector with fold-regulation estimates (not log2!!!).
#' Elements are named with Entrez IDs.
#' @param err a numeric vector with an error estimate such as SD or SEM
#' (not log2!!!), named with Entrez IDs.
#' @param n number of values to be drawn.
#' @param seed seed of the random number generator.
#' @param .parallel logical, should the function work in parallel? Works only if
#' TRUE and a parallel backend compatible with furrr::future_map() is provided.
#' @param return a list with n numeric values corresponding to the estimates,
#' named with the Entrez IDs.

  draw_norm <- function(x, err, n, seed = NULL, .parallel = FALSE) {

    ## entry control

    stopifnot(!is.null(names(x)))
    stopifnot((!is.null(names(err))))
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(err))
    stopifnot(is.numeric(n))

    n <- integer(n)

    stopifnot(n < 1)

    cmm_genes <- intersect(names(x), names(err))

    if(length(cmm_genes) == 0) {

      stop('Names of x and err do not fit. Sure, they are the right vectors?',
           call. = FALSE)

    }

    x <- x[cmm_genes]
    err <- err[cmm_genes]

    if(!is.null(seed)) {

      set.seed(seed)

    }

    ## generating the distributions

    if(.parallel) {

      distr_lst <-
        furrr::future_map2(x, err,
                           ~as.list(rnorm(n, mean = .x, sd = .y)),
                           .options = furrr::furrr_options(seed = TRUE))

    } else {

      distr_lst <- purrr::map2(x, err,
                               ~as.list(rnorm(n, mean = .x, sd = .y)))

    }

    distr_lst <-
      purrr::map(distr_lst,
                 ~rlang::set_names(.x, paste0('draw_', 1:length(.x))))

    distr_lst <- purrr::transpose(distr_lst)

    purrr::map(distr_lst, unlist)

  }

# Testing of regulation significance -----

#' Regulation significance based on MC simulation results.
#'
#' @description Tests for gene regulation significance based on
#' Monte Carlo (MC) simulation results.
#' @param fold_reg MC-simulated fold-regulation values.
#' @param mu the reference value.
#' @return a numeric p value.

  test_mc <- function(fold_reg, mu = 1) {

    fold_reg <- fold_reg[!is.na(fold_reg)]

    n <- length(fold_reg)

    expected <- mean(fold_reg, na.rm = TRUE)


    if(is.na(expected)) {

      return(NA)

    }

    if(expected == mu) {

      return(1)

    }

    if(expected > mu) {

      fail_n <- length(fold_reg[fold_reg <= mu])

    } else {

      fail_n <- length(fold_reg[fold_reg >= mu])

    }

    if(fail_n == 0) fail_n <- 1

    return(fail_n/n)

  }

# Calculation of the of the regulation estimates -------

#' Calculate reaction regulation estimated based on gene expression.
#'
#' @description Computes regulation estimates for the reactions present
#' in the database based on the regulation estimates for gene expression.
#' As described in the seminal BiGGR package paper
#' (DOI 10.1371/journal.pone.0119016), it is assumed that the change
#' in gene expression corresponds directly to the magnitude of modulation of
#' the pathway, e.g. two-fold increase in expression of an enzyme leads to a
#' two-fold increase in activity of the respective pathway
#' (gene - protein - reaction or GPR principle). In case, more genes
#' are mapped to a given pathway, the fold regulation is computed accordingly
#' to the logical operators between the genes.
#' If a given reaction is not associate with any genes, its regulation estimate
#' is set as specified by the x_default argument.
#' If a vector with regulation estimate errors is provided, the cumulative error
#' will be returned as well - see Details.
#' @details Works in a comparable fashion to \code{\link[BiGGR]{gprMapping}}
#' but considerably faster.
#' When a vector of errors for gene regulation is provided, there are two ways
#' to compute errors for the pathway regulation:
#' * __from the normal distribution__: the cumulative reaction error is
#' calculated as a sum of errors for the genes mapped to the reaction,
#' z statistic is computed as (fold-regulation - 1)/error,
#' raw p values and 95% confidence intervals are obtained from the theoretical
#' distribution of z statistics and corrected for multiple testing with the
#' FDR/Benjamini-Hochberg method.
#' While this method is computationally fast, the error, t statistic,
#' confidence intervals and p values are just raw estimates of significance.
#' * __by Monte Carlo method__: this option is still pending.
#' @references Gavai et al. \url{https://doi.org/10.1371/journal.pone.0119016}.
#' @return a data frame with the following variables: reaction identifiers
#' and the fold regulation estimates. If err is not NULL, an error estimate per
#' reaction is included in the data frame as well along with the z statistic
#' (fold-regulation/error), 95% confidence intervals and
#' raw p-values (normal distribution) and FDR-corrected p-values.
#' If return_genes is set to TRUE, a gene-mapping table is returned as well in
#' a list with the fold regulation data frame.
#' @param x a numeric vector with fold-regulation estimates (not log2!!!).
#' Elements are named with Entrez IDs.
#' @param err a numeric vector with an error estimate such as SD or SEM
#' (not log2!!!), named with Entrez IDs.
#' @param database an object of class SBMLDocument.
#' @param or_fun one of 'mean' (default), 'median', 'min' or 'max'. Specifies the name of
#' the function used to handle the OR operator between the genes.
#' @param and_fun one of 'mean', 'median', 'min' (default) or 'max'. Specifies the name of
#' the function used to handle the OR operator between the genes.
#' @param x_default the default numeric value or NA for regulation of the genes
#' present in the database and absent from the regulation data. If NA,
#' the genes missing from x won't be included in subsequent SBML model creation.
#' @param return_genes logical, should a data frame with assignment of the
#' genes to reactions be returned? Defaults to FALSE.
#' @param return_mc logical, should all Monte Carlo-determined regulation
#' estimates for reactions be returned? Defaults to FALSE.
#' @param err_method specifies the method of calculation of errors for reactions
#' as described in Details: 'norm' (from normal distribution, default) or
#' 'mc' (Monte Carlo simulation).
#' @param n_iter number of iterations for the Monte carlo simulation of
#' the reaction errors.
#' @param ci_method method of calculation of confidence intervals:
#' percentile ('perc', default) or BCA ('bca').
#' @param seed a seed for the rangom number generator.

  get_regulation <- function(x,
                             err = NULL,
                             database,
                             or_fun = c('mean', 'median', 'min', 'max'),
                             and_fun = c('min', 'max', 'mean', 'median'),
                             x_default = 1,
                             return_genes = FALSE,
                             return_mc = FALSE,
                             err_method = c('norm', 'mc'),
                             n_iter = 1000,
                             ci_method = c('perc', 'bca'),
                             seed = NULL,
                             .parallel = TRUE) {

    ## entry control --------

    start_time <- Sys.time()

    on.exit(future::plan('sequential'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)), add = TRUE)

    if(!is.numeric(x)) {

      stop('x has to be a numeric vector with  Entrez IDs as names.',
           call. = FALSE)

    }

    if(is.null(names(x))) {

      stop('x has to be a numeric vector with  Entrez IDs as names.',
           call. = FALSE)

    }

    if(any(x < 0)) {

      stop('All elements of x have to be > 0.
           Are you sure, your values are in linear scale?',
           call. = FALSE)

    }

    or_fun <- match.arg(or_fun[1], c('mean', 'median', 'min', 'max'))
    and_fun <- match.arg(and_fun[1], c('min', 'max', 'mean', 'median'))
    err_method <- match.arg(err_method[1], c('norm', 'mc'))
    ci_method <- match.arg(ci_method[1], c('perc', 'bca'))

    stopifnot(is.numeric(n_iter))

    n_iter <- as.integer(n_iter)

    stopifnot(is.logical(return_genes))
    stopifnot(is.logical(return_mc))

    ## gene identifier data frame for the database -------

    message('Reading the gene mapping')

    gene_tbl <- extract_genes(database)

    ## building a parent evaluation environment -------

    message('Building an evaluation environment')

    o_fun <- switch(or_fun,
                    min = function(x, y) min(c(x, y), na.rm = TRUE),
                    max = function(x, y) max(c(x, y), na.rm = TRUE),
                    mean = function(x, y) mean(c(x, y), na.rm = TRUE),
                    median = function(x, y) median(c(x, y), na.rm = TRUE))

    a_fun <- switch(and_fun,
                    min = function(x, y) min(c(x, y), na.rm = TRUE),
                    max = function(x, y) max(c(x, y), na.rm = TRUE),
                    mean = function(x, y) mean(c(x, y), na.rm = TRUE),
                    median = function(x, y) median(c(x, y), na.rm = TRUE))


    all_genes <- purrr::reduce(gene_tbl$entrez_id, union)

    miss_genes <- all_genes[!all_genes %in% names(x)]

    miss_gene_lst <- rlang::set_names(rep(x_default, length(miss_genes)),
                                      miss_genes)
    fun_ev <-
      rlang::new_environment(data = c(list(`%AND%` = a_fun,
                                           `%OR%` = o_fun,
                                           `(` = `(`),
                                      as.list(miss_gene_lst)))

    ## serial evaluation: determination of regulation estimates ------

    message('Calculation of reaction regulation estimates')

    reg_lst <- eval_rule(x = x,
                         rule_exp = gene_tbl$exprs,
                         parent_env = fun_ev)

    reg_tbl <- tibble::tibble(react_id = names(reg_lst),
                              fold_reg = unlist(reg_lst))

    reg_tbl <- dplyr::filter(reg_tbl, !is.na(fold_reg))
    reg_tbl <- dplyr::filter(reg_tbl, !is.nan(fold_reg))
    reg_tbl <- dplyr::filter(reg_tbl, !is.infinite(fold_reg))

    ## filling with the default value if a reaction is absent
    ## from the regulation estimate table.

    all_reactions <- names(database@model@reactions)

    miss_reactions <- all_reactions[!all_reactions %in% reg_tbl$react_id]

    reg_tbl <- rbind(reg_tbl,
                     tibble::tibble(react_id = miss_reactions,
                                    fold_reg = x_default))

    if(is.null(err)) {

      if(!return_genes) {

        return(reg_tbl)

      } else {

        return(list(gene_map = gene_tbl,
                    reg = reg_tbl))

      }

    }

    ## error calculation: entry control ------

    if(!is.numeric(err)) {

      stop('err has to be a numeric vector with  Entrez IDs as names.',
           call. = FALSE)

    }

    if(is.null(names(err))) {

      stop('err has to be a numeric vector with  Entrez IDs as names.',
           call. = FALSE)

    }

    if(any(err < 0)) {

      stop('All elements of err have to be > 0.
           Are you sure, your values are in linear scale?',
           call. = FALSE)

    }

    ## computing the sum error from the normal distribution -------

    if(err_method == 'norm') {

      message('Calulation of the errors: normal distribution')

      err_lst <- purrr::map(gene_tbl$entrez_id, ~err[.x])

      err_lst <- purrr::map(err_lst, sum, na.rm = TRUE)

      err_tbl <- tibble::tibble(react_id = names(err_lst),
                                error = unlist(err_lst))

      err_tbl <- dplyr::filter(err_tbl, !is.na(error))

      reg_tbl <- dplyr::left_join(reg_tbl, err_tbl, by = 'react_id')

      reg_tbl <-
        dplyr::mutate(reg_tbl,
                      z = (fold_reg - 1)/error,
                      z = ifelse(is.infinite(z), NA, z),
                      lower_ci = fold_reg + error * qnorm(0.025),
                      upper_ci = fold_reg + error * qnorm(0.975),
                      p_value = pnorm(abs(z), lower.tail = FALSE) * 2,
                      p_adjusted = p.adjust(p_value, 'BH'))

      reg_tbl <-
        reg_tbl[c('react_id', 'fold_reg', 'error',
                  'lower_ci', 'upper_ci',
                  'z', 'p_value', 'p_adjusted')]

      if(!return_genes & !return_mc) {

        return(reg_tbl)

      } else {

        if(!return_mc) {

          return(list(gene_map = gene_tbl,
                      reg = reg_tbl))

        } else {

          return(list(gene_map = gene_tbl,
                      reg = reg_tbl,
                      mc = NULL))

        }

      }

    }

    ## computing the errors: MC -------

    message('Simulation of gene regulation estimates')

    ### CI computing functions and parallel backend

    perc_fun <- function(x) {

      c(error = sd(x, na.rm = TRUE),
        rlang::set_names(quantile(x,
                                  c(0.025, 0.975),
                                  na.rm = TRUE),
                         c('lower_ci', 'upper_ci')),
        p_value = test_mc(x))

    }

    bca_fun <- function(x) {

      c(error = sd(x, na.rm = TRUE),
        rlang::set_names(coxed::bca(x[!is.na(x)],
                                    conf.level = 0.95),
                         c('lower_ci', 'upper_ci')),
        p_value = test_mc(x))

    }

    ci_fun <-
      switch(ci_method,
             perc = perc_fun,
             bca = bca_fun)

    if(.parallel) {

      future::plan('multisession')

    }

    ## simulation of the gene regulation values

    if(n_iter > 500 & .parallel) {

      gene_est_lst <- draw_norm(x = x,
                                err = err,
                                n = n_iter,
                                seed = seed,
                                .parallel = TRUE)

    } else {

      gene_est_lst <- draw_norm(x = x,
                                err = err,
                                n = n_iter,
                                seed = seed,
                                .parallel = FALSE)

    }

    ## calculation of the errors

    message('Error calculation')

    if(.parallel) {

      react_est_lst <-
        furrr::future_map(gene_est_lst,
                          biggrExtra:::eval_rule,
                          rule_exp = gene_tbl$exprs,
                          parent_env = fun_ev,
                          .options = furrr::furrr_options(seed = TRUE,
                                                          globals = c('eval_rule',
                                                                      'fun_ev',
                                                                      'gene_tbl'),
                                                          packages = c('generics',
                                                                       'rlang',
                                                                       'biggrExtra')))

    } else {

      react_est_lst <-
        purrr::map(gene_est_lst,
                   eval_rule,
                   rule_exp = gene_tbl$exprs,
                   parent_env = fun_ev)

    }

    react_est_lst <- purrr::map(react_est_lst,
                                ~ifelse(.x < 0, 0, .x))

    react_est_mtx <- do.call('rbind', react_est_lst)

    if(n_iter > 1000 & .parallel) {

      err_tbl <-
        furrr::future_map(colnames(react_est_mtx),
                          ~ci_fun(react_est_mtx[, .x]),
                          .options = furrr::furrr_options(seed = TRUE,
                                                          packages = c('biggrExtra',
                                                                       'rlang',
                                                                       'generics'),
                                                          globals = c('react_est_mtx',
                                                                      'ci_fun',
                                                                      'test_mc')))

    } else {

      err_tbl <- purrr::map(colnames(react_est_mtx),
                            ~ci_fun(react_est_mtx[, .x]))

    }

    err_tbl <- do.call('rbind', err_tbl)

    rownames(err_tbl) <- colnames(react_est_mtx)

    err_tbl <-
      tibble::rownames_to_column(as.data.frame(err_tbl), 'react_id')

    reg_tbl <- left_join(reg_tbl, err_tbl, by = 'react_id')

    reg_tbl <- dplyr::mutate(reg_tbl,
                             p_adjusted = p.adjust(p_value, 'BH'))

    reg_tbl <-
      reg_tbl[c('react_id', 'fold_reg', 'error',
                'lower_ci', 'upper_ci',
                'p_value', 'p_adjusted')]

    if(!return_genes & !return_mc) {

      return(reg_tbl)

    } else {

      if(!return_mc) {

        return(list(gene_map = gene_tbl,
                    reg = reg_tbl))

      } else {

        return(list(gene_map = gene_tbl,
                    reg = reg_tbl,
                    mc = react_est_mtx))

      }

    }

  }

# END ------
