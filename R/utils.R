# Functions used in extraction of the gene and reaction identifiers
# and calculation of regulation estimates and their errors.

# Escape numeric in a string -----

#' Escape EntrezID in a string.
#'
#' @description
#' Finds EntrezID in a string and escapes them with '``.
#'
#' @return a sting.
#'
#' @param str a string.

  escape_numbers <- function(str) {

    stopifnot(is.character(str))

    numbs <- unlist(stri_extract_all(str, regex = '\\d+\\.\\d+'))

    for(i in numbs) {

      rpl <- stri_replace(i, regex = '\\.\\d+', replacement = '')

      str <-
        stri_replace_first_fixed(str,
                                 pattern = i,
                                 replacement = paste0('`', rpl, '`'))

    }

    str

  }

# Evaluation of the gene mapping rules ----

#' Evaluate gene mapping rules.
#'
#' @description
#' Evaluates a list of gene mapping rules stored as expression.
#' Such rules are retrieved from the SBML database e.g.
#' with \code{\link{extract_genes}}.
#'
#' @return a numeric vector with the reaction regulation estimates.
#' Named with reaction IDs.
#'
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


    map_dbl(rule_exp, eval, envir = eval_env)


  }

# Simulation of gene regulation estimate distribution -------

#' Draw values from the gene regulation estimate distribution.
#'
#' @description
#' Draws values from the gene regulation estimate
#' normal distribution given the standard deviation and expected value.
#'
#' @param x a numeric vector with fold-regulation estimates.
#' Elements are named with Entrez IDs.
#' @param err a numeric vector with an error estimate such as SD or SEM,
#' named with Entrez IDs.
#' @param n number of values to be drawn.
#' @param scale name of the scale of the provided regulation estimates and
#' their errors: 'identity' (default), 'log2' or 'log'.
#' @param seed seed of the random number generator.
#'
#' @return a numeric matrix: genes are columns, theri regulation estimates are
# rows.

  draw_norm <- function(x,
                        err,
                        n,
                        scale = c('identity', 'log2', 'log'),
                        seed = NULL) {

    ## entry control -----

    stopifnot(!is.null(names(x)))
    stopifnot((!is.null(names(err))))
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(err))
    stopifnot(is.numeric(n))

    stopifnot(n > 1)

    cmm_genes <- intersect(names(x), names(err))

    if(length(cmm_genes) == 0) {

      stop('Names of x and err do not fit. Sure, they are the right vectors?',
           call. = FALSE)

    }

    x <- x[cmm_genes]
    err <- err[cmm_genes]

    scale <- match.arg(scale[1], c('identity', 'log2', 'log'))

    if(!is.null(seed)) set.seed(seed)

    ## generating the distributions -----------

    distr_mtx <- drawNorm(mu = x,
                          sd = err,
                          n = n)

    dimnames(distr_mtx) <-
      list(x = paste0('iter_', 1:n),
           y = names(x))

    if(scale == 'identity') return(distr_mtx)

    trans_fun <-
      switch(scale,
             log2 = function(x) 2^x,
             log = exp)

    trans_fun(distr_mtx)

  }

# Calculation of the of the regulation estimates -------

#' Calculate reaction regulation estimated based on gene expression.
#'
#' @description
#' Intended for internal use, see: \code{\link{build_geneSBML}}.
#'
#' @return
#' a data frame with the following variables: reaction identifiers
#' and the fold regulation estimates. If err is not NULL, an error estimate per
#' reaction is included in the data frame as well along with the z statistic
#' (fold-regulation/error), 95% confidence intervals and
#' raw p-values (normal distribution) and FDR-corrected p-values.
#' If return_genes is set to TRUE, a gene-mapping table is returned as well in
#' a list with the fold regulation data frame.
#'
#' @param x a numeric vector with fold-regulation estimates. Elements are named
#' with Entrez IDs.
#' @param err a numeric vector with an error estimate such as SD or SEM, named
#' with Entrez IDs.
#' @param database an object of class `SBMLDocument` providing the gene -
#' reaction association rules.
#' @param scale specifies the scale of the provided expression regulation and
#' error estimates.
#' @param or_fun one of 'mean' (default), 'median', 'min' or 'max'. Specifies
#' the name of the function used to handle the `OR` operator between the genes in
#' the association rule.
#' @param and_fun one of 'mean', 'median', 'min' (default) or 'max'. Specifies
#' the name of the function used to handle the `OR` operator between the genes in
#' the association rule.
#' @param x_default the default numeric value or `NA` for regulation of the genes
#' present in the database and absent from the regulation data. If `NA`,
#' the genes missing from x won't be included in subsequent SBML model.
#' @param return_genes logical, should a data frame with assignment of the
#' genes to reactions be returned? Defaults to `FALSE`.
#' @param return_mc logical, should all Monte Carlo-determined regulation
#' estimates for reactions be returned? Defaults to `FALSE`.
#' @param err_method specifies the method of calculation of errors for reactions
#' as described in Details: 'norm' (from normal distribution, default) or
#' 'mc' (Monte Carlo simulation).
#' @param n_iter number of iterations for the Monte Carlo simulation of
#' the reaction errors.
#' @param mc_estimate statistic used to compute the reaction activity estimate
#' in Monte Carlo simulation: either mean or median over algorithm iterations.
#' Ignored if no errors were provided or `err_method = 'norm'`.
#' @param burn_in number of Monte Carlo iterations to be discarded prior to
#' computation of activity estimates and activity regulation errors. Ignored
#' if no errors were provided or `err_method = 'norm'`.
#' @param ci_method method of calculation of confidence intervals:
#' percentile ('perc', default) or BCA ('bca').
#' @param seed a seed for the random number generator.
#' @param .parallel logical, should th computation be run in parallel?

  get_regulation <- function(x,
                             err = NULL,
                             database = biggrExtra::Recon2D,
                             scale = c('identity', 'log2', 'log'),
                             or_fun = c('mean', 'median', 'min', 'max'),
                             and_fun = c('min', 'max', 'mean', 'median'),
                             x_default = 1,
                             return_genes = FALSE,
                             return_mc = FALSE,
                             err_method = c('norm', 'mc'),
                             n_iter = 1000,
                             mc_estimate = c('mean', 'median'),
                             burn_in = 0,
                             ci_method = c('perc', 'bca'),
                             seed = NULL,
                             .parallel = TRUE) {

    ## entry control --------

    start_time <- Sys.time()

    fold_reg <- NULL
    error <- NULL
    z <- NULL
    p_value <- NULL

    on.exit(plan('sequential'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)), add = TRUE)

    if(!is.numeric(x)) {

      stop('x has to be a numeric vector with  Entrez IDs as names.',
           call. = FALSE)

    }

    if(is.null(names(x))) {

      stop('x has to be a numeric vector with  Entrez IDs as names.',
           call. = FALSE)

    }

    scale <- match.arg(scale[1], c('identity', 'log2', 'log'))

    scale_fun <-
      switch(scale,
             identity = identity,
             log2 = function(x) 2^x,
             log = exp)

    nat_x <- x

    x <- scale_fun(x)

    if(scale == 'identity') {

      if(any(x < 0)) {

        stop(paste('All elements of x have to be >= 0. Are you sure,',
                   'your values are on the identity scale?'),
             call. = FALSE)

      }

    }

    or_fun <- match.arg(or_fun[1], c('mean', 'median', 'min', 'max'))
    and_fun <- match.arg(and_fun[1], c('min', 'max', 'mean', 'median'))
    err_method <- match.arg(err_method[1], c('norm', 'mc'))
    ci_method <- match.arg(ci_method[1], c('perc', 'bca'))

    stopifnot(is.numeric(n_iter))

    n_iter <- as.integer(n_iter)

    mc_estimate <- match.arg(mc_estimate[1], c('mean', 'median'))

    stopifnot(is.numeric(burn_in))
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

    all_genes <- reduce(gene_tbl$entrez_id, union)

    miss_genes <- all_genes[!all_genes %in% names(x)]

    miss_gene_lst <- set_names(rep(x_default, length(miss_genes)),
                                      miss_genes)
    fun_ev <-
      rlang::new_environment(data = c(list(`%AND%` = a_fun,
                                           `%OR%` = o_fun,
                                           `(` = `(`),
                                      as.list(miss_gene_lst)))

    ## serial evaluation: determination of regulation estimates ------
    ## only in case there are no errors or errors computed from the normal
    ## distribution

    if(is.null(err) | err_method == 'norm') {

      message('Calculation of reaction regulation estimates')

      reg_lst <- eval_rule(x = x,
                           rule_exp = gene_tbl$exprs,
                           parent_env = fun_ev)

      reg_tbl <- tibble(react_id = names(reg_lst),
                        fold_reg = unlist(reg_lst))

      reg_tbl <- filter(reg_tbl, !is.na(fold_reg))
      reg_tbl <- filter(reg_tbl, !is.nan(fold_reg))
      reg_tbl <- filter(reg_tbl, !is.infinite(fold_reg))

      ## filling with the default value if a reaction is absent
      ## from the regulation estimate table.

      all_reactions <- names(database@model@reactions)

      miss_reactions <- all_reactions[!all_reactions %in% reg_tbl$react_id]

      reg_tbl <- rbind(reg_tbl,
                       tibble(react_id = miss_reactions,
                              fold_reg = x_default))

      if(is.null(err)) {

        if(!return_genes) {

          return(reg_tbl)

        } else {

          return(list(gene_map = gene_tbl,
                      reg = reg_tbl))

        }

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

    if(scale == 'identity') {

      if(any(err < 0)) {

        stop(paste('All elements of err have to be >= 0.',
                   'Are you sure, your values are on the identity scale?'),
             call. = FALSE)

      }

    }

    ## computing the sum error from the normal distribution -------

    if(err_method == 'norm') {

      message('Calulation of the errors: normal distribution')

      err <- scale_fun(err)

      err_lst <- map(gene_tbl$entrez_id, ~err[.x])

      err_lst <- map(err_lst, sum, na.rm = TRUE)

      err_tbl <- tibble(react_id = names(err_lst),
                        error = unlist(err_lst))

      err_tbl <- filter(err_tbl, !is.na(error))

      reg_tbl <- left_join(reg_tbl, err_tbl, by = 'react_id')

      reg_tbl <-
        mutate(reg_tbl,
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

    ## computing the regulation estimates and errors: MC -------

    message('Regulation estimate and error calculation')

    ### CI computing functions
    ### simulation of the gene regulation values

    sim_est <- draw_norm(nat_x,
                         err = err,
                         n = n_iter,
                         scale = scale,
                         seed = seed)

    iters <- rownames(sim_est)

    ## calculation of the estimates and errors

    if(.parallel) plan('multisession')

    exp_globals <- c('eval_rule',
                     'fun_ev',
                     'gene_tbl')

    exp_packages <- c('generics',
                      'rlang',
                      'biggrExtra')

    react_est_lst <-
      future_map(iters,
                 ~eval_rule(x = sim_est[.x, ],
                            rule_exp = gene_tbl$exprs,
                            parent_env = fun_ev),
                 .options = furrr_options(seed = TRUE,
                                          globals = exp_globals,
                                          packages = exp_packages))


    react_est_lst <- map(react_est_lst, ~ifelse(.x < 0, 0, .x))

    react_est_mtx <- do.call('rbind', react_est_lst)

    if(burn_in > 0) {

      react_est_mtx <-
        react_est_mtx[(burn_in + 1):nrow(react_est_mtx), ]

    }

    ci_fun = switch(ci_method,
                    perc = perci,
                    bca = bca)

    ## calculation of MC estimates, SD and confidence intervals

    mc_est <- estMC(react_est_mtx,
                    ci_fun = ci_fun,
                    conf_level = 0.95)

    ## significance testing: H0 - fold regulation is mu = 0

    mc_p_vals <- simPval(react_est_mtx,
                         mu = 1,
                         median_estimate = mc_estimate == 'median')

    rownames(mc_est) <- colnames(react_est_mtx)
    names(mc_p_vals) <- names(react_est_mtx)

    if(mc_estimate == 'median') {

      mc_est <- mc_est[, c(1, 3:5)]

    } else {

      mc_est <- mc_est[, -1]

    }

    colnames(mc_est) <- c('fold_reg', 'error', 'lower_ci', 'upper_ci')

    mc_est <- rownames_to_column(as.data.frame(mc_est), 'react_id')

    mc_est[['p_value']] <- mc_p_vals

    mc_est[['p_adjusted']] <- p.adjust(mc_p_vals, 'BH')

    mc_est <- as_tibble(mc_est)

    if(!return_genes & !return_mc) {

      return(mc_est)

    } else {

      if(!return_mc) {

        return(list(gene_map = gene_tbl,
                    reg = mc_est))

      } else {

        return(list(gene_map = gene_tbl,
                    reg = mc_est,
                    mc = react_est_mtx))

      }

    }

  }

# END ------
