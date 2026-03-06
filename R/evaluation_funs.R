# Functions used for evaluation of gene - reaction association rules

# Evaluation of the gene mapping rules ----

#' Evaluate gene mapping rules.
#'
#' @description
#' Evaluates a list of gene mapping rules stored as expression.
#' Such rules are retrieved from data frames of reaction annotations
#' with \code{\link{extract_genes}} or \code{\link{as_reactDB}}.
#'
#' @return a numeric vector with the reaction regulation estimates.
#' Named with reaction IDs.
#'
#' @param x a vector with gene regulation estimates named
#' with Entrez ID identifiers.
#' @param database an object of class \code{\link{reactDB}} created with
#' functions \code{\link{extract_genes}} or \code{\link{as_reactDB}}.
#' @param parent_env an environments containing evaluation functions and values
#' for genes missing from the x gene regulation estimate vector.

  eval_rules <- function(x, database, parent_env) {

    ## entry controls ------

    stopifnot(is.numeric(x))
    stopifnot(!is.null(names(x)))

    stopifnot(is_reactDB(database))

    stopifnot(is_environment(parent_env))

    ## serial evaluation --------

    exprs_lst <- set_names(database[["exprs"]],
                           database[["id"]])

    eval_env <-
      new_environment(data = as.list(x),
                      parent = parent_env)

    map_dbl(compact(exprs_lst), eval, envir = eval_env)

  }

# Calculation of the regulation estimates with differential gene expression data -------

#' Calculate reaction activity estimates based on differential gene expression.
#'
#' @description
#' Intended for internal use.
#' Function `get_regulation()` computes fold-regulation estimates of metabolic
#' reaction activity given estimates of differential gene-expression with their
#' errors and gene - reaction association rules in a metabolic reaction annotation
#' database object..
#' If the errors are provided, statistical significance and confidence intervals
#' of the reaction activity are computed with two methods: (1) derived from normal
#' distribution of the error sums or (2) computed by serial evaluation of the
#' gene - reaction association rules for random draws from normal distributions
#' of the estimates of differential gene expression.
#' In the later case, the differential gene expression estimate serves as the
#' distribution mean and the error of differential gene expression serves as
#' the standard deviation.
#'
#' @return
#' The function returns a list with the following elements:
#'
#' * __reg__ : a data frame with reaction identifiers (`id`),
#' fold-regulation estimates of reaction activity (`fold_reg`), estimates of
#' error and 95% confidence intervals of reaction activity (only if `err` vector was
#' provided, columns: `error`, `lower_ci`, `upper_ci`), testing statistic Z (`z`, only
#' if `err` was provided and `err_method = "norm"`), raw p values and p values
#' adjusted for multiple testing with the Benjamini-Hochberg method
#' (only if `err` vector was provided, columns: `p_value` and `p_adjusted`).
#'
#' * __mc__: returned only if `err_method = "mc"` and `return_mc = TRUE`, a numeric
#' matrix with estimates of fold-regulation of reaction activity in iterations of
#' the simulation algorithm. Reactions are in columns, iterations are in rows.
#'
#' @inheritParams eval_rules
#' @param err a vector with errors gene regulation estimates named
#' with Entrez ID identifiers.
#' @param scale specifies the scale of the provided expression regulation and
#' error estimates: "identity" (default), "log2", or "log".
#' @param or_fun one of "mean" (default), "median", "min", or "max".
#' Name of the function used to handle the `OR` operator between the genes in
#' the association rule.
#' @param and_fun one of "mean", "median", "min" (default), or "max".
#' Name of the function used to handle the `OR` operator between the genes in
#' the association rule.
#' @param x_default the default numeric value or `NA` for regulation of the genes
#' present in the database and absent from the regulation data. If `NA`,
#' the genes missing from `x` won't be used in calculation of the reaction
#' activity estimates.
#' @param err_method specifies the method of calculation of errors for reactions
#' as described in Details: "norm" (from normal distribution, default) or
#' "mc" (Monte Carlo simulation).
#' @param return_mc logical, should all Monte Carlo-determined regulation
#' estimates for reactions be returned? Defaults to `FALSE`.
#' @param n_iter number of iterations for the Monte Carlo simulation of
#' the reaction errors.
#' @param mc_estimate statistic used to compute the reaction activity estimate
#' in Monte Carlo simulation: either mean or median over algorithm iterations.
#' Ignored if no errors were provided or `err_method = "norm"`.
#' @param burn_in number of Monte Carlo iterations to be discarded prior to
#' computation of activity estimates and activity regulation errors. Ignored
#' if no errors were provided or `err_method = "norm"`.
#' @param ci_method method of calculation of confidence intervals:
#' BCA ("bca", default) or percentile ("perc").
#' @param seed a seed for the random number generator.
#' @param .parallel logical, should th computation be run in parallel?
#' @param ... additional argument, currently none.
#'
#' @md
#' @export

  get_regulation <- function(x,
                             err = NULL,
                             database,
                             scale = c("identity", "log2", "log"),
                             or_fun = c("mean", "median", "min", "max"),
                             and_fun = c("min", "max", "mean", "median"),
                             x_default = 1,
                             err_method = c("norm", "mc"),
                             return_mc = FALSE,
                             n_iter = 1000,
                             mc_estimate = c("mean", "median"),
                             burn_in = 0,
                             ci_method = c("bca", "perc"),
                             seed = NULL,
                             .parallel = TRUE, ...) {

    ## benchmarking and exit handlers --------

    start_time <- Sys.time()

    on.exit(plan("sequential"))
    on.exit(message(paste("Elapsed:", Sys.time() - start_time)), add = TRUE)

    ## variable symbol definitions --------

    id <- NULL
    fold_reg <- NULL
    error <- NULL
    z <- NULL
    p_value <- NULL

    ## entry control --------

    x_err <- "`x` has to be a numeric vector with  Entrez IDs as names."

    if(!is.numeric(x)) stop(x_err, call. = FALSE)
    if(is.null(names(x))) stop(x_err, call. = FALSE)

    if(!is.null(err)) {

      err_err <- "`err` has to be a numeric vector with  Entrez IDs as names."

      if(!is.numeric(err)) stop(err_err, call. = FALSE)
      if(is.null(names(err))) stop(err_err, call. = FALSE)

    }

    if(!is_reactDB(database)) {

      stop(paste("`database` has to be a `reactDB` object.",
                 "Please consult functions `as_reactDB()` and `extract_genes()`."),
           call. = FALSE)

    }

    scale <- match.arg(scale[1], c("identity", "log2", "log"))

    scale_fun <-
      switch(scale,
             identity = identity,
             log2 = function(x) 2^x,
             log = exp)

    identity_x <- x

    x <- scale_fun(x)

    if(scale == "identity") {

      if(any(x < 0)) {

        stop(paste("All elements of `x` have to be >= 0. Are you sure,",
                   "your values are on the identity scale?"),
             call. = FALSE)

      }

      if(!is.null(err)) {

        if(any(err < 0)) {

          stop(paste("All elements of `err` have to be >= 0.",
                     "Are you sure, your values are on the identity scale?"),
               call. = FALSE)

        }

      }

    }

    or_fun <- match.arg(or_fun[1], c("mean", "median", "min", "max"))
    and_fun <- match.arg(and_fun[1], c("min", "max", "mean", "median"))
    err_method <- match.arg(err_method[1], c("norm", "mc"))
    ci_method <- match.arg(ci_method[1], c("bca", "perc"))

    stopifnot(is.numeric(n_iter))
    n_iter <- as.integer(n_iter)

    mc_estimate <- match.arg(mc_estimate[1], c("mean", "median"))

    stopifnot(is.numeric(burn_in))
    burn_in <- as.integer(burn_in[1])

    stopifnot(is.logical(return_mc))
    return_mc <- return_mc[1]

    if(!is.na(x_default)) {

      if(!is.numeric(x_default)) stop("`x_default` has to be a number or NA.", call. = FALSE)

      x_default <- x_default[1]

    }

    ## building a parent evaluation environment -------

    message("Building an evaluation environment")

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

    all_genes <- reduce(database$entrez_id, union)

    miss_genes <- all_genes[!all_genes %in% names(x)]

    miss_gene_lst <- set_names(rep(x_default, length(miss_genes)),
                               miss_genes)
    fun_ev <-
      new_environment(data = c(list(`%AND%` = a_fun,
                                    `%OR%` = o_fun,
                                    `(` = `(`),
                               as.list(miss_gene_lst)))

    ## rule evaluation: determination of regulation estimates ------

    if(is.null(err) | err_method == "norm") {

      message("Calculation of reaction regulation estimates")

      reg_lst <- eval_rules(x = x,
                            database = database,
                            parent_env = fun_ev)

      reg_tbl <- tibble(id = names(reg_lst),
                        fold_reg = unlist(reg_lst))

      reg_tbl <- filter(reg_tbl,
                        !is.na(fold_reg),
                        !is.nan(fold_reg),
                        !is.infinite(fold_reg))

      ## filling with the default value if a reaction is absent
      ## from the regulation estimate table.

      all_reactions <- database[["id"]]

      miss_reactions <-
        setdiff(all_reactions, reg_tbl[["id"]])

      reg_tbl <- rbind(reg_tbl,
                       tibble(id = miss_reactions,
                              fold_reg = x_default))

      if(is.null(err)) return(list(reg = reg_tbl))

    }

    ## computing the sum error from the normal distribution -------

    if(err_method == "norm") {

      message("Calulation of the errors: normal distribution")

      err <- scale_fun(err)

      err_lst <- map(set_names(database[["entrez_id"]],
                               database[["id"]]),
                     ~err[.x])

      err_lst <- map(err_lst, sum, na.rm = TRUE)

      err_tbl <- tibble(id = names(err_lst),
                        error = unlist(err_lst))

      err_tbl <- filter(err_tbl, !is.na(error))

      reg_tbl <- left_join(reg_tbl, err_tbl, by = "id")

      reg_tbl <-
        mutate(reg_tbl,
               z = (fold_reg - 1)/error,
               z = ifelse(is.infinite(z), NA, z),
               lower_ci = fold_reg + error * qnorm(0.025),
               upper_ci = fold_reg + error * qnorm(0.975),
               p_value = pnorm(abs(z), lower.tail = FALSE) * 2,
               p_adjusted = p.adjust(p_value, "BH"))

      reg_tbl <-
        reg_tbl[c("id", "fold_reg", "error",
                  "lower_ci", "upper_ci",
                  "z", "p_value", "p_adjusted")]

      return(list(reg = reg_tbl))

    }

    ## computing the regulation estimates and errors: Monte Carlo -------

    message("Simulations: regulation estimate and error calculation")

    ### CI computing functions
    ### simulation of the gene regulation values

    sim_est <- draw_norm(identity_x,
                         err = err,
                         n = n_iter,
                         scale = scale,
                         seed = seed)

    iters <- rownames(sim_est)

    ## calculation of the estimates and errors

    if(.parallel) plan("multisession")

    exp_globals <- c("eval_rules",
                     "fun_ev",
                     "database")

    exp_packages <- c("generics",
                      "rlang",
                      "biggrExtra")

    react_est_lst <-
      future_map(iters,
                 ~eval_rules(x = sim_est[.x, ],
                             database = database,
                             parent_env = fun_ev),
                 .options = furrr_options(seed = TRUE,
                                          globals = exp_globals,
                                          packages = exp_packages))


    react_est_lst <- map(react_est_lst, ~ifelse(.x < 0, 0, .x))

    react_est_mtx <- do.call("rbind", react_est_lst)

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
                         median_estimate = mc_estimate == "median")

    rownames(mc_est) <- colnames(react_est_mtx)
    names(mc_p_vals) <- names(react_est_mtx)

    if(mc_estimate == "median") {

      mc_est <- mc_est[, c(1, 3:5)]

    } else {

      mc_est <- mc_est[, -1]

    }

    colnames(mc_est) <- c("fold_reg", "error", "lower_ci", "upper_ci")

    mc_est <- rownames_to_column(as.data.frame(mc_est), "id")

    mc_est[["p_value"]] <- mc_p_vals

    mc_est[["p_adjusted"]] <- p.adjust(mc_p_vals, "BH")

    mc_est <- as_tibble(mc_est)

    if(!return_mc) {

      return(list(reg = mc_est))

    } else {

      return(list(reg = mc_est,
                  mc = react_est_mtx))

    }

  }

# Reaction activity scores with bare expression data --------

#' Calculate reaction activity scores based on gene expression data.
#'
#' @description
#' Function `get_activity()` calculates activity scores of metabolic reactions
#' by evaluation of gene - reaction association rules in a reaction annotation
#' database given gene expression data.
#'
#' @details
#' Function `get_activity()` takes a numeric matrix or a numeric data frame `x`
#' with gene expression metrics to calculate activity scores of metabolic reactions.
#' The scores are calculated by simple evaluation of gene - reaction association
#' rules stored in a reaction annotation database data frame of class
#' \code{\link{reactDB}} (created with \code{\link{as_reactDB}} or
#' \code{\link{extract_genes}}).
#'
#' Of note, the gene expression metrics in `x` are not pre-processed prior to
#' calculation of the reaction activity scores; the pre-processing steps are
#' intended to be done by the user.
#' Because evaluation of the gene - reaction association rules uses simple
#' arithmetic with means, medians, minima or maxima of gene expression values,
#' it's highly recommended to bring the expression values in `x` approximate
#' equal scales, e.g. by computing Z-scores with \code{\link[microViz]{zScores}}.
#'
#' @return a numeric matrix or a data frame (if `as_data_frame = TRUE`) with
#' reaction activity scores.
#' Metabolic reactions named with their identifiers are in columns,
#' samples are in rows.
#' In the data frame output, sample identifiers are stored in the first column
#' named `sample_id`.
#'
#' @inheritParams get_regulation
#' @param x a numeric matrix or a data frame with gene expression metrics; genes
#' in columns, samples in rows.
#' @param as_data_frame logical, should the matrix of reaction activity scores
#' be coerced to a data frame?
#' @param id_col `NULL` or a character string specifying a variable serving as
#' a sample identifier. If `id_col = NULL`, the sample identifiers will be
#' extracted from row names of `x`.
#' @param variables `NULL` or a character vector specifying gene expression
#' variables used for calculation of reaction activity scores. If `NULL`, all
#' variables in `x` except of the sample identifier will be used.
#' @param ... additional arguments, currently none.
#'
#' @export

  get_activity <- function(x, ...) UseMethod("get_activity")

#' @rdname get_activity
#' @export

  get_activity.matrix <- function(x,
                                  database,
                                  or_fun = c("mean", "median", "min", "max"),
                                  and_fun = c("min", "max", "mean", "median"),
                                  x_default = 1,
                                  as_data_frame = TRUE,
                                  .parallel = FALSE, ...) {

    ## benchmarking and exit handlers --------

    start_time <- Sys.time()

    on.exit(plan("sequential"))
    on.exit(message(paste("Elapsed:", Sys.time() - start_time)), add = TRUE)

    ## entry control --------

    stopifnot(is.matrix(x))

    if(!is.numeric(x)) stop("`x` must be a numeric matrix.", call. = FALSE)

    if(is.null(colnames(x))) stop("Column names of `x` must be specified.", call. = FALSE)
    if(is.null(rownames(x))) stop("Row names of `x` must be specified.", call. = FALSE)

    if(!is_reactDB(database)) {

      stop(paste("`database` has to be a `reactDB` object.",
                 "Please consult functions `as_reactDB()` and `extract_genes()`."),
           call. = FALSE)

    }

    or_fun <- match.arg(or_fun[1], c("mean", "median", "min", "max"))
    and_fun <- match.arg(and_fun[1], c("min", "max", "mean", "median"))

    if(!is.na(x_default)) {

      if(!is.numeric(x_default)) stop("`x_default` has to be a number or NA.", call. = FALSE)

      x_default <- x_default[1]

    }

    stopifnot(is.logical(as_data_frame))
    as_data_frame <- as_data_frame[1]

    ## evaluation environment ---------

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

    all_genes <- reduce(database$entrez_id, union)

    miss_genes <- all_genes[!all_genes %in% colnames(x)]

    miss_gene_lst <- set_names(rep(x_default, length(miss_genes)),
                               miss_genes)
    fun_ev <-
      new_environment(data = c(list(`%AND%` = a_fun,
                                    `%OR%` = o_fun,
                                    `(` = `(`),
                               as.list(miss_gene_lst)))

    ## calculation of the reaction activity scores ---------

    row_data <- map(1:nrow(x), function(idx) x[idx, ])

    row_data <- set_names(row_data,
                          rownames(x))

    if(.parallel) plan("multisession")

    exp_globals <- c("eval_rules",
                     "fun_ev",
                     "database")

    exp_packages <- c("generics",
                      "rlang",
                      "biggrExtra")

    activity_mtx <-
      future_map(row_data,
                 eval_rules,
                 database = database,
                 parent_env = fun_ev,
                 .options = furrr_options(seed = TRUE,
                                          globals = exp_globals,
                                          packages = exp_packages))

    activity_mtx <- reduce(activity_mtx, rbind)

    rownames(activity_mtx) <- names(row_data)

    if(!as_data_frame) return(activity_mtx)

    return(as_tibble(rownames_to_column(as.data.frame(activity_mtx), "sample_id")))

  }

#' @rdname get_activity
#' @export

  get_activity.data.frame <- function(x,
                                      database,
                                      id_col = NULL,
                                      variables = NULL,
                                      or_fun = c("mean", "median", "min", "max"),
                                      and_fun = c("min", "max", "mean", "median"),
                                      x_default = 1,
                                      as_data_frame = TRUE,
                                      .parallel = FALSE, ...) {

    ## minimal input control --------

    stopifnot(is.data.frame(x))

    if(!is.null(id_col)) {

      stopifnot(is.character(id_col))
      id_col <- id_col[1]

      if(!id_col %in% names(x)) {

        stop("Variable specified by `id_col` is missing from `x`",
             call. = FALSE)

      }

      x <- column_to_rownames(x, id_col)

    }

    if(is.null(variables)) variables <- names(x)

    if(!is.character(variables)) {

      stop("`variables` has to be NULL or a character vector.", call. = FALSE)

    }

    missing_cols <- setdiff(variables, names(x))

    if(length(missing_cols) > 0) {

      if(length(missing_cols) > 20) {

        missing_txt <- paste0(paste(missing_cols[1:20], collapse = ", "), ", ...")

      } else {

        missing_txt <- paste(missing_cols, collapse = ", ")

      }

      stop(paste("Some variables specified by the user are missing from `x`:",
                 missing_txt),
           call. = FALSE)

    }

    x <- x[, variables, drop = FALSE]

    class_check <- map_lgl(x, is.numeric)

    if(any(!class_check)) {

      error_vars <- names(x)[!class_check]

      if(length(error_vars) > 20) {

        err_var_txt <- paste0(paste(error_vars[1:20], collapse = ", "), ", ...")

      } else {

        err_var_txt <- paste(error_vars, collapse = ", ")

      }

      stop(paste("Some variables in `x` are not nomeric:", err_var_txt), call. = FALSE)

    }

    ## calculation of the activity scores ---------

    get_activity(as.matrix(x),
                 database = database,
                 or_fun = or_fun,
                 and_fun = and_fun,
                 x_default = x_default,
                 as_data_frame = as_data_frame,
                 .parallel = .parallel, ...)

  }

# END --------
