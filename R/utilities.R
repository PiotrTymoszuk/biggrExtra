# Utilities

# Escape numerics in a string -----

#' Escape Entrez ID numbers in a string.
#'
#' @description
#' Finds Entrez ID in a string and escapes them with '``.
#'
#' @return a sting.
#'
#' @param str a string.

  escape_numbers <- function(str) {

    stopifnot(is.character(str))

    numbs <- sort(unique(unlist(stri_extract_all(str, regex = "\\d+"))))

    for(i in numbs) {

      search_rex <- paste0("(?<!\\d{1})", i, "(?!\\d+)")

      str <-
        stri_replace_all(str,
                         regex = search_rex,
                         replacement = paste0(" `", i, "` "))

    }

    return(str)

  }

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
      new_environment(data = as.list(x), parent = parent_env)

    map_dbl(exprs_lst, eval, envir = eval_env)

  }

# Simulation of gene regulation estimate distribution -------

#' Draw values from the gene regulation estimate distribution.
#'
#' @description
#' Draws values from normal distributions of gene regulation estimates.
#' The gene expression regulation estimate serves as the mean of the distribution,
#' the gene expression regulation error serves as the standard deviation of
#' the distribution.
#'
#' @return a numeric matrix: genes are columns, their regulation estimates are
#  rows.
#'
#' @param x a numeric vector with fold-regulation estimates.
#' Elements are named with Entrez IDs.
#' @param err a numeric vector with an error estimate such as SD or SEM,
#' named with Entrez IDs.
#' @param n number of values to be drawn.
#' @param scale name of the scale of the provided regulation estimates and
#' their errors: "identity" (default), "log2", or "log".
#' @param seed seed of the random number generator.

  draw_norm <- function(x,
                        err,
                        n,
                        scale = c("identity", "log2", "log"),
                        seed = NULL) {

    ## entry control -----

    stopifnot(is.numeric(x))
    stopifnot(!is.null(names(x)))

    stopifnot(is.numeric(err))
    stopifnot((!is.null(names(err))))

    stopifnot(is.numeric(n))
    n <- as.integer(n[1])
    stopifnot(n > 1)

    cmm_genes <- intersect(names(x), names(err))

    if(length(cmm_genes) == 0) {

      stop('Names of `x` and `err` do not fit. Right vectors?',
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

# END ------
