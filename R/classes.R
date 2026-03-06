# S3 class definitions

# `reactDB` class --------

#' Database of reaction annotation: `reactDB` class.
#'
#' @description
#' Constructs an instance of a reaction annotation database data frame with
#' processed gene - reaction association rules and gene - reaction rule evaluation
#' expressions.
#'
#' @details
#' The function conducts basic validation for conformity with reaction
#' evaluation tools.
#' The input and output data frames have the following columns:
#' * __id__: reaction identifier beginning with `"R_"` string
#' * __name__: character strings with reaction names
#' * __subsystem__: character strings with names of Recon subsystems
#' * __gene_association__: character strings with gene - reaction association rules
#' * __entrez_id__: a column with lists of Entrez IDs of genes associated with
#' the reactions
#' * __exprs__: a column with `NULL`, R symbols or R language expressions used
#' to evaluate the gene - reaction association rules
#'
#' @return `reactDB()`: a data frame of class `reactDB` described in Details;
#' `is_reactDB()`: a logical value indicating it an object is an instance of `reactDB` class.
#'
#' @param x a data frame with columns specified in Details or an R object.
#' @param ... additional arguments, currently none.
#'
#' @md
#' @export

  reactDB <- function(x, ...) {

    ## input controls --------

    if(is_reactDB(x)) return(x)

    if(!is.data.frame(x)) stop("`x` has to be a data frame.", call. = FALSE)

    fix_cols <-
      c("id", "name", "subsystem", "gene_association", "entrez_id", "exprs" )

    missing_cols <- setdiff(fix_cols, names(x))

    if(length(missing_cols) > 0) {

      stop(paste("The following obligatory columns are missing from `x`:",
                 paste(missing_cols, collapse = ", ")),
           call. = FALSE)

    }

    char_cols <- c("id", "name", "subsystem", "gene_association")

    class_check <- map_lgl(x[char_cols], is.character)

    if(any(!class_check)) {

      stop(paste("The following obligatory column in `x` must be of character type:",
                 paste(char_cols[!class_check]), collapse = ", "),
           call. = FALSE)

    }

    if(!is.list(x[["entrez_id"]])) stop("Column `entrez_id` must be a list of Entrez IDs.", call. = FALSE)

    if(!is.list(x[["exprs"]])) stop("Column `exprs` must be a list.", call. = FALSE)

    class_check <-
      map_lgl(x[["exprs"]], function(x) is.call(x) | is.name(x) | is.null(x))

    if(any(!class_check)) {

      stop(paste("Unrecognized objects in `exprs` columns.",
                 "The allowed formats are NULL, R calls, or names."),
           call. = FALSE)

    }

    ## the structure -------

    structure(x, class = c("reactDB", class(x)))

  }

#' @rdname reactDB
#' @export

  is_reactDB <- function(x) inherits(x, "reactDB")

# Storage of reaction activity data `actiData` class ---------

#' Storage of estimates of regulation of activity of metabolic reactions: `actiData` class.
#'
#' @description
#' A list for storage of estimates of regulation of activity of metabolic reactions
#' and, optionally, results of Monte Carlo simulations of reaction activity.
#'
#' @details
#' `actiData` class objects are lists with the obligatory element `reg`,
#' a data frame which stores the reaction regulation estimates, and `mc`,
#' a numeric matrix of regulation estimates in single iterations of Monte Carlo
#' simulations.
#' The `reg` data frame is expected to have the following obligatory columns:
#'
#' * __id__ and __name__: character variables storing identifiers of metabolic reactions
#' * __subsystem__: a character variable with names of Recon subsystems
#' * __fold_reg__: a numeric variable with estimates of fold-regulation of reaction
#' activity.
#'
#' The optional columns in `reg` are:
#'
#' * __error__, __lower_ci__, and __upper_ci__: numeric variables with errors,
#' lower and upper bounds of confidence intervals of fold-regulation of reaction activity.
#' * __z__: a numeric variable with values of the test statistic Z
#' * __p_value__ and __p_adjusted__:  numeric variables with raw p values and
#' p-values adjusted for multiple testing.
#'
#'
#' @return `actiData()`: a list of class `actiData` as described in Details,
#' `is_actiData()`: a logical value indicating if an object is an instance of
#' `actiData` class.
#'
#' @param x an object
#' @param reg a data frame with estimates of regulation of metabolic reactions
#' as described in Details.
#' @param mc `NULL` or a numeric matrix with estimates of activity of metabolic
#' reactions obtained in Monte-Carlo simulations.
#' @param ... additional arguments, currently none.
#'
#' @md
#' @export

  actiData <- function(reg, mc = NULL, ...) {

    ## input controls: mandatory columns in `reg` ---------

    if(!is.data.frame(reg)) stop("`reg` has to be a data frame.", call. = FALSE)

    fix_cols <- c("id", "name", "subsystem", "fold_reg")

    missing_cols <- setdiff(fix_cols, names(reg))

    if(length(missing_cols) > 0) {

      stop(paste("The following obligatory columns are missing from `x`:",
                 paste(missing_cols, collapse = ", ")),
           call. = FALSE)

    }

    char_cols <- c("id", "name", "subsystem")

    char_check <- map_lgl(reg[, char_cols], is.character)

    if(any(!char_check)) {

      stop(paste("The following columns in `x` must be character variables:",
                 paste(char_cols[!char_check], collapse = ", ")),
           call. = FALSE)

    }

    if(!is.numeric(reg[["fold_reg"]])) {

      stop("`fold_reg` column in `x` has to be numeric.", call = FALSE)

    }

    ## input controls: optional columns in `reg` ------

    if("z" %in% names(reg)) stopifnot(is.numeric(reg[["z"]]))
    if("error" %in% names(reg)) stopifnot(is.numeric(reg[["error"]]))
    if("lower_ci" %in% names(reg)) stopifnot(is.numeric(reg[["lower_ci"]]))
    if("upper_ci" %in% names(reg)) stopifnot(is.numeric(reg[["upper_ci"]]))
    if("p_value" %in% names(reg)) stopifnot(is.numeric(reg[["p_value"]]))
    if("p_adjusted" %in% names(reg)) stopifnot(is.numeric(reg[["p_adjusted"]]))

    ## input controls: `mc` matrix --------

    if(!is.null(mc)) {

      mc_err <- "`mc` has to be a numeric matrix."

      if(!is.matrix(mc)) stop(mc_err, call. = FALSE)
      if(!is.numeric(mc)) stop(mc_err, call. = FALSE)

      reaction_ids <- reg[["id"]]
      mc_ids <- colnames(mc)

      stopifnot(!is.null(mc_ids))

      if(any(!reaction_ids %in% mc_ids) | any(!mc_ids %in% reaction_ids)) {

        stop("Reactions in `reg` and `mc` do not match.", call. = FALSE)

      }

    }

    ## the structure --------

    if(is.null(mc)) {

      return(structure(list(reg = reg),
                       class = "actiData"))

    }

    return(structure(list(reg = reg,
                          mc = mc),
                     class = "actiData"))

  }

#' @rdname actiData
#' @export

  is_actiData <- function(x) inherits(x, "actiData")

# END ------
