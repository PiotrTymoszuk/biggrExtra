# S3 class definitions

# `geneSBML` class ------

#' Create a `geneSBML` object.
#'
#' @description
#' Creates a `geneSBML` object given a SBML model object and a
#' data frame with reaction regulation estimates and, optionally, regulation
#' errors, and a data frame with gene - reaction mapping.
#'
#' @details
#' In `model = NULL` and `mc = NULL`, a mamory saver object is returned
#' (subclass `memoSaver`).#'
#'
#' @return an instance of `geneSBML` class.
#'
#' @param x an object.
#' @param model `NULL` or a `SBML` model.
#' @param reg a data frame with reaction IDs (`react_id`), regulation estimates
#' (`fold_reg`) and, optionally, regulation errors (`error`).
#' @param gene_map a data frame with reaction IDs (`react_id`),
#' list of Entrez ID identifiers (`react_id`)
#' and evaluation expressions (`react_id`).
#' @param mc optional, a matrix with reaction regulation estimates in each step
#' of Monte Carlo simulation. Columns are reactions, rows are subsequent runs.
#'
#' @export

  geneSBML <- function(model,
                       reg,
                       gene_map,
                       mc = NULL) {

    ## entry control -------

    if(!is.null(model)) {

      if(!inherits(model, "Model")) {

        stop("model has to be a valid SBML model.", call. = FALSE)

      }

    }

    reg_error_txt <-
      "`reg` has to be a data frame with `react_id` and `fold_reg` columns."

    if(!is.data.frame(reg)) stop(reg_error_txt, call. = FALSE)

    if(any(!c("react_id", "fold_reg") %in% names(reg))) stop(reg_error_txt, call. = FALSE)

    map_error_txt <-
      paste("`gene_map` has to be a data frame with `react_id`,",
            "`entrez_id` and `exprs` columns.")

    if(!is.data.frame(gene_map)) stop(map_error_txt)

    if(any(!c("react_id", "entrez_id", "exprs") %in% names(gene_map))) {

      stop(map_error_txt, call. = FALSE)

    }

    if(!is.null(mc)) {

      mc_error_txt <- "mc has to be a numeric matrix."

      if(!is.numeric(mc)) stop(mc_error_txt, call. = FALSE)

      if(!is.matrix(mc)) stop(mc_error_txt, call. = FALSE)

    }

    ## the output object ---------

    out_lst <- list(model = model,
                    reg = reg,
                    gene_map = gene_map,
                    mc = mc)

    if(is.null(mc) & is.null(model)) {

      out_lst <-
        structure(out_lst, class = c("memoSaver", "geneSBML"))

    } else {

      out_lst <-
        structure(out_lst, class = "geneSBML")

    }

    out_lst

  }

#' @rdname geneSBML
#' @export

  is_geneSBML <- function(x) inherits(x, "geneSBML")

#' @rdname geneSBML
#' @export

  is_memoSaver <- function(x) inherits(x, "geneSBML") & inherits(x, "memoSaver")

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
#' * __id__: reaction identifier beginning with `"R_"` string.
#' * __name__: character strings with reaction names.
#' * __subsystem__: character strings with names of Recon subsystems.
#' * __gene_association__: character strings with gene - reaction association rules.
#' * __entrez_id__: a column with lists of Entrez IDs of genes associated with
#' the reactions.
#' * __exprs__: a column with `NULL`, R symbols or R language expressions used
#' to evaluate the gene - reaction association rules.
#'
#' @return a data frame of class `reactDB`.
#'
#' @param x a data frame with columns specified in Details.
#'
#' @md
#' @export

  reactDB <- function(x) {

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


# END ------
