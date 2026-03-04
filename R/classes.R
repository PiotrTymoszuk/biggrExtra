# S3 class definitions

# geneSBML class ------

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

# END ------
