# S3 methods for the `reactDB` class

# appearance ---------

#' Appearance for `reactDB` objects.
#'
#' @description
#' `print()` method for `reactDB` objects.
#'
#' @return nothing, called for its side effects.
#'
#' @param x an \code{\link{reactDB}} object.
#' @param ... additional arguments, currently none.
#'
#' @export

  print.reactDB <- function(x, ...) {

    stopifnot(is_reactDB(x))

    ## header of the printed text

    n_reactions <- nrow(x)
    n_rules <- nrow(filter(x, !is.na(.data[["gene_association"]])))

    cat(paste("`reactDB` object with", n_reactions,
              "reactions and", n_rules, "gene - reaction association rules"))
    cat("\n\n")
    NextMethod("print", x)

  }

# summary ----------

#' Summary method for `reactDB` objects.
#'
#' @description
#' `summary()` method for `reactDB` objects.
#'
#' @return
#' A data frame with numbers of reactions, subsystems, genes, and gene - reaction
#' association rules.
#'
#' @param object an \code{\link{reactDB}} object.
#' @param ... additional arguments, currently none.
#'
#' @export

  summary.reactDB <- function(object, ...) {

    stopifnot(is_reactDB(object))

    n_numbers <-
      map_dbl(object[, c("id", "subsystem", "gene_association")],
              ~length(unique(na.omit(.x))))

    n_numbers <- c(n_numbers,
                   c("genes" = length(unique(unlist(object[["entrez_id"]])))))

    feature <- NULL
    n <- NULL

    tibble(feature = c("subsystems",
                       "reactions",
                       "genes",
                       "association rules"),
           n = n_numbers[c("subsystem",
                           "id",
                           "genes",
                           "gene_association")])

  }

# END --------
