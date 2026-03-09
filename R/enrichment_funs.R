# functions for subsystem enrichment analyses

#' Subsystem enrichment analysis for differentially regulated reactions.
#'
#' @description
#' Function `suba()` performs enrichment analyses for member reactions of
#' metabolic subsystems among significantly activated and inhibited reactions in
#' an \code{\link{actiData}} object.
#' Please note, that prior to the analyses, significantly activated and
#' inhibited reactions need to by determined by calling
#' \code{\link{identify_regulated}}.
#'
#' @details
#' The analyses are done with \code{\link[fastTest]{f_enrichment}} with two methods:
#' Fisher's exact test (`type = "fisher`) or by random drawing from the entire
#' reaction pool (`type = "random"`).
#'
#' @return a data frame with enrichment analysis results:
#' numbers of reactions in the activated and inhibited reaction sets and in
#' the subsystems, enrichment magnitude measured by odds ratio (OR) statistic,
#' raw p values, and p values adjusted for multiple testing.
#'
#' @param x a \code{\link{actiData}} object.
#' @param type type of the hypothesis test: Fisher's exact test or random drawing
#' from the whole reaction pool.
#' @param ... additional arguments passed to \code{\link[fastTest]{f_enrichment}}
#' which specify e.g. type of confidence intervals for the OR statistic and
#' number of random draws, and serial/parallel execution of the testing algorithm.
#'
#' @seealso [get_regulation()]
#'
#' @export

  suba <- function(x,
                   type = c("fisher", "random"), ...) {

    ## the input control -------

    if(!is_actiData(x)) {

      stop("`x` has to be an `actiData` object created e.g. with `get_regulation()`",
           call. = FALSE)

    }

    type <- match.arg(type[[1]], c("fisher", "random"))

    if(!"regulation" %in% names(x[["reg"]])) {

      stop(paste("No regulation status in the `x` object.",
                 "Please call `identify_regulated()`",
                 "prior to launching `suba()`."),
           call. = FALSE)

    }

    ## the testing data ---------

    ### vectors of significantly regulated reactions

    signif_reactions <-
      filter(x$reg,
             .data[["regulation"]] %in% c("activated", "inhibited"))

    if(nrow(signif_reactions) == 0) {

      stop("No significantly regulated reactions.", call. = FALSE)

    }

    signif_reactions <-
      split(signif_reactions[["id"]],
            f = signif_reactions[["regulation"]])

    signif_reactions <-
      compact(map(signif_reactions,
                  function(x) if(length(x) == 0) NULL else x))

    ### dictionary of reactions in the subsystems
    ### vector with identifiers of all reactions

    reaction_dict <- split(x[["reg"]][["id"]],
                           x[["reg"]][["subsystem"]])

    all_reactions <- unique(x[["reg"]][["id"]])

    ## testing and customizing the results ---------

    result <- list()

    for(i in names(signif_reactions)) {

      result[[i]] <-  f_enrichment(x = signif_reactions[[i]],
                                   type = type,
                                   dict = reaction_dict,
                                   all = all_reactions,
                                   as_data_frame = TRUE,
                                   adj_method = "BH", ...)

    }

    reaction_set <- NULL

    result <-
      map2_dfr(result, names(result),
               ~mutate(.x,
                       reaction_set = factor(.y,
                                             c("activated", "inhibited"))))

    n_total_reaction_set <- NULL
    subsystem <- NULL
    n_total_subsystem <- NULL
    n_intersect <- NULL
    or <- NULL
    lower_ci <- NULL
    upper_ci <- NULL
    p_value <- NULL
    p_adjusted <- NULL

    transmute(as_tibble(result),
              reaction_set = .data[["reaction_set"]],
              n_total_reaction_set = .data[["n_x_total"]],
              subsystem = .data[["entry_name"]],
              n_total_subsystem = .data[["n_entry"]],
              n_intersect = .data[["n_intersect"]],
              or = .data[["or"]],
              lower_ci = if(type == "fisher") NULL else .data[["lower_ci"]],
              upper_ci = if(type == "fisher") NULL else .data[["upper_ci"]],
              p_value = .data[["p_value"]],
              p_adjusted = .data[["p_adjusted"]])

  }

# END ----------
