# S3 methods for the `actiData` class

# Appearance --------

#' Print method for `actiData` objects.
#'
#' @description
#' Print method for `actiData` objects.
#'
#' @return nothing, called for its side effects.
#'
#' @param x an object of class \code{\link{actiData}}.
#' @param ... additional arguments, currently none.
#'
#' @export

  print.actiData <- function(x, ...) {

    ## header contents

    n_reactions <- nrow(x[["reg"]])
    n_subsystems <- length(unique(na.omit(x[["reg"]][["subsystem"]])))

    cat(paste("`reactDB` object with activity estimates of", n_reactions,
              "in", n_subsystems, "subsystems"))
    cat("\n\n")
    print(x[["reg"]])

  }

# Accessing the object's components ---------

#' Access components of `actiData` objects.
#'
#' @description
#' Extracts the `reg` and `mc` components of the object.
#'
#' @return the output specified by `type` argument:
#' a data frame with reaction identifiers, names, subsystems and regulation
#' estimates (`type = "reg"`), a numeric matrix with reaction activity estimates
#' in Monte Carlo simulations (`type = "mc"`).
#'
#' @param object a \code{\link{actiData}} object.
#' @param type type of the output, see above.
#' @param ... extra arguments, currently none.
#'
#' @export

  components.actiData <- function(object, type = c("reg", "mc"), ...) {

    stopifnot(is_actiData(object))

    type <- match.arg(type[1], c("reg", "mc"))

    object[[type]]

  }

# Identification of significantly activated and inhibited reactions --------

#' Identity significantly regulated reactions.
#'
#' @description
#' Identifies significantly activated and inhibited reactions in
#' a \code{\link{actiData}} objects based on a cutoff of raw/adjusted p value
#' and activity estimate.
#'
#' @details
#' This method works only with \code{\link{actiData}} which store statistical
#' significance information (i.e. columns `p_value` and `p_adjusted`).
#' If this information is absent a warning is raised and regulation status
#' is inferred from the fold-regulation estimate only.
#'
#' The `fold_cutoff` arguments defines the fold-regulation cutoff used for
#' identification of significantly regulated reactions.
#' Reactions with regulation estimates lower than `1/fold_cutoff` and passing
#' the significance criterion are identified as inhibited.
#' Reactions with regulation estimates higher than `fold_cutoff` and passing
#' the significance criterion are identified as inhibited.
#'
#' @return
#' a \code{\link{actiData}} object.
#' The `reg` element of the object is appended with additional columns:
#' `p_cutoff` (cutoff of the p-value), `p_type` (type of the p value used to
#' identify significant effects), `fold_cutoff` (cutoff of the fold-regulation estimate),
#' `regulation` (regulation status of the reaction: activated, inhibited, or
#' non-significant/ns).
#'
#' @param x a \code{\link{actiData}} object.
#' @param p_type type of p-values to be used for definition of
#' regulated reactions: "raw" or "adjusted" (FDR-corrected, default).
#' @param p_cutoff cutoff of the p value used for definition of regulated
#' reactions, see Details.
#' @param fold_cutoff fold-regulation level cutoff used for definition
#' of regulated reactions.
#' @param ... extra arguments passed to methods.
#'
#' @seealso [get_regulation()]
#'
#' @export

  identify_regulated <- function(x, ...) UseMethod("identify_regulated")

#' @rdname identify_regulated
#' @export

  identify_regulated.actiData <- function(x,
                                          p_type = c("adjusted", "raw"),
                                          p_cutoff = 0.05,
                                          fold_cutoff = 1, ...) {

    ## input control --------

    stopifnot(is_actiData(x))

    signif_present <- TRUE

    if(any(!c("p_value", "p_adjusted") %in% names(x[["reg"]]))) {

      warning("No significance information in the object.", call. = FALSE)

      signif_present <- FALSE

    }

    p_type <- match.arg(p_type[1], c("adjusted", "raw"))

    p_err <- "`p_cutoff` must be a numeric value in (0, 1) range."

    if(!is.numeric(p_cutoff)) stop(p_err, call. = FALSE)
    p_cutoff <- p_cutoff[1]
    if(p_cutoff <= 0 | p_cutoff >= 1) stop(p_err, call. = FALSE)

    fold_err <- "`fold_cutoff` must be a numeric value >= 1."

    if(!is.numeric(fold_cutoff)) stop(fold_err, call. = FALSE)
    fold_cutoff <- fold_cutoff[1]
    if(fold_cutoff < 1) stop(fold_err, call. = FALSE)

    p_variable <- switch(p_type,
                         adjusted = "p_adjusted",
                         raw = "p_value")

    ## identification of the significant effects ------

    if(!signif_present) {

      p_type <- NA
      p_cutoff <- NA

    }

    x[["reg"]][["p_type"]] <- p_type
    x[["reg"]][["p_cutoff"]] <- p_cutoff
    x[["reg"]][["fold_cutoff"]] <- fold_cutoff

    x[["reg"]][["regulation"]] <-
      ifelse(x[["reg"]][["fold_reg"]] >= fold_cutoff,
             "activated",
             ifelse(x[["reg"]][["fold_reg"]] <= 1/fold_cutoff, "inhibited", "ns"))

    if(signif_present) {

      x[["reg"]][["regulation"]] <-
        ifelse(x[["reg"]][[p_variable]] >= p_cutoff,
               "ns",
               x[["reg"]][["regulation"]])

    }

    x[["reg"]][["regulation"]] <-
      factor(x[["reg"]][["regulation"]],
             c("activated", "inhibited", "ns"))

    return(x)

  }

# Transformation of regulation activity estimates --------

#' Transform of regulation activity estimates.
#'
#' @description
#' Transformation of estimates of regulation of reaction activity with a
#' function.
#'
#' @details
#' Function `transform_estimates` transforms estimates of reaction activity
#' along with their errors and confidence intervals (if present) with an
#' user-provided function (`log2` by default).
#' Note, no checks of validity of the function output are conducted!
#'
#' @return
#' a \code{\link{actiData}} object.
#' The `reg` element of the object is appended with additional with transformed
#' fold-regulation estimates, their errors, and bounds of the confidence intervals.
#' Names of these new column are preceded with a character string specified by
#' `prefix`
#'
#' @param x a \code{\link{actiData}} object.
#' @param fun transformation function.
#' @param prefix a character string preceding names of columns storing the
#' transformed estimates, errors, and confidence interval's bounds.
#' @param ... additional arguments passed to methods.
#'
#' @seealso [get_regulation()]
#'
#' @export

  transform_estimates <- function(x, ...)  UseMethod("transform_estimates")

#' @rdname transform_estimates
#' @export

  transform_estimates.actiData <- function(x, fun = log2, prefix = "log2_", ...) {

    ## input control --------

    stopifnot(is_actiData(x))

    if(!is_function(fun)) stop("`fun` has to be a function.", call. = FALSE)

    if(!is.character(prefix)) stop("`prefix` has to be a character string.", call. = FALSE)

    ## transformation ---------

    for(i in c("fold_reg", "error", "lower_ci", "upper_ci")) {

      if(!i %in% names(x[["reg"]])) next

      new_var_name <- paste0(prefix, i)

      x[["reg"]][[new_var_name]] <- fun(x[["reg"]][[i]])

    }

    return(x)

  }

# Counts of significantly regulated reactions -------

#' Count regulated reactions per subsystem.
#'
#' @description
#' Computes counts of significantly activated and inhibited reactions:
#' in the whole database and for particular subsystems.
#' Please make sure that the object contains activity regulation status
#' information.
#' If not, please call \code{\link{identify_regulated}} prior to
#' counting the significant effects.
#'
#' @details
#' If the object has no regulation status information (column `regulation` in the
#' `reg` element), an error is raised.
#'
#' @return a data frame with reaction counts according
#' to their regulation status.
#'
#' @param x a \code{\link{actiData}} object.
#' @param ... additional arguments passed to methods.
#'
#' @seealso [get_regulation()]
#'
#' @export

  count_regulated <- function(x, ...) UseMethod("count_regulated")

#' @rdname count_regulated
#' @export

  count_regulated.actiData <- function(x, ...) {

    ## input control --------

    stopifnot(is_actiData(x))

    if(!"regulation" %in% names(x[["reg"]])) {

      stop(paste("No regulation status information in the object.",
                 "Please call `identify_regulated()` prior to counting."),
           call. = FALSE)

    }

    ## counting of regulated reactions ---------

    reg <- x[["reg"]]

    reg[["subsystem"]] <-
      ifelse(!is.na(reg[["subsystem"]]),
             reg[["subsystem"]],
             "unknown")

    sub_reg <-
      c(list(all = reg),
        split(reg, reg[["subsystem"]]))

    n_numbers <- map(sub_reg, count, .data[["regulation"]], .drop = FALSE)

    n_total <- NULL
    percent <- NULL
    subsystem <- NULL
    regulation <- NULL

    n_numbers <- map(n_numbers,
                     mutate,
                     n_total = sum(.data[["n"]]),
                     percent = 100 * .data[["n"]]/sum(.data[["n"]]))

    n_numbers <- map2_dfr(n_numbers, names(n_numbers),
                          ~mutate(.x, subsystem = .y))

    relocate(filter(n_numbers,
                    .data[["regulation"]] != "ns",
                    !is.na(.data[["regulation"]])),
             regulation)

  }

# Sub-setting: selection of reactions and subsystems ---------

#' Select reactions or subsystems.
#'
#' @description
#' Select reactions or subsystems in an \code{\link{actiData}}.
#'
#' @return an \code{\link{actiData}} object with reactions and/or subsystems
#' of interest.
#'
#' @param .data a \code{\link{actiData}} object.
#' @param reactions a character vector with BiGG IDs of metabolic reactions.
#' @param subsystems a character vector with names of metabolic subsystems.
#' @param ... additional arguments passed to methods.
#'
#' @seealso [get_regulation()]
#'
#' @export select.actiData
#' @export

  select.actiData <- function(.data,
                              reactions = NULL,
                              subsystems = NULL, ...) {

    ## input control ---------

    stopifnot(is_actiData(.data))

    reacts_present <- unique(.data[["reg"]][["id"]])
    subs_present <- unique(.data[["reg"]][["subsystem"]])

    if(!is.null(reactions)) {

      stopifnot(is.character(reactions))

      reacts_missing <- setdiff(reactions, reacts_present)

      if(length(reacts_missing) > 0) {

        if(length(reacts_missing) > 20) {

          react_txt <-
            paste0(paste(reacts_missing[1:20], collapse = ", "), ", ...")

        } else {

          react_txt <- paste(reacts_missing, collapse = ", ")

        }

        warning(paste("Some reaction were not fount in the object:",
                      react_txt),
                call. = FALSE)

      }

    } else {

      reactions <- reacts_present

    }

    if(!is.null(subsystems)) {

      stopifnot(is.character(subsystems))

      subs_missing <- setdiff(subsystems, subs_present)

      if(length(subs_missing) > 0) {

        if(length(subs_missing) > 20) {

          sub_txt <-
            paste0(paste(subs_missing[1:20], collapse = ", "), ", ...")

        } else {

          sub_txt <- paste(subs_missing, collapse = ", ")

        }

        warning(paste("Some reaction were not fount in the object:",
                      sub_txt),
                call. = FALSE)

      }

    } else {

      subsystems <- subs_present

    }

    ## filtering -----------

    sub_reactions <- filter(.data[["reg"]],
                            .data[["subsystem"]] %in% subsystems)$id

    reactions <- intersect(sub_reactions, reactions)

    .data[["reg"]] <- filter(.data[["reg"]],
                             .data[["id"]] %in% reactions)

    if("mc" %in% names(.data)) .data[["mc"]] <- .data[["mc"]][, reactions]


    return(.data)

  }

# Plotting -----------

#' Diagnostic plots and visualization of reaction activity in `actiData` objects.
#'
#' @description
#' Diagnostic plots and visualization of estimates of regulation of reaction
#' activity - the plotting method for \code{\link{actiData}} objects.
#'
#' @details
#' The `plot()` method for \code{\link{actiData}} objects generates the following
#' types of plots specified by `plot_type` argument:
#'
#' * `plot_type = "errors"`: scatter plot of regulation estimates and their errors
#' made with \code{\link{plot_errors}}.
#'
#' * `plot_type = "numbers"`: bar/stack plots of numbers or percentages of
#' differentially regulated
#' reactions in subsystems generated with \code{\link{plot_numbers}}. Please note
#' that the regulation status of metabolic reactions in `x` has to be determined
#' prior to plotting by calling \code{\link{identify_regulated}}.
#'
#' * `plot_type == "mc"`: plots of reactions regulation estimates in single
#' iterations of the Monte Carlo simulation with \code{\link{plot_mc}}.
#'
#' @return a `ggplot` graphic object.
#'
#' @param x an \code{\link{actiData}} object.
#' @param plot_type plot type, see Details.
#' @param ... additional arguments passed to the plotting functions specified
#' in Details.
#'
#' @seealso [get_regulation()]
#'
#' @md
#'
#' @export plot.actiData
#' @export

  plot.actiData <- function(x,
                            plot_type = c("errors",
                                          "numbers",
                                          "mc"), ...) {

    ## input control --------

    stopifnot(is_actiData(x))

    plot_type <- match.arg(plot_type[1],
                           c("errors", "numbers", "mc"))

    ## re-routing ----------

    switch(plot_type,
           errors = plot_errors(x, ...),
           numbers = plot_numbers(x, ...),
           mc = plot_mc(x, ...))

  }

# END ---------
