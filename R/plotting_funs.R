# functions for plotting

# Plotting of errors and regulation estimates ---------

#' Scatter plot of errors and regulation estimates for metabolic reactions.
#'
#' @description
#' Function `plot_errors()` takes an \code{\link{actiData}} object with estimates
#' of regulation of activity of metabolic reactions and visualizes the reaction
#' activity estimates and their errors in a scatter plot.
#'
#' @details
#' The function works only for \code{\link{actiData}} objects with errors of
#' the regulation estimates obtained from the normal distribution or by Monte
#' Carlo simulations, i.e. objects generated with \code{\link{get_regulation}}
#' with the non-`NULL` `err` argument.
#'
#' @return a `ggplot` scatter plot.
#'
#' @seealso [get_regulation()]
#'
#' @param x an \code{\link{actiData}} object.
#' @param reactions a character vector with BiGG IDs of metabolic reactions
#' to be presented in the plot. If `NULL`, all reactions will be displayed.
#' @param subsystems a character vector with names of metabolic subsystems to
#' be presented in the plot. If `NULL`, all subsystems will be displayed.
#' @param fun a function used to transform the regulation estimates ad their
#' errors prior to plotting. Defaults to `identity`.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param x_lab X axis label.
#' @param y_lab Y axis label.
#' @param show_trend logical, should a trend line be plotted
#' (with \code{\link[ggplot2]{geom_smooth}}).
#' @param line_color color of the trend line.
#' @param line_alpha alpha/opacity of the trend line.
#' @param show_vline logical: should a vertical line representing non-regulated
#' reactions (`fold_reg = 1`) be displayed in the plot?
#' @param vline_color color of the line representing no regulation.
#' @param vline_type type of the line representing no regulation.
#' @param vline_alpha alpha/opacity of the line representing no regulation.
#' @param show_hline logical: should a vertical line representing zero errors
#' (`error = 0`) be displayed in the plot?
#' @param hline_color color of the line representing zero errors.
#' @param hline_type type of the line representing zero errors.
#' @param hline_alpha alpha/opacity of the line representing zero errors.
#' @param point_color fill color of data points.
#' @param point_alpha alpha/opacity of data points.
#' @param point_size size of data points.
#' @param cust_theme `NULL` or a custom `ggplot` theme.
#' @param ... additional arguments passed to \code{\link[ggplot2]{geom_smooth}}.
#'
#' @export

  plot_errors <- function(x,
                          reactions = NULL,
                          subsystems = NULL,
                          fun = identity,
                          point_color = "steelblue",
                          point_size = 2,
                          point_alpha = 0.75,
                          show_trend = TRUE,
                          line_color = "black",
                          line_alpha = 1,
                          show_vline = TRUE,
                          vline_color = "black",
                          vline_type = "dashed",
                          vline_alpha = 1,
                          show_hline = TRUE,
                          hline_color = "black",
                          hline_type = "dashed",
                          hline_alpha = 1,
                          cust_theme = microViz::theme_micro(),
                          plot_title = NULL,
                          plot_subtitle = NULL,
                          x_lab = "regulation estimate",
                          y_lab = "error", ...) {

    ## input control ------

    if(!is_actiData(x)) stop("`x` has to be an `actiData` object.", call. = FALSE)

    if(!"error" %in% names(x$reg)) {

      warning("No errors in the object.", call. = FALSE)

      return(NULL)

    }

    if(!is_function(fun)) stop("`fun` has to be a function.", call. = FALSE)

    stopifnot(is.numeric(point_size))
    point_size <- point_size[1]

    stopifnot(is.numeric(point_alpha))
    point_alpha <- point_alpha[1]

    stopifnot(is.logical(show_trend))
    show_trend <- show_trend[1]

    stopifnot(is.numeric(line_alpha))
    line_alpha <- line_alpha[1]

    stopifnot(is.logical(show_vline))
    show_vline <- show_vline[1]

    stopifnot(is.numeric(vline_alpha))
    vline_alpha <- vline_alpha[1]

    stopifnot(is.logical(show_hline))
    show_hline <- show_hline[1]

    stopifnot(is.numeric(hline_alpha))
    hline_alpha <- hline_alpha[1]

    if(!is.null(cust_theme)) {

      if(!is.theme(cust_theme)) {

        stop("`cust_theme` has to be a `ggplot` theme object.", call. = FALSE)

      }

    }

    ## plotting ----------

    x <- select(x, reactions, subsystems)

    x <- transform_estimates(x, fun = fun, prefix = "transf_")

    err_plot <- ggplot(x$reg,
                       aes(x = .data[["transf_fold_reg"]],
                           y = .data[["transf_error"]]))

    if(show_vline) {

      err_plot <- err_plot +
        geom_vline(xintercept = fun(1),
                   color = vline_color,
                   alpha = vline_alpha,
                   linetype = vline_type)

    }

    if(show_hline) {

      err_plot <- err_plot +
        geom_hline(yintercept = fun(0),
                   color = hline_color,
                   alpha = hline_alpha,
                   linetype = hline_type)

    }

    err_plot <- err_plot +
      geom_point(shape = 21,
                 color = "black",
                 fill = point_color,
                 size = point_size,
                 alpha = point_alpha) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    if(!is.null(cust_theme)) err_plot <- err_plot + cust_theme

    if(show_trend) {

      err_plot <- err_plot +
        geom_smooth(color = line_color,
                    alpha = line_alpha, ...)

    }

    return(err_plot)

  }

# Plotting regulation estimates in Monte Carlo simulations --------

#' Plot reaction regulation estimates in iterations of Monte-Carlo simulations.
#'
#' @description
#' Function `plot_mc()` creates line plots of reaction regulation estimates in
#' iterations of the simulation algorithm, i.e. subsequent draws from normal
#' distribution of estimates of differential gene expression.
#'
#' @details
#' The function works only for \code{\link{actiData}} objects with Monte Carlo
#' estimates (i.e. created with \code{\link{get_regulation}} with
#' non `NULL` `err` vector, `err_method = "mc`, and `return_mc = TRUE`).
#'
#' Two graph forms are available:
#'
#' * line/path plots (`type = "line"`), where subsequent iterations of are presented
#' in the X axis and regulation estimates are presented in the Y axis,
#' estimates of the same reaction are connected with lines.
#'
#' * violin plots (`type = "violin"`), where reactions are presented in the
#' Y axis and estimates are presented in the X axis. Distribution of the
#' estimates is presented in violins, estimates in single reactions are shown
#' as points.
#'
#' @return a `ggplot` graphic.
#'
#' @inheritParams plot_errors
#' @param fun a function used to transform the regulation estimates ad their
#' errors prior to plotting. Defaults to `identity`.
#' @param type type of the plot, see Details.
#' @param line_alpha alpha/opacity of lines representing single reactions
#' when `type = "line"`.
#' @param line_width width of lines representing single reactions when `type = "line"`.
#' @param show_zero_line logical should lines representing no regulation be
#' displayed in the plots (`fold_reg = 1`)?
#' @param zero_line_color color of the lines of no regulation.
#' @param zero_line_type type of the no regulation line.
#' @param zero_line_alpha alpha/opacity of the no regulation lines.
#' @param point_wjitter width of jittering of data points.
#' @param point_hjitter height of jittering of data points.
#' @param ... additional arguments passed to \code{\link[ggplot2]{geom_line}} and
#' \code{\link[ggplot2]{geom_violin}}.
#'
#' @export

  plot_mc <- function(x,
                      reactions = NULL,
                      subsystems = NULL,
                      fun = identity,
                      type = c("line", "violin"),
                      point_color = "steelblue",
                      point_size = 2,
                      point_alpha = 0.75,
                      point_wjitter = 0,
                      point_hjitter = 0.1,
                      line_alpha = 1,
                      line_width = 0.25,
                      show_zero_line = TRUE,
                      zero_line_color = "black",
                      zero_line_type = "dashed",
                      zero_line_alpha = 1,
                      cust_theme = microViz::theme_micro(),
                      plot_title = NULL,
                      plot_subtitle = NULL,
                      x_lab = NULL,
                      y_lab = NULL, ...) {

    ## entry control -----------

    if(!is_actiData(x)) stop("`x` has to be an `actiData` object.", call. = FALSE)

    if(!"mc" %in% names(x)) {

      warning("No simulation results in the object.", call. = FALSE)

      return(NULL)

    }

    if(!is_function(fun)) stop("`fun` has to be a function.", call. = FALSE)

    type <- match.arg(type[1], c("line", "violin"))

    stopifnot(is.numeric(point_size))
    point_size <- point_size[1]

    stopifnot(is.numeric(point_alpha))
    point_alpha <- point_alpha[1]

    stopifnot(is.numeric(point_wjitter))
    point_wjitter <- point_wjitter[1]

    stopifnot(is.numeric(point_hjitter))
    point_hjitter <- point_hjitter[1]

    stopifnot(is.numeric(line_alpha))
    line_alpha <- line_alpha[1]

    stopifnot(is.numeric(line_width))
    line_width <- line_width[1]

    stopifnot(is.logical(show_zero_line))
    show_zero_line <- show_zero_line[1]

    stopifnot(is.numeric(zero_line_alpha))
    zero_line_alpha <- zero_line_alpha[1]

    if(!is.null(cust_theme)) {

      if(!is.theme(cust_theme)) {

        stop("`cust_theme` has to be a `ggplot` theme object.", call. = FALSE)

      }

    }

    ## the plotting data ----------

    x <- select(x, reactions, subsystems)

    reactions <- colnames(x[["mc"]])

    plot_data <- as.data.frame(x[["mc"]])

    plot_data[["iter_number"]] <- paste0("iter_", 1:nrow(plot_data))

    plot_data <- pivot_longer(plot_data,
                              cols = all_of(reactions),
                              names_to = "id",
                              values_to = "fold_reg")

    plot_data[["fold_reg"]] <- fun(plot_data[["fold_reg"]])

    ## line plots ---------

    if(type == "line") {

      mc_plot <- ggplot(plot_data,
                        aes(x = .data[["iter_number"]],
                            y = .data[["fold_reg"]],
                            color = .data[["id"]])) +
        geom_line(aes(group = .data[["id"]]),
                  alpha = line_alpha,
                  linewidth = line_width, ...)

      if(show_zero_line) {

       mc_plot <- mc_plot +
         geom_hline(yintercept = fun(1),
                    linetype = zero_line_type,
                    alpha = zero_line_alpha,
                    color = zero_line_color)

      }

      if(is.null(x_lab)) x_lab <- "iteration number"
      if(is.null(y_lab)) y_lab <- "activity estimate"

    }

    ## violin plots --------

    if(type == "violin") {

      mc_plot <- ggplot(plot_data,
                        aes(x = .data[["fold_reg"]],
                            y = reorder(.data[["id"]],
                                        .data[["fold_reg"]]))) +
        geom_violin(...) +
        geom_point(shape = 21,
                   size = point_size,
                   color = "black",
                   fill = point_color,
                   alpha = point_alpha,
                   position = position_jitter(width = point_wjitter,
                                              height = point_hjitter))

      if(show_zero_line) {

        mc_plot <- mc_plot +
          geom_vline(xintercept = fun(1),
                     linetype = zero_line_type,
                     alpha = zero_line_alpha,
                     color = zero_line_color)

      }

      if(is.null(x_lab)) x_lab <- "activity estimate"
      if(is.null(y_lab)) y_lab <- "reaction"

    }

    ## the output -----------

    if(!is.null(cust_theme)) mc_plot <- mc_plot + cust_theme

    mc_plot +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

  }

# END --------
