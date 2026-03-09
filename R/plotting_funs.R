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
                          cust_theme = theme_classic(),
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
                      cust_theme = theme_classic(),
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

# Plotting of numbers of significantly regulate reactions ---------

#' Bar/stack plots with numbers of significantly activated and inhibited reactions.
#'
#' @description
#' Representation of numbers or percentages of differentially regulated
#' reactions in subsystems as bar or stack plot.
#'
#' @details
#' Please make sure that the object contains activity regulation status
#' information.
#' If not, please call \code{\link{identify_regulated}} prior to
#' counting the significant effects.
#'
#' The plot type is determined by `type` argument:
#' `"plus_minus"` (bar plot, inhibited reaction as "negative" numbers),
#' `"bar"` (bar plot, dodged bars for inhibited and activated reactions),
#' `"stack"` plot.
#'
#' @return a `ggplot` graphic.
#'
#' @inheritParams plot_errors
#' @param scale type of the frequency statistic to be shown in the plot:
#' count (default) or percentage of reactions in the subsystem.
#' @param type type of the plot, see Details.
#' @param show_all logical, should frequencies of all regulated reactions be
#' presented in the plot together with counts/percentages for the subsystems?
#' @param palette color palette: a named character vector with at least two elements
#' specifying the colors for activated and inhibited reactions.
#' @param bar_rim_color color of the bar's rim/line.
#' @param show_n_total logical: if `TRUE`, total numbers of reactions in subsystems
#' are displayed in the Y axis along with the subsystem names.
#' @param labeller_fun a function used to transform subsystem names, e.g. into
#' abbreviations, which are shown in the Y axis of the plot.
#' @param fill_lab title of the fill scale.
#' @param ... additional arguments passed to \code{\link[ggplot2]{geom_bar}}.
#'
#' @export

  plot_numbers <- function(x,
                           reactions = NULL,
                           subsystems = NULL,
                           scale = c("count", "percent"),
                           type = c("plus_minus", "bar", "stack"),
                           palette = c("activated" = "firebrick",
                                       "inhibited" = "steelblue"),
                           bar_rim_color = "black",
                           show_all = FALSE,
                           labeller_fun = identity,
                           show_n_total = TRUE,
                           cust_theme = theme_classic(),
                           plot_title = NULL,
                           plot_subtitle = NULL,
                           x_lab = NULL,
                           y_lab = "subsystem",
                           fill_lab = "regulation\nstatus", ...) {

    ## input control --------

    if(!is_actiData(x)) stop("`x` has to be an `actiData` object.", call. = FALSE)

    scale <- match.arg(scale[1], c("count", "percent"))

    type <- match.arg(type[1], c("plus_minus", "bar", "stack"))

    stopifnot(is.logical(show_all))
    show_all <- show_all[1]

    stopifnot(is.logical(show_n_total))
    show_n_total <- show_n_total[1]

    stopifnot(is.character(palette))

    if(length(palette) < 2) {

      stop("`palette` has to have at least two color names.", call. = FALSE)

    }

    palette <- palette[1:2]

    if(!is_function(labeller_fun)) {

      stop("`labeller_fun` has to be a function.", call = FALSE)

    }

    if(!is.null(cust_theme)) {

      if(!is.theme(cust_theme)) {

        stop("`cust_theme` has to be a `ggplot` theme object.", call. = FALSE)

      }

    }

    ## plotting data --------

    x <- select(x, reactions, subsystems)

    count_data <- count_regulated(x)

    if(!show_all) {

      count_data <- filter(count_data, .data[["subsystem"]] != "all")

    }

    plot_variable <- switch(scale,
                            count = "n",
                            percent = "percent")

    count_data[["axis_label"]] <- labeller_fun(count_data[["subsystem"]])

    if(show_n_total) {

      count_data[["axis_label"]] <-
        paste(count_data[["axis_label"]],
              count_data[["n_total"]],
              sep = "\nn = ")

    }

    if(is.null(x_lab)) {

      x_lab <- switch(scale,
                      count = "number of reactions",
                      percent = "% of subsystem reactions")

    }

    ## the plots ----------

    if(type %in% c("bar", "stack")) {

      pos_txt <- switch(type,
                        bar = "dodge",
                        stack = "stack")

      number_plot <- ggplot(count_data,
                            aes(x = .data[[plot_variable]],
                                y = reorder(.data[["axis_label"]],
                                            .data[[plot_variable]]),
                                fill = .data[["regulation"]])) +
        geom_bar(stat = "identity",
                 color = bar_rim_color,
                 position = pos_txt, ...)

    } else {

      number_plot <- ggplot(count_data,
                            aes(x = ifelse(.data[["regulation"]] == "activated",
                                           .data[[plot_variable]],
                                           -.data[[plot_variable]]),
                                y = reorder(.data[["axis_label"]],
                                            ifelse(.data[["regulation"]] == "activated",
                                                   .data[[plot_variable]],
                                                   -.data[[plot_variable]])),
                                fill = .data[["regulation"]])) +
        geom_bar(stat = "identity",
                 color = bar_rim_color, ...) +
        scale_x_continuous(labels = abs)

    }

    if(!is.null(cust_theme)) number_plot <- number_plot + cust_theme

    number_plot +
      geom_vline(xintercept = 0) +
      scale_fill_manual(values = palette,
                        drop = FALSE) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab,
           fill = fill_lab)

  }

# END --------
