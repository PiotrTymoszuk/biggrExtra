# Methods and functions for the geneSBML class

# Class testing -----

#' Test for the geneSBML and SBML class instance.
#'
#' @description
#' The function test whether the object is an instance of the geneSBML class
#' or SBML class.
#'
#' @param x an object.
#' @return a logical value.
#'
#' @export

  is_geneSBML <- function(x) inherits(x, 'geneSBML')

#' @rdname is_geneSBML
#' @export

  is_SBML <- function(x) inherits(x, 'SBML')

# Appearance -----

#' Appearance for the geneSBML class.
#'
#' @description Print method for the geneSBML class.
#' @param x an object of the geneSBML class.
#' @param ... extra arguments, currently none.
#' @export

  print.geneSBML <- function(x, ...) {

    cat(paste('A geneSBML object for', nrow(x$reg), 'reactions\n'))
    print(x$reg)

    invisible(NULL)

  }

# Calculating reaction rates ------

#' Reaction regulation rates.
#'
#' @description
#' Extracts regulation rates and their errors from and geneSBML object.
#'
#' @param x a geneSBML object.
#' @param react_id a character vector of reaction IDs.
#' @param ... extra arguments, currently none.
#' @return a data frame with reaction regulation estimates and, optionally,
#' their errors.
#' @export calculate.geneSBML
#' @export

  calculate.geneSBML <- function(x, react_id = x$reg$react_id, ...) {

    stopifnot(is_geneSBML(x))

    filter(x$reg, .data[['react_id']] %in% .env[['react_id']])

  }

# Accessing model components ------

#' Extract geneSBML model components.
#'
#' @description
#' Extracts the following geneSBML components:
#' SBML model (`type = 'model'`, default),
#' reaction regulation estimates ('regulation'),
#' gene mapping ('gene_map'),
#' all reaction IDs ('reactions'),
#' all metabolite IDs ('metabolites')
#' all gene IDs ('genes') or a data frame
#' with subsystem mapping ('subsystems').
#'
#' @return the output specified by `type` argument.
#'
#' @param object a geneSBML object.
#' @param type type of the output, see above.
#' @param ... extra arguments, currently none.
#'
#' @export components.geneSBML
#' @export

  components.geneSBML <- function(object,
                                  type = c('model', 'regulation',
                                           'gene_map', 'reactions',
                                           'metabolites', 'genes',
                                           'subsystems'),
                                  ...) {

    stopifnot(is_geneSBML(object))

    type <- match.arg(type[1],
                      c('model', 'regulation',
                        'gene_map', 'reactions',
                        'metabolites', 'genes',
                        'subsystems'))

    switch(type,
           model = object$model,
           regulation = object$reg,
           gene_map = object$gene_map,
           reactions = object$reg$react_id,
           metabolites = names(object$model@species),
           genes = reduce(object$gene_map$entrez_id, union),
           subsystems = extract_subsystems(object, as_list = FALSE))

  }

# Plotting -------

#' Plot features of a geneSBML object.
#'
#' @description
#' Plots features of the geneSBML object as specified by the `type`
#' argument.
#' For `type` set to 'regulation' a bar plot of the
#' regulation values is created, 'errors' returns a point plot of regulation
#' estimates and regulation errors if errors can be retrieved from the model.
#' For `type` set to 'top' a bar plot with top up- and downregulated reactions
#' is returned.
#' For 'volcano' a Volcano plot of raw or FDR-corrected p values and reaction
#' regulation estimates is returned.
#' The 'top_forest' draws a Forest plot of regulation estimates
#' along with 95% confidence intervals for the top activated and
#' inhibited reactions.
#' The 'forest' option draws a Forest plot of regulation estimates
#' along with 95% confidence intervals.
#' The 'bar' option generates a bar plot of regulation estimates with bar color
#' coding for significance and regulation status.
#' Finally, the type argument set to 'mc' draws Monte Carlo estimates of
#' pathway regulation in subsequent runs; this options is available only for
#' geneSBML model with errors estimated by the Monte Carlo method without
#' memory saving option.
#'
#' @return a ggplot2 graphic object.
#'
#' @param x an object of the geneSBML class.
#' @param type type of the plot, as specified above.
#' @param relevant.reactions a vector of reaction IDs to be presented
#' in the plots, defaults to all available reactions.
#' @param n_top number of top up- and downregulated reactions to be plotted
#' (type 'top' and 'forest') or top significant reactions to be labeled in
#' the Volcano plot.
#' @param show_fit logical, should a trend line be plotted? valid only for
#' error and regulation plots and ignored otherwise.
#' @param point_size size of the points in the errors, Volcano or Forest plots.
#' Ignored for other plot types'.
#' @param point_color color of the points in the error plots.
#' Ignored for other plot types.
#' @param point_alpha alpha of the points in the error and volcano plots.
#' Ignored for other plot types.
#' @param line_alpha alpha of the lines connecting estimates in the MC plots.
#' Ignored for other plot types.
#' @param signif_type significance to be plotted in Volcano or Forest plots:
#' 'raw' or 'fdr' (FDR-corrected, default).
#' @param regulation_level fold-regulation level cutoff to be presented in
#' Forest and Volcano plots, defaults to 1.
#' @param label_names logical, should the reactions be labeled with their
#' full names retrieved from the model? Defaults to FALSE.
#' @param order_reactions logical, should the reactions be ordered in the plot
#' by their regulation estimates? Applies only to the regulation bar plot
#' and defaults to FALSE. If FALSE, the reactions are plotted as specified in
#' 'relavant.reactions'.
#' @param cust_theme custom ggplot2 theme.
#' @param ... extra arguments: for 'errors' arguments controlling the
#' appearance of the trend line passed to \code{\link[ggplot2]{geom_smooth}};
#' for Volcano: arguments passed to \code{\link[microViz]{plot_volcano}};
#' for Forest plots: arguments passed to  \code{\link[microViz]{plot_top}} or
#' \code{\link[microViz]{plot_forest}}.
#'
#' @export plot.geneSBML
#' @export

  plot.geneSBML <- function(x,
                            type = c('regulation',
                                     'errors',
                                     'top',
                                     'volcano',
                                     'forest',
                                     'bar',
                                     'top_forest',
                                     'mc'),
                            relevant.reactions = components(x, 'reactions'),
                            n_top = 10,
                            show_fit = TRUE,
                            point_size = 2,
                            point_color = 'black',
                            point_alpha = 0.25,
                            line_alpha = 1,
                            signif_type = c('fdr', 'raw'),
                            regulation_level = 1,
                            label_names = FALSE,
                            order_reactions = FALSE,
                            cust_theme = ggplot2::theme_classic(), ...) {

    ## entry control -------

    stopifnot(is_geneSBML(x))
    stopifnot(is.logical(show_fit))
    stopifnot(is.numeric(n_top))
    stopifnot(is.logical(show_fit))
    stopifnot(is.numeric(point_size))
    stopifnot(is.numeric(point_alpha))
    stopifnot(is.numeric(regulation_level))
    stopifnot(is.logical(label_names))
    stopifnot(is.numeric(line_alpha))
    stopifnot(is.logical(order_reactions))

    fold_reg <- NULL
    react_id <- NULL
    reg_sign <- NULL
    significant <- NULL
    lower_ci <- NULL
    upper_ci <- NULL
    reg_split <- NULL
    error <- NULL
    run <- NULL

    type <- match.arg(type[1],
                      c('regulation',
                        'errors',
                        'top',
                        'volcano',
                        'forest',
                        'bar',
                        'top_forest',
                        'mc'))

    n_top <- as.integer(n_top)

    if(!ggplot2::is.theme(cust_theme)) {

      stop('cust_theme has to be a valid ggplot2 theme object.',
           call. = FALSE)

    }

    signif_type <- match.arg(signif_type[1],
                             c('raw', 'fdr'))

    signif_var <- switch(signif_type,
                         raw = 'p_value',
                         fdr = 'p_adjusted')

    p_lab <- switch(signif_type,
                    raw = expression('-log'[10] * ' p'),
                    fdr = expression('-log'[10] * ' pFDR'))

    if(label_names) {

      react_names <-
        annotate_bigg(bigg_id = components(x, 'reactions'),
                      value = 'name',
                      annotation_db = x)

    }

    miss_reactions <-
      relevant.reactions[!relevant.reactions %in% components(x, 'reactions')]

    if(length(miss_reactions) > 0) {

      warning('Some reactions missing from the model.',
              call. = FALSE)

    }

    if(type == 'mc' & inherits(x, 'memoSaver')) {

      stop(paste('Plots of Monte Carlo simulation results are not avalable",
                  "for models with the memonry saving option.'),
            call. = FALSE)

    }

    ## plotting table ------

    plot_tbl <- dplyr::mutate(x$reg,
                              fold_reg = log2(fold_reg))

    plot_tbl <- dplyr::filter(plot_tbl,
                              react_id %in% relevant.reactions)

    if(!'error' %in% names(plot_tbl)) {

      plot_tbl <-
        mutate(plot_tbl,
               reg_sign = ifelse(fold_reg > log2(regulation_level),
                                 'activated',
                                 ifelse(fold_reg < -log2(regulation_level),
                                        'inhibited', 'constant')),
               reg_sign = factor(reg_sign,
                                 c('activated', 'inhibited', 'constant')))

      fill_vals <- c(activated = 'firebrick',
                     inhibited = 'steelblue',
                     constant = 'gray60')

    } else {

      plot_tbl <-
        mutate(plot_tbl,
               significant = ifelse(.data[[signif_var]] < 0.05,
                                    'yes', 'no'),
               reg_sign = ifelse(significant == 'no',
                                 'ns',
                                 ifelse(fold_reg > log2(regulation_level),
                                        'activated',
                                        ifelse(fold_reg < -log2(regulation_level),
                                               'inhibited', 'ns'))),
               reg_sign = factor(reg_sign,
                                 c('activated', 'inhibited', 'ns')),
               lower_ci = log2(lower_ci),
               upper_ci = log2(upper_ci),
               lower_ci = ifelse(is.nan(lower_ci), fold_reg, lower_ci),
               upper_ci = ifelse(is.nan(upper_ci), fold_reg, upper_ci))

      fill_vals <- c(activated = 'firebrick',
                     inhibited = 'steelblue',
                     ns = 'gray60')

    }

    ## plotting: top regulated plot -------

    if(type == 'top') {

      plot_tbl <- filter(plot_tbl, fold_reg != 0)

      plot_tbl <-
        mutate(plot_tbl, reg_split = sign(fold_reg))

      plot_tbl <- group_by(plot_tbl, reg_split)

      plot_tbl <- top_n(plot_tbl, n = n_top, abs(fold_reg))

      plot_tbl <- ungroup(plot_tbl)

      top_plot <-
        ggplot(plot_tbl,
               aes(x = fold_reg,
                   y = reorder(react_id, fold_reg),
                   fill = reg_sign)) +
        geom_vline(xintercept = 0,
                   linetype = 'dashed') +
        geom_bar(stat = 'identity',
                 color = 'black') +
        scale_fill_manual(values = fill_vals,
                          name = 'Reaction status') +
        cust_theme +
        theme(axis.title.y = element_blank()) +
        labs(title = paste('Top', n_top, 'regulated reactions'),
             x = expression('log'[2] * ' fold regulation'))

      if(label_names) {

        top_plot <- top_plot +
          scale_y_discrete(labels = react_names)

      }

      return(top_plot)

    }

    ## information of total reaction numbers and their regulation

    total_pathways <- nrow(plot_tbl)

    total_tag <- paste('reations: total: n =', total_pathways)

    reg_counts <-
      count(filter(plot_tbl,
                   reg_sign %in% c('activated', 'inhibited')),
            reg_sign)

    reg_tag <- map2_chr(reg_counts[[1]],
                        reg_counts[[2]],
                        paste, sep = ': n = ')

    plot_cap = paste(total_tag,
                     paste(reg_tag, collapse = ', '),
                     sep = ', ')

    ## regulation plot -------

    if(type == 'regulation') {

      reg_plot <-
        ggplot(plot_tbl,
                        aes(x = reorder(react_id, fold_reg),
                                     y = fold_reg)) +
        geom_bar(aes(fill = reg_sign,
                              color = reg_sign),
                          stat = 'identity') +
        scale_fill_manual(values = fill_vals,
                                   name = 'Reaction status') +
        scale_color_manual(values = fill_vals,
                                    name = 'Reaction status') +
        cust_theme +
        theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()) +
        labs(title = 'Reaction regulation',
                      subtitle = plot_cap,
                      y = expression('log'[2] * ' fold regulation'),
                      x = 'Reaction')

      if(show_fit) {

        reg_plot <- reg_plot +
          geom_smooth(show.legend = FALSE, ...)

      }

      return(reg_plot)

    }

    ## errors plot: fold regulation versus error -------

    if(type == 'errors') {

      err_plot <-
        ggplot(plot_tbl,
               aes(x = fold_reg,
                   y = log2(error))) +
        geom_vline(xintercept = 0,
                   linetype = 'dashed') +
        geom_point(shape = 16,
                   size = point_size,
                   color = point_color,
                   alpha = point_alpha) +
        cust_theme +
        labs(title = 'Fold regulation and regulation error',
             subtitle = plot_cap,
             x = expression('log'[2] * ' fold regulation'),
             y = expression('log'[2] * ' regulation error'))

      if(show_fit) {

        err_plot <- err_plot +
          geom_smooth(show.legend = FALSE, ...)

      }

      return(err_plot)

    }

    ## volcano plot -------

    if(type == 'volcano') {

      if(!'error' %in% names(plot_tbl)) {

        warning(paste('Volcano plots are available only for',
                      'geneSBML objects with error estimates.'),
                call. = FALSE)

        return(NULL)

      }

      return(microViz::plot_volcano(data = plot_tbl,
                                    regulation_variable = 'fold_reg',
                                    p_variable = signif_var,
                                    signif_level = 0.05,
                                    regulation_level = log2(regulation_level),
                                    y_lab = p_lab,
                                    x_lab = expression('log'[2] * ' fold-regulation'),
                                    label_variable = 'react_id',
                                    top_significant = n_top,
                                    fill_title = 'Reaction status',
                                    plot_title = 'Reaction regulation and significance',
                                    cust_theme = cust_theme, ...))

    }

    ## Top Forest plot -------

    if(type == 'top_forest') {

      if(!'error' %in% names(plot_tbl)) {

        warning(paste('Forest plots are available only for geneSBML',
                      'objects with error estimates.'),
                call. = FALSE)

        return(NULL)

      }

      forest_plot <-
        microViz::plot_top(data = plot_tbl,
                           regulation_variable = 'fold_reg',
                           p_variable = signif_var,
                           signif_level = 0.05,
                           regulation_level = log2(regulation_level),
                           lower_ci_variable = 'lower_ci',
                           upper_ci_variable = 'upper_ci',
                           top_regulated = n_top,
                           plot_title = paste('Top',
                                              n_top,
                                              'regulated reactions'),
                           x_lab = expression('log'[2] * ' fold-regulation'),
                           label_variable = 'react_id',
                           fill_title = 'Reaction status',
                           cust_theme = cust_theme, ...)

      if(label_names) {

        forest_plot <- forest_plot +
          scale_y_discrete(labels = react_names)

      }

      return(forest_plot)

    }

    ## General Forest plot -----

    if(type == 'forest') {

      if(!'error' %in% names(plot_tbl)) {

        warning(paste('Forest plots are available only for geneSBML objects',
                      'with error estimates.'),
                call. = FALSE)

        return(NULL)

      }

      forest_plot <-
        microViz::plot_forest(data = plot_tbl,
                              regulation_variable = 'fold_reg',
                              label_variable = 'react_id',
                              p_variable = signif_var,
                              signif_level = 0.05,
                              regulation_level = log2(regulation_level),
                              lower_ci_variable = 'lower_ci',
                              upper_ci_variable = 'upper_ci',
                              plot_title = paste('Reaction regulation',
                                                 'estimates'),
                              x_lab = expression('log'[2] * ' fold-regulation'),
                              fill_title = 'Reaction status',
                              cust_theme = cust_theme, ...)

      if(label_names) {

        forest_plot <- forest_plot +
          scale_y_discrete(labels = react_names)

      }

      return(forest_plot)

    }

    ## regulation bar plot -------

    if(type == 'bar') {

      if(order_reactions) {

        bar_plot <-
          ggplot(plot_tbl,
                 aes(x = fold_reg,
                     y = reorder(react_id, fold_reg),
                     fill = reg_sign))

      } else {

        plot_tbl <- mutate(plot_tbl,
                           react_id = factor(react_id,
                                             relevant.reactions))

        bar_plot <-
          ggplot(plot_tbl,
                 aes(x = fold_reg,
                     y = react_id,
                     fill = reg_sign))

      }

      bar_plot <- bar_plot +
        geom_bar(stat = 'identity',
                 color = 'black') +
        scale_fill_manual(values = fill_vals,
                          name = 'Status') +
        cust_theme +
        theme(axis.title.y = element_blank()) +
        labs(title = 'Reaction regulation estimates',
             x = expression('log'[2] * ' fold-regulation'))

      if(label_names) {

        bar_plot <- bar_plot +
          scale_y_discrete(labels = react_names)

      }

      return(bar_plot)

    }

    ## Monte Carlo runs -------

    if(type == 'mc') {

      if(is.null(x$mc)) {

        warning('MC plots are available only for model with Monte Carlo errors.',
                call = FALSE)

      }

      relevant.reactions <-
        relevant.reactions[relevant.reactions %in% colnames(x$mc)]

      plot_tbl <- x$mc[, relevant.reactions]

      plot_tbl <- mutate(as.data.frame(plot_tbl),
                         run = 1:nrow(plot_tbl))

      plot_tbl <- tidyr::pivot_longer(plot_tbl,
                                      cols = all_of(relevant.reactions),
                                      names_to = 'react_id',
                                      values_to = 'fold_reg')

      run_plot <-
        ggplot(plot_tbl,
               aes(x = run,
                   y = fold_reg,
                   color = react_id))  +
        geom_line(aes(group = react_id),
                  alpha = line_alpha) +
        cust_theme +
        labs(title = 'Reaction regulation, Monte Carlo simulation',
             y = 'Fold-regulation',
             x = 'MC iteration')

      if(show_fit) {

        run_plot <- run_plot +
          geom_smooth(se = FALSE, ...)

      }

      return(run_plot)

    }

  }

# Counts of regulated reactions ------

#' Count regulated reactions per subsystem.
#'
#' @description
#' Computes counts of significantly activated, inhibited and
#' all reactions per subsystem.
#'
#' @return a data frame with reaction counts according
#' to their regulation status.
#'
#' @param x an object of the geneSBML class.
#' @param signif_type type of significance to be used for definition of
#' regulated reactions:
#' 'raw' or 'fdr' (FDR-corrected, default).
#' Ignored if the model contains no error information.
#' @param regulation_level fold-regulation level cutoff used for definition
#' of regulated reactions.
#' @param ... extra arguments, currently none.
#'
#' @export count.geneSBML
#' @export

  count.geneSBML <- function(x,
                             signif_type = c('fdr', 'raw'),
                             regulation_level = 1, ...) {

    ## entry control -------

    fold_reg <- NULL
    react_id <- NULL

    stopifnot(is_geneSBML(x))
    stopifnot(is.numeric(regulation_level))

    signif_type <- match.arg(signif_type[1],
                             c('raw', 'fdr'))

    signif_var <- switch(signif_type,
                         raw = 'p_value',
                         fdr = 'p_adjusted')

    reg_tbl <- mutate(x$reg, fold_reg = log2(fold_reg))

    ## definition of regulated reactions ------

    ### working with log2-transformation of the regulation estimates
    ### simplifies the coding

    reg_sign <- NULL
    significant <- NULL
    subsystem <- NULL

    if(!'error' %in% names(reg_tbl)) {

      reg_tbl <-
        mutate(reg_tbl,
               reg_sign = ifelse(fold_reg > log2(regulation_level),
                                 'activated',
                                 ifelse(fold_reg < -log2(regulation_level),
                                        'inhibited', 'constant')),
               reg_sign = factor(reg_sign,
                                 c('activated', 'inhibited', 'constant')))

    } else {

      reg_tbl <-
        mutate(reg_tbl,
               significant = ifelse(.data[[signif_var]] < 0.05,
                                    'yes', 'no'),
               reg_sign = ifelse(significant == 'no',
                                 'ns',
                                 ifelse(fold_reg > log2(regulation_level),
                                        'activated',
                                        ifelse(fold_reg < -log2(regulation_level),
                                               'inhibited', 'ns'))),
               reg_sign = factor(reg_sign,
                                 c('activated', 'inhibited', 'ns')))

    }

    ## counting: all reactions and as per subsystem --------

    n <- NULL

    all_counts <- count(reg_tbl, reg_sign, .drop = FALSE)

    all_counts <- mutate(all_counts, subsystem = 'All reactions')

    sub_counts <- count(reg_tbl, subsystem, reg_sign, .drop = FALSE)

    sub_counts <- rbind(all_counts, sub_counts)

    ## computing total reaction counts per subsystem

    total_counts <- group_by(sub_counts, subsystem)

    total_counts <- summarise(total_counts, n_total = sum(n))

    sub_counts <-
      left_join(sub_counts,
                total_counts,
                by = 'subsystem')

    ## output -------

    set_names(sub_counts[c('subsystem', 'reg_sign', 'n', 'n_total')],
              c('subsystem', 'status', 'n', 'n_total'))

  }

# Fitting a LIM model -------

#' Fit a LIM based on geneSBML gene regulation data.
#'
#' @description
#' Fits a linear inverse model (LIM) to the SBML metabolic
#' reaction network model and reaction regulation estimates stored in the
#' geneSBML object.
#' Technically, a wrapper around \code{\link[BiGGR]{createLIMFromSBML}} with an
#' additional entry check. See the documentation of the genuine BiGGR function
#' for details. Reaction regulation estimates are fed into the LIM as
#' 'equations' argument.
#'
#' @return A model file with with extension ".lim" is created.
#'
#' @param object a geneSBML object.
#' @param maximize a reaction or a series thereof to be maximized provided as
#' a single string or a character vector. Each reaction specified has
#' to be present in the regulation data stored in the geneSBML model.
#' @param inequalities a list specifying inequality constraints on the system,
#' see: \code{\link[BiGGR]{createLIMFromSBML}}.
#' @param constraints model constraints.
#' @param externals a character vector of metabolites for which the flux
#' balance analysis (FBA) has to be performed.
#' @param file.name a path and name of the LIM file.
#' @param ... extra arguments, currently none.
#'
#'
#' @export fit.geneSBML
#' @export

  fit.geneSBML <- function(object,
                           maximize,
                           inequalities,
                           constraints,
                           externals,
                           file.name = 'model.lim', ...) {

    ## entry control

    stopifnot(is_geneSBML(object))

    if(inherits(object, 'memoSaver')) {

      stop('Unable to construct LIM from a memory-saving model.',
           call. = FALSE)

    }

    max_reactions <- stri_extract_all(maximize, regex = 'R_\\w+')

    max_reactions <- unlist(max_reactions)

    if(length(max_reactions) == 0) {

      stop('No reactions specified in maximized were mapped to the model data.',
           call. = FALSE)

    }

    missing_max <- max_reactions[!max_reactions %in% object$reg$react_id]

    if(length(missing_max) > 0) {

      stop(paste('The following reactions specified in maximized are missing from the model data:',
                 paste(missing_max, sep = ', ')),
           call. = FALSE)

    }

    ## fitting the model

    reg_tbl <-
      filter(object$reg[c('react_id', 'fold_reg')],
             complete.cases(object$reg[c('react_id', 'fold_reg')]))

    eqn <- as.list(reg_tbl)

    createLIMFromSBML(model = object$model,
                      maximize = maximize,
                      equations = eqn,
                      inequalities = inequalities,
                      constraints = constraints,
                      externals = externals,
                      file.name = file.name)

  }

# Visualization -------

#' Create a hypergraph from a geneSBML object.
#'
#' @description
#' Generates a hypergraph (instance of the 'RagraphBPH' class)
#' for the requested reaction and metabolites.
#'
#' @details
#' Technically, a wrapper around \code{\link[BiGGR]{sbml2hyperdraw}}.
#' Rates for particular reactions and their errors are retrieved from the
#' regulation table stored within the object.
#' Reactions and their rates are presented in the edges. In addition, those
#' reaction labels may be appended by user-specified suffixes such as
#' confidence intervals or p values.
#' Edge width corresponds to absolute log~2~ fold regulation. The edge color is
#' specified by raw or FDR-corrected significance and the regulation sign.
#'
#' @param x a geneSBML object.
#' @param rate_sep a character separating the reaction identifier
#' and its regulation rate in the graph.
#' @param suffixes additional information to be displayed in the graph along
#' with the reactions and their rates: none ('none', default),
#' error ('error'), raw p value ('p_raw')
#' or FDR-adjusted p value ('p_fdr').
#' @param suffix_sep a character separating the reaction identifier/rate
#' from the suffix in the in the graph.
#' @param colors a named vector of colors for the reactions
#' (names: 'activated', 'inhibited', 'ns').
#' @param node_color default node color.
#' @param signif_type significance to be plotted in Volcano or Forest plots:
#' 'raw' (default) or 'fdr' (FDR-corrected).
#' @param signif_digits significant digits used for rounding of numeric values
#' presented in the plot.
#' @param fontsize size of the font, in cex-units.
#' @param relevant.species metabolites to be plotted. If not specified by
#' the user, all metabolites participating in the reactions specified by
#' 'relevant.reactions' are displayed.
#' @param relevant.reactions BiGG ID of the reactions to be plotted.
#' With or without leading 'R_' string.
#' @param layoutType is a character string representing the layout engine to
#' be used for visualization. Current supported layouts are "dot", "twopi",
#' "neato","fdp","sfdp" and "circo". Defaults to "dot".
#' @param lwd.max a numeric given the maximum edge width. Defaults to 3.
#' @param lwd.min a numeric given the minimum edge width. Defaults to 0.5
#' @param plt.margins A numerical vector of the form
#' c(bottom, left, top, right) giving additional white space around the graph
#' (in case long node or edge labels fall outside the plotting region).
#' Defaults to c(150,150,150,150).
#' @param ... extra arguments, currently none.
#'
#'
#' @export visualize.geneSBML
#' @export

  visualize.geneSBML <- function(x,
                                 rate_sep = ':',
                                 suffixes = c('none', 'error',
                                              'p_raw', 'p_fdr'),
                                 suffix_sep = ' ',
                                 colors = c(activated = 'firebrick',
                                            inhibited = 'steelblue',
                                            ns = 'gray70'),
                                 node_color = 'gray30',
                                 signif_type = c('raw', 'fdr'),
                                 signif_digits = 2,
                                 fontsize = 0.75,
                                 relevant.species = components(x, 'metabolites'),
                                 relevant.reactions = components(x, 'reactions'),
                                 layoutType = 'dot',
                                 lwd.max = 3,
                                 lwd.min = 0.5,
                                 plt.margins = c(150, 150, 150, 150), ...) {

    ## entry control -------

    stopifnot(is_geneSBML(x))

    if(inherits(x, 'memoSaver')){

      stop('Unable to visualize metabolic pathways in a memory-saving model.',
           call. = FALSE)

    }

    error <- NULL
    lower_ci <- NULL
    upper_ci <- NULL
    p_value <- NULL
    p_adjusted <- NULL
    fold_reg <- NULL
    color_var <- NULL


    suffixes <- match.arg(suffixes[1],
                          c('none', 'error',
                            'p_raw', 'p_fdr'))

    if(any(!names(colors) %in% c('activated', 'inhibited', 'ns'))) {

      stop('colors have to consist of three elements named activated, inhibited and ns',
           call. = FALSE)

    }

    signif_type <- match.arg(signif_type[1],
                             c('raw', 'fdr'))

    signif_var <- switch(signif_type,
                         raw = 'p_value',
                         fdr = 'p_adjusted')

    ## relevant species and metabolites ------

    relevant.species <-
      ifelse(stri_detect(relevant.species, regex = '^M_'),
             relevant.species,
             paste0('M_', relevant.species))

    relevant.reactions <-
      ifelse(stri_detect(relevant.reactions, regex = '^R_'),
             relevant.reactions,
             paste0('R_', relevant.reactions))

    ## rates, suffixes and colors ------

    rate_vec <- set_names(x$reg$fold_reg, x$reg$react_id)

    suff_tbl <-
      mutate(x$reg,
             error = signif(error, signif_digits),
             ci = paste0('(', signif(lower_ci, signif_digits),
                         ' - ', signif(upper_ci, signif_digits), ')'),
             p_raw = signif(p_value, signif_digits),
             p_fdr = signif(p_adjusted, signif_digits),
             color_var = ifelse(.data[[signif_var]] >= 0.05,
                                'ns',
                                ifelse(fold_reg > 1,
                                       'activated', 'inhibited')),
             color_var = ifelse(is.na(color_var), 'ns', color_var),
             color_var = colors[color_var])

    if(suffixes == 'none') {

      suff_vec <- set_names(rep('', nrow(suff_tbl)),
                            suff_tbl$react_id)

      suffix_sep <- ' '



    } else {

      suff_vec <- set_names(suff_tbl[[suffixes]], suff_tbl$react_id)

    }

    color_vec <- set_names(suff_tbl$color_var, suff_tbl$react_id)

    ## graph object

    sbml_to_hd(model = x$model,
               rates = rate_vec,
               rate_sep = rate_sep,
               suffixes = suff_vec,
               suffix_sep = suffix_sep,
               colors = color_vec,
               node_color = node_color,
               signif_digits = signif_digits,
               fontsize = fontsize,
               relevant.species = relevant.species,
               relevant.reactions = relevant.reactions,
               layoutType = layoutType,
               lwd.max = lwd.max,
               lwd.min = lwd.min,
               plt.margins = plt.margins)

  }

# END ------
