# Utility functions (no-exported) for testing of significant
# enrichment of significantly regulated reactions within metabolic subsystems.

# Fisher's exact test ------

#' Test for enrichment in a subsystem with Fisher's exact test.
#'
#' @description Tests for significant enrichment of generally significantly
#' regulated, significantly activated and significantly inhibited reactions in
#' a metabolic subsystem.
#' The significance is determined by Fisher's exact test.
#' @details The method is biased towards subsystems with multiple reactions!
#' Technically, Fisher's exact tests are computed with
#' \code{\link[rstatix]{fisher_test}}.
#' @param data a data frame with the following variables:
#' 'subsystem', 'status', 'n' (number of regulated reactions) and
#' 'n_total' (total number of reactions in a subsystem).
#' @param alternative indicates the alternative hypothesis
#' and must be one of "two.sided", "greater" or "less". Defaults to "greater".
#' @param ... extra arguments passed to \code{\link[rstatix]{fisher_test}}.
#' @return a data frame with subsystem name and testing results
#' for regulated, activated and inhibited reactions.

  sub_fisher_ <- function(data, alternative = 'greater', ...) {

    ## entry control ------

    ## done be the upstream function!

    ## analysis data ------

    ### in case there are multiple subsystems in 'data':

    data <- dplyr::filter(data, subsystem == data$subsystem[[1]])

    data <- rbind(data,
                  tibble::tibble(subsystem = data$subsystem[[1]],
                                 status = 'regulated',
                                 n = sum(data$n),
                                 n_total = data$n_total[[1]]))

    data <- dplyr::mutate(data,
                          frac_total = n/n_total)

    ##  test matrices ------

    tst_mtx <- plyr::dlply(data, 'status', dplyr::select, n, n_total)

    tst_mtx <- purrr::map(tst_mtx,
                          function(x) rbind(c(activated = x$n[1],
                                              constant = x$n_total[1] - x$n[1]),
                                            c(activated = 0,
                                              constant = x$n_total[1])))

    for(i in names(tst_mtx)) {

      rownames(tst_mtx[[i]]) <- c('model', 'expected')

    }

    ## Fisher tests -------

    test_res <- purrr::map(tst_mtx,
                           ~rstatix::fisher_test(xtab = as.table(.x),
                                                 alternative = alternative,
                                                 ...))

    test_res <-
      purrr::map2_dfr(test_res,
                      names(test_res),
                      ~dplyr::mutate(.x, status = .y))

    ## output --------

    test_res <- dplyr::mutate(test_res, p_value = p)

    test_res <- dplyr::select(test_res, -n, -p, -p.signif)

    dplyr::left_join(data, test_res, by = 'status')

  }

#' Test for enrichment in subsystems with Fisher's exact test.
#'
#' @description Tests for significant enrichment of generally significantly
#' regulated, significantly activated and significantly inhibited reactions in
#' metabolic subsystems.
#' The significance is determined by Fisher's exact test
#' (H0: no significantly regulated reactions)
#' @details The method is biased towards subsystems with multiple reactions!
#' Technically, Fisher's exact tests are computed with
#' \code{\link[rstatix]{fisher_test}}.
#' @inheritParams sub_fisher
#' @param .parallel logical, should the computation be run in parallel?
#' Defaults to FALSE.
#' @return a data frame with subsystem name and testing results
#' for regulated, activated and inhibited reactions
#' as well as FDR-adjusted p values.

  sub_fisher <- function(data,
                         alternative = 'greater',
                         .parallel = FALSE, ...) {

    ## entry control ------

    start_time <- Sys.time()
    on.exit(future::plan('sequential'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    no_subsystems <- length(unique(data$subsystem))

    message(paste('Fisher testing for enrichment in n =',
                  no_subsystems, 'subsystems'))

    alternative <- match.arg(alternative[[1]],
                             c('greater', 'less', 'two.sided'))

    if(!is.data.frame(data)) {

      stop("'data' hast to be a data frame with the 'subsystem', 'status', 'n' and 'n_total' variables.",
           call. = FALSE)

    }

    if(!all(c('subsystem', 'status', 'n', 'n_total') %in% names(data))) {

      stop("'data' hast to be a data frame with the 'subsystem', 'status', 'n' and 'n_total' variables.",
           call. = FALSE)

    }

    stopifnot(is.logical(.parallel))

    data <- dplyr::filter(data, status %in% c('activated', 'inhibited'))

    if(any(!unique(data$status) %in% c('activated', 'inhibited'))) {

      stop("The 'status' variable of 'data' has to have 'activated' and 'inhibited' levels.",
           call. = FALSE)

    }

    ## serial testing -------

    data <- dplyr::mutate(data,
                          subsystem = factor(subsystem,
                                             unique(data$subsystem)))

    sub_lst <- plyr::dlply(data, 'subsystem')

    if(.parallel) {

      future::plan('multisession')

      tst_results <-
        furrr::future_map(sub_lst,
                          ~sub_fisher_(data = .x,
                                       alternative = alternative,
                                       ...),
                          .options = furrr::furrr_options(seed = TRUE))

    } else {

      tst_results <-
        map(sub_lst,
            ~sub_fisher_(data = .x,
                         alternative = alternative,
                         ...))

    }

    tst_results <- do.call('rbind', tst_results)

    dplyr::mutate(dplyr::as_tibble(tst_results),
                  frac_total = n/n_total,
                  p_adjusted = p.adjust(p_value, 'BH'))

  }

# Random reaction draws ------

#' Random reaction draws.
#'
#' @description Draws randomly i, j and k reactions from all reactions
#' available in the model, where i, j and k correspond to the numbers
#' of activated, inhibited and regulated reactions.
#' @return a list of data frames with counts of activated, inhibited
#' and regulated reactions in each subsystem.
#' @param model a geneSBML model.
#' @param signif_type type of significance to be used for definition of
#' regulated reactions:
#' 'raw' (default) or 'fdr' (FDR-corrected).
#' Ignored if the model contains no error information.
#' @param regulation_level fold-regulation level cutoff used for definition
#' of regulated reactions.
#' @param n_iter number of random draws from the reaction set.
#' @param aggregate logical, should the results be aggregated
#' to a single data frame?
#' @param .parallel logical, should the computation be run in parallel?


  draw_reactions <- function(model,
                             signif_type = c('fdr', 'raw'),
                             regulation_level = 1,
                             n_iter = 1000,
                             aggregate = FALSE,
                             .parallel = FALSE) {

    ## entry control ------

    on.exit(future::plan('sequential'))

    if(!is_geneSBML(model)) {

      stop('A valid geneSBML model is required!')

    }

    signif_type <- match.arg(signif_type[[1]],
                             c('fdr', 'raw'))

    stopifnot(is.numeric(n_iter))
    stopifnot(is.logical(aggregate))
    stopifnot(is.logical(.parallel))

    n_iter <- as.integer(n_iter)

    ## reaction vector and counts ------

    all_reacts <- components(model, 'reactions')

    react_map <- components(model, 'regulation')

    react_map <- rlang::set_names(react_map$subsystem,
                                  react_map$react_id)

    reg_counts <- count(model,
                        signif_type = signif_type,
                        regulation_level = regulation_level)

    reg_counts <- dplyr::filter(reg_counts,
                                subsystem == 'All reactions')

    reg_counts <- dplyr::filter(reg_counts,
                                status %in% c('activated', 'inhibited'))

    reg_counts <- reg_counts[c('status', 'n')]

    total_counts <- tibble::tibble(status = 'regulated',
                                   n = sum(reg_counts$n))

    reg_counts <- rbind(reg_counts,
                        total_counts)

    reg_counts <- plyr::dlply(reg_counts, 'status', function(x) x$n)

    ## random draws ------

    draw_names <- paste0('draw_', 1:n_iter)

    draws <-
      purrr::map(reg_counts,
                 function(status) purrr::map(draw_names,
                                             ~sample(all_reacts,
                                                     size = status,
                                                     replace = FALSE)))

    draws <- purrr::map(draws,
                        rlang::set_names,
                        draw_names)

    ## data frames with random draws and subsystem counts ------

    if(.parallel) future::plan('multisession')

    draws <-
      furrr::future_map(draws,
                        ~purrr::map(.x,
                                    ~tibble::tibble(react_id = .x,
                                                    subsystem = react_map[.x])))
    draws <-
      furrr::future_map(draws,
                        ~purrr::map(.x,
                                    mutate,
                                    subsystem = factor(subsystem,
                                                       sort(unique(unname(react_map))))))

    draws <- furrr::future_map(draws,
                               ~purrr::map(.x,
                                           count, subsystem, .drop = FALSE))

    if(!aggregate) {

      return(draws)

    }

    ## aggregating: via matrix fo speed -----

    draws <- purrr::map(draws,
                        ~purrr::map(.x,
                                    tibble::column_to_rownames,
                                    'subsystem'))

    draws <- purrr::map(draws,
                        ~do.call('cbind', .x))

    draws <- purrr::map(draws,
                        tibble::rownames_to_column,
                        'subsystem')

    draws <- purrr::map(draws,
                        rlang::set_names,
                        c('subsystem', draw_names))

    purrr::map(draws, dplyr::as_tibble)

  }


#' Test for enrichment against random reaction draws.
#'
#' @description Tests for enrichment of significantly regulated reactions
#' within metabolic subsystems by comparing the actual counts with
#' n random draws from the reaction set.
#' @param model a geneSBML model.
#' @param signif_type type of significance to be used for definition of
#' regulated reactions:
#' 'raw' (default) or 'fdr' (FDR-corrected).
#' Ignored if the model contains no error information.
#' @param regulation_level fold-regulation level cutoff used for definition
#' of regulated reactions.
#' @param n_iter number of random draws from the reaction set.
#' @param aggregate logical, should the results be aggregated
#' to a single data frame?
#' @param .parallel logical, should the computation be run in parallel?

  sub_simulation <- function(model,
                             signif_type = c('fdr', 'raw'),
                             regulation_level = 1,
                             n_iter = 1000,
                             .parallel = FALSE) {

    ## entry control -------

    start_time <- Sys.time()
    on.exit(future::plan('sequential'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    on.exit(future::plan('sequential'))

    if(!is_geneSBML(model)) {

      stop('A valid geneSBML model is required!')

    }

    signif_type <- match.arg(signif_type[[1]],
                             c('fdr', 'raw'))

    stopifnot(is.numeric(n_iter))
    stopifnot(is.logical(.parallel))

    n_iter <- as.integer(n_iter)

    ## counts of regulated reactions ------

    react_counts <- count(model,
                          signif_type = signif_type,
                          regulation_level = regulation_level)

    react_counts <- dplyr::filter(react_counts,
                                  subsystem != 'All reactions')

    react_counts <- dplyr::filter(react_counts,
                                  status %in% c('activated', 'inhibited'))

    reg_counts <- react_counts[c('subsystem', 'status', 'n')]

    all_subs <- sort(unique(reg_counts$subsystem))

    message(paste('Simulation testing for enrichment in n =',
                  length(all_subs), 'subsystems'))

    reg_lst <- plyr::dlply(reg_counts, 'status', dplyr::select, -status)

    reg_lst$regulated <- dplyr::group_by(reg_counts, subsystem)

    reg_lst$regulated <- dplyr::summarise(reg_lst$regulated,
                                          n = sum(n))

    reg_lst <- purrr::map(reg_lst,
                          tibble::column_to_rownames, 'subsystem')

    reg_lst <- purrr::map(reg_lst,
                          ~.x[all_subs, ])

    ## random draws -------

    draws <- draw_reactions(model = model,
                            signif_type = signif_type,
                            regulation_level = regulation_level,
                            n_iter = n_iter,
                            aggregate = TRUE,
                            .parallel = .parallel)

    draws <- purrr::map(draws,
                        tibble::column_to_rownames, 'subsystem')

    draws <- purrr::map(draws,
                        ~.x[all_subs, ])

    ## comparing the observed and random counts ------

    comp_results <- list()

    for(i in names(reg_lst)) {

      comp_results[[i]] <- purrr::map_dfc(draws[[i]],
                                          ~.x >= reg_lst[[i]])

    }

    comp_results <- purrr::map(comp_results, rowSums)

    comp_results <- purrr::map(comp_results, ~.x/n_iter)

    comp_results <- purrr::map(comp_results,
                               ~tibble::tibble(subsystem = all_subs,
                                               p_value = .x))

    comp_results <- purrr::map2_dfr(comp_results,
                                    names(comp_results),
                                    ~dplyr::mutate(.x, status = .y))

    ## output -------

    reg_counts <- plyr::dlply(react_counts,
                              'subsystem',
                              function(x) tibble::tibble(status = 'regulated',
                                                         n = sum(x$n),
                                                         n_total = x$n_total[1]))

    reg_counts <- purrr::map2_dfr(reg_counts,
                                  names(reg_counts),
                                  ~dplyr::mutate(.x, subsystem = .y))

    react_counts <- rbind(react_counts,
                          reg_counts)

    react_counts <- dplyr::mutate(react_counts,
                                  frac_total = n/n_total)

    react_counts <- dplyr::left_join(react_counts,
                                     comp_results,
                                     by = c('subsystem', 'status'))

    react_counts <- dplyr::mutate(react_counts,
                                  p_value = ifelse(p_value == 0,
                                                   1/n_iter, p_value),
                                  p_adjusted = p.adjust(p_value, 'BH'),
                                  subsystem = factor(subsystem, all_subs))

    dplyr::arrange(react_counts, subsystem, status)

  }


# END -------
