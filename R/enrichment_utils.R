# Utility functions (no-exported) for testing of significant
# enrichment of significantly regulated reactions within metabolic subsystems.

# Fisher's exact test ------

#' Test for enrichment in a subsystem with Fisher's exact test.
#'
#' @description
#' Tests for significant enrichment of generally significantly
#' regulated, significantly activated and significantly inhibited reactions in
#' a metabolic subsystem.
#' The significance is determined by Fisher's exact test (H0: no significantly
#' regulated reactions)
#'
#' @details The method is biased towards subsystems with multiple reactions!
#' Technically, Fisher's exact tests are computed with
#' \code{\link[rstatix]{fisher_test}}.
#'
#' @param data a data frame with the following variables:
#' 'subsystem', 'status', 'n' (number of regulated reactions) and
#' 'n_total' (total number of reactions in a subsystem).
#' @param alternative indicates the alternative hypothesis
#' and must be one of "two.sided", "greater" or "less". Defaults to "greater".
#' @param ... extra arguments passed to \code{\link[rstatix]{fisher_test}}.
#'
#' @return a data frame with subsystem name and testing results
#' for regulated, activated and inhibited reactions
#' as well as, for, `sub_fisher()` FDR-adjusted p values.

  sub_fisher_ <- function(data, alternative = 'greater', ...) {

    ## entry control done be the upstream function!

    ## analysis data ------

    ### in case there are multiple subsystems in 'data':

    subsystem <- NULL

    data <- filter(data, subsystem == data$subsystem[[1]])

    data <- rbind(data,
                  tibble(subsystem = data$subsystem[[1]],
                         status = 'regulated',
                         n = sum(data$n),
                         n_total = data$n_total[[1]]))

    frac_total <- NULL
    n <- NULL
    n_total <- NULL

    data <- mutate(data, frac_total = n/n_total)

    ## test matrices ------

    tst_mtx <- split(data[c('n', 'n_total')],
                     f = data[['status']])

    tst_mtx <-
      map(tst_mtx,
          function(x) rbind(c(activated = x$n[1],
                              constant = x$n_total[1] - x$n[1]),
                            c(activated = 0,
                              constant = x$n_total[1])))

    for(i in names(tst_mtx)) {

      rownames(tst_mtx[[i]]) <- c('model', 'expected')

    }

    ## Fisher tests -------

    test_res <- map(tst_mtx,
                    ~rstatix::fisher_test(xtab = as.table(.x),
                                          alternative = alternative,
                                          ...))

    status <- NULL

    test_res <-
      map2_dfr(test_res,
               names(test_res),
               ~mutate(.x, status = .y))

    ## output --------

    p_value <- NULL
    n <- NULL
    p <- NULL
    p.signif <- NULL

    test_res <- mutate(test_res, p_value = p)

    test_res <- select(test_res, -n, -p, -p.signif)

    left_join(data, test_res, by = 'status')

  }

#' @rdname sub_fisher_

  sub_fisher <- function(data,
                         alternative = 'greater', ...) {

    ## entry control ------

    start_time <- Sys.time()
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)),
            add = TRUE)

    unique_subsystems <- unique(data$subsystem)

    no_subsystems <- length(unique_subsystems)

    message(paste('Fisher testing for enrichment in n =',
                  no_subsystems, 'subsystems'))

    alternative <- match.arg(alternative[[1]],
                             c('greater', 'less', 'two.sided'))

    if(!is.data.frame(data)) {

      stop(paste("'data' hast to be a data frame with the 'subsystem',",
                 "'status', 'n' and 'n_total' variables."),
           call. = FALSE)

    }

    if(!all(c('subsystem', 'status', 'n', 'n_total') %in% names(data))) {

      stop(paste("'data' hast to be a data frame with the 'subsystem',",
                 "'status', 'n' and 'n_total' variables."),
           call. = FALSE)

    }

    data <- filter(data, status %in% c('activated', 'inhibited'))

    if(any(!unique(data$status) %in% c('activated', 'inhibited'))) {

      stop(paste("The 'status' variable of 'data' has to have 'activated'",
                 "and 'inhibited' levels."),
           call. = FALSE)

    }

    ## serial testing -------

    subsystem <- NULL
    status <- NULL
    n <- NULL
    n_total <- NULL
    frac_total <- NULL
    p_value <- NULL
    p_adjusted <- NULL

    data <-
      mutate(data,
             subsystem = factor(subsystem, unique_subsystems),
             status = droplevels(status))

    sub_lst <- split(data, data[['subsystem']])

    tst_results <-
      map(sub_lst,
          sub_fisher_,
          alternative = alternative, ...)

    tst_results <- do.call('rbind', tst_results)

    mutate(as_tibble(tst_results),
           frac_total = n/n_total,
           p_adjusted = p.adjust(p_value, 'BH'))

  }

# Random reaction draws ------

#' Random reaction draws.
#'
#' @description
#' Draws randomly i, j and k reactions from all reactions
#' available in the model, where i, j and k correspond to the numbers
#' of activated, inhibited and regulated reactions.
#'
#' @return
#' a list of matrices with counts of activated, inhibited
#' and regulated reactions in each subsystem.
#'
#' @param model a geneSBML model.
#' @param signif_type type of significance to be used for definition of
#' regulated reactions:
#' 'raw' (default) or 'fdr' (FDR-corrected).
#' Ignored if the model contains no error information.
#' @param regulation_level fold-regulation level cutoff used for definition
#' of regulated reactions.
#' @param n_iter number of random draws from the reaction set.
#' @param .parallel logical, should the computation be run in parallel?


  draw_reactions <- function(model,
                             signif_type = c('fdr', 'raw'),
                             regulation_level = 1,
                             n_iter = 1000,
                             .parallel = FALSE) {

    ## entry control: done by the up-level function ------

    n <- NULL
    n_total <- NULL
    p_value <- NULL
    subsystem <- NULL

    on.exit(plan('sequential'))

    ## reaction vector and counts ------

    all_subs <- components(model, 'regulation')[['subsystem']]

    subsystem <- NULL
    status <- NULL
    n <- NULL

    reg_counts <- count(model,
                        signif_type = signif_type,
                        regulation_level = regulation_level)

    reg_counts <- filter(reg_counts, subsystem == 'All reactions')

    reg_counts <- filter(reg_counts, status %in% c('activated', 'inhibited'))

    reg_counts <- reg_counts[c('status', 'n')]

    total_counts <- tibble(status = 'regulated',
                           n = sum(reg_counts$n))

    reg_counts <- rbind(reg_counts, total_counts)

    reg_counts <- mutate(reg_counts, status = droplevels(status))

    reg_counts <- split(reg_counts$n, reg_counts$status)

    ## random draws ------

    draw_names <- paste0('draw_', 1:n_iter)

    if(.parallel) plan('multisession')

    draws <-
      future_map(reg_counts,
                 ~vecDraws(all_subs,
                           size = .x,
                           n_iter = n_iter),
                 .options = furrr_options(seed = TRUE))

    return(draws)

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
#' @param .parallel logical, should the computation be run in parallel?

  sub_simulation <- function(model,
                             signif_type = c('fdr', 'raw'),
                             regulation_level = 1,
                             n_iter = 1000,
                             .parallel = FALSE) {

    ## entry control -------

    subsystem <- NULL
    status <- NULL
    n <- NULL
    n_total <- NULL
    p_value <- NULL

    start_time <- Sys.time()
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

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

    react_counts <- filter(react_counts,
                           subsystem != 'All reactions')

    react_counts <- filter(react_counts,
                           status %in% c('activated', 'inhibited'))

    react_counts[['status']] <- droplevels(react_counts[['status']])

    reg_counts <- react_counts[c('subsystem', 'status', 'n')]

    all_subs <- sort(unique(reg_counts$subsystem))

    message(paste('Simulation testing for enrichment in n =',
                  length(all_subs), 'subsystems'))

    reg_lst <- split(reg_counts[c('subsystem', 'n')],
                     f = reg_counts[['status']])

    reg_lst$regulated <- group_by(reg_counts, subsystem)

    reg_lst$regulated <- summarise(reg_lst$regulated, n = sum(n))

    reg_lst <- map(reg_lst, column_to_rownames, 'subsystem')

    reg_lst <- map(reg_lst, ~.x[all_subs, ])

    ## random draws -------

    draws <- draw_reactions(model = model,
                            signif_type = signif_type,
                            regulation_level = regulation_level,
                            n_iter = n_iter,
                            .parallel = .parallel)

    draws <- map(draws, ~.x[all_subs, ])

    ## comparing the observed and random counts ------

    comp_results <- list()

    for(i in names(reg_lst)) {

      comp_results[[i]] <- compareCounts(draws[[i]], reg_lst[[i]])

      comp_results[[i]] <- comp_results[[i]][, -1]

      colnames(comp_results[[i]]) <-
        c('mean_n_expected', 'OR', 'n_failed_draws', 'p_value')

      comp_results[[i]] <- mutate(as.data.frame(comp_results[[i]]),
                                  status = i,
                                  subsystem = names(reg_lst[[i]]))

      comp_results[[i]] <-
        rownames_to_column(comp_results[[i]], 'subsystem')

    }

    ## output -------

    reg_counts <- split(react_counts, react_counts[['subsystem']])

    reg_counts <-
      map(reg_counts,
          ~tibble(status = 'regulated',
                  n = sum(.x$n),
                  n_total = .x$n_total[1]))

    reg_counts <-
      map2_dfr(reg_counts, names(reg_counts),
               ~mutate(.x, subsystem = .y))

    react_counts <- rbind(react_counts, reg_counts)

    react_counts <- mutate(react_counts,
                           frac_total = n/n_total)

    react_counts <- left_join(react_counts,
                              do.call('rbind', comp_results),
                              by = c('subsystem', 'status'))

    react_counts <- mutate(react_counts,
                           p_value = ifelse(p_value == 0,
                                            1/n_iter, p_value),
                           p_adjusted = p.adjust(p_value, 'BH'),
                           subsystem = factor(subsystem, all_subs))

    arrange(react_counts, subsystem, status)

  }

# END -------
