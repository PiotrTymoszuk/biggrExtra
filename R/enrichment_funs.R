# Functions for testing enrichment of significantly regulated reactions
# in metabolic subsystems

# Enrichment testing -------

#' Test for enrichment of significantly regulated reactions
#' in metabolic subsystems.
#'
#' @description Tests for enrichment of significantly regulated reactions
#' within metabolic subsystems. This is done with one of two methods:
#' 'fisher' using Fisher's exact test or 'simulation' where the observed
#' enrichment is compared with random draws from the reaction set.
#' @details Technically, Fisher's tests are done with
#' \code{\link[rstatix]{fisher_test}}. This method is much faster but biased
#' towards subsystems with multiple reactions mapped to them
#' (like fatty acid oxidation).
#' @return a data frame with subsystem name and testing results
#' for regulated, activated and inhibited reactions
#' as well as FDR-adjusted p values.
#' @param model a geneSBML model.
#' @param signif_type type of significance to be used for definition of
#' regulated reactions:
#' 'raw' (default) or 'fdr' (FDR-corrected).
#' Ignored if the model contains no error information.
#' @param regulation_level fold-regulation level cutoff used for definition
#' of regulated reactions.
#' @param method 'fisher' or 'simulation' as described above.
#' Defaults to 'fisher'.
#' @param alternative indicates the alternative hypothesis
#' and must be one of "two.sided", "greater" or "less". Defaults to "greater".
#' Applies only to the 'fisher' method.
#' @param n_iter number of random draws from the reaction set.
#' @param .parallel logical, should the computation be run in parallel?
#' Defaults to FALSE.
#' @param ... extra arguments passed to \code{\link[rstatix]{fisher_test}}.
#' @export

  suba <- function(model,
                   signif_type = c('raw', 'fdr'),
                   regulation_level = 1,
                   method = c('fisher', 'simulation'),
                   alternative = 'greater',
                   n_iter = 1000,
                   .parallel = FALSE, ...) {

    ## entry control ----------

    if(!is_geneSBML(model)) {

      stop('A valid geneSBML model is required!', call. = FALSE)

    }

    signif_type <- match.arg(signif_type[[1]],
                             c('raw', 'fdr'))

    method <- match.arg(method[[1]],
                        c('fisher', 'simulation'))

    stopifnot(is.logical(.parallel))
    stopifnot(is.numeric(n_iter))

    n_iter <- as.integer(n_iter)

    ## subsystem counts ------

    sub_counts <- count(model,
                        signif_type = signif_type,
                        regulation_level = regulation_level)

    sub_counts <- dplyr::filter(sub_counts,
                                subsystem != 'All reactions')

    ## testing ------

    switch(method,
           fisher = sub_fisher(data = sub_counts,
                               alternative = alternative,
                               .parallel = .parallel, ...),
           simulation = sub_simulation(model = model,
                                       signif_type = signif_type,
                                       regulation_level = regulation_level,
                                       n_iter = n_iter,
                                       .parallel = .parallel))

  }

# END -------
