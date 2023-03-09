# Modeling functions

# Build SBML model from gene IDs ------

#' Create a SBML model from gene Entrez IDs.
#'
#' @description Builds a SBML model given gene Entrez ID identifiers.
#' Technically, an up-stream function of \code{\link{get_regulation}}.
#' @details see: \code{\link{get_regulation}}.
#' @return a geneSBML object containing the following:
#' * a rsbml Model built with the reaction IDs associated with the provided genes,
#' as described for \code{\link[BiGGR]{buildSBMLFromReactionIDs}}.
#' * a data frame with reaction regulation estimates and their errors as
#' specified for \code{\link{get_regulation}}.
#' @inheritParams get_regulation
#' @export

  build_geneSBML <- function(x,
                             err = NULL,
                             database,
                             or_fun = c('mean', 'median', 'min', 'max'),
                             and_fun = c('min', 'max', 'mean', 'median'),
                             x_default = 1,
                             err_method = c('norm', 'mc'),
                             n_iter = 1000,
                             ci_method = c('perc', 'bca'),
                             .parallel = TRUE) {

    ## entry control is accomplished by downstream functions
    ## reaction regulation estimates

    res <- get_regulation(x = x,
                          err = err,
                          database = database,
                          or_fun = or_fun,
                          and_fun = and_fun,
                          x_default = x_default,
                          err_method = err_method,
                          n_iter = n_iter,
                          .parallel = .parallel,
                          return_genes = TRUE,
                          return_mc = TRUE)

    reg <- dplyr::filter(res$reg, !is.na(fold_reg))

    gene_map <- dplyr::filter(res$gene_map,
                              react_id %in% reg$react_id)

    sub_map <- extract_subsystems(database)

    reg <- dplyr::left_join(reg,
                            sub_map,
                            by = 'react_id')

    geneSBML(model = buildSBMLFromReactionIDs(reaction.ids = reg$react_id,
                                              database = database),
             reg = reg,
             gene_map = gene_map,
             mc = res$mc)

  }

# END -----
