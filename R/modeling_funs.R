# Modeling functions

# Build SBML model from gene IDs ------

#' Create a SBML model from gene Entrez IDs.
#'
#' @description
#' Computes regulation estimates for the reactions present
#' in the database based on the regulation estimates for gene expression.
#' As described in the seminal BiGGR package paper
#' (DOI 10.1371/journal.pone.0119016), it is assumed that the change
#' in gene expression corresponds directly to the magnitude of modulation of
#' the pathway, e.g. two-fold increase in expression of an enzyme leads to a
#' two-fold increase in activity of the respective pathway
#' (gene - protein - reaction or GPR principle). In case, more genes
#' are mapped to a given pathway, the fold regulation is computed accordingly
#' to the logical operators between the genes.
#' If a given reaction is not associated with any genes, its regulation estimate
#' is set as specified by the `x_default` argument.
#' If a vector with regulation estimate errors is provided, the cumulative error
#' will be returned as well - see Details.
#'
#' @details
#' Works in a comparable fashion to \code{\link[BiGGR]{gprMapping}} but
#' considerably faster.
#' When a vector of errors for gene regulation is provided, there are two ways
#' to compute estimates and errors for the pathway regulation:
#'
#' * __from the normal distribution__: the regulation estimates are computed
#' by evaluation of the gene - reaction association rules directly with
#' expression regulation estimates provided in `x`. The cumulative reaction
#' error is calculated as a sum of errors for the genes mapped to the reaction,
#' z statistic is computed as (fold-regulation - 1)/error,
#' raw p values and 95% confidence intervals are obtained from the theoretical
#' distribution of z statistics and corrected for multiple testing with the
#' FDR/Benjamini-Hochberg method.
#' While this method is computationally fast, the error, t statistic,
#' confidence intervals and p values are just raw estimates of significance.
#'
#' * __by Monte Carlo method__: given the gene expression regulation estimates
#' in `x` and their errors provided in `err`, a bunch of random values is drawn
#' from the normal distribution of with means specified by `x` and standard
#' deviations specified by `err`. Subsequently, activity regulation estimates
#' of the reactions are computed by evaluation of the gene - reaction
#' association rules for each set of simulated expression estimates. The
#' expected values of reaction activity regulation activity are derived as means
#' or medians across the algorithm iteration.
#'
#' @references
#' Gavai AK, Supandi F, Hettling H, Murrell P, Leunissen JAM, Van Beek JHGM.
#' Using Bioconductor Package BiGGR for Metabolic Flux Estimation Based on
#' Gene Expression Changes in Brain. PLoS One (2015) 10:e0119016.
#' doi:10.1371/JOURNAL.PONE.0119016
#'
#' @references
#' King ZA, Lu J, Dräger A, Miller P, Federowicz S, Lerman JA, Ebrahim A,
#' Palsson BO, Lewis NE. BiGG Models: A platform for integrating,
#' standardizing and sharing genome-scale models.
#' Nucleic Acids Res (2016) 44:D515–D522. doi:10.1093/NAR/GKV1049
#'
#' @return a geneSBML object containing the following:
#'
#' * `model`: a rsbml Model built with the reaction IDs associated with the
#' provided genes, as described for
#' \code{\link[BiGGR]{buildSBMLFromReactionIDs}}.
#'
#' * `reg`: a data frame with reaction regulation estimates and their errors as
#' specified for \code{\link{get_regulation}}.
#'
#' * `gene_map`: a data frame with mapping of the reaction to gene identifiers
#' and association rules
#'
#' * `mc`: for models with regulation estimates and errors computed with the
#' Monte Carlo method, this component stores a matrix with results of
#' the simulation, i.e. activity regulation estimates for particular reactions
#' at algorithm iterations. Omitted when `save_memory = TRUE`.
#'
#' @inheritParams get_regulation
#' @param save_memory logical, should an abridged object without results of
#' single iterations of the Monte Carlo algorithm and without the SBML model
#' be returned?
#'
#' @export

  build_geneSBML <- function(x,
                             err = NULL,
                             database = biggrExtra::Recon2D,
                             scale = c('identity', 'log2', 'log'),
                             or_fun = c('mean', 'median', 'min', 'max'),
                             and_fun = c('min', 'max', 'mean', 'median'),
                             x_default = 1,
                             err_method = c('norm', 'mc'),
                             n_iter = 1000,
                             mc_estimate = c('mean', 'median'),
                             burn_in = 0,
                             ci_method = c('perc', 'bca'),
                             seed = NULL,
                             save_memory = FALSE,
                             .parallel = TRUE) {

    ## entry control is accomplished by downstream functions
    ## reaction regulation estimates

    stopifnot(is.logical(save_memory))

    res <- get_regulation(x = x,
                          err = err,
                          database = database,
                          scale = scale,
                          or_fun = or_fun,
                          and_fun = and_fun,
                          x_default = x_default,
                          err_method = err_method,
                          n_iter = n_iter,
                          mc_estimate = mc_estimate,
                          burn_in = burn_in,
                          ci_method = ci_method,
                          .parallel = .parallel,
                          return_genes = TRUE,
                          return_mc = !save_memory)

    fold_reg <- NULL
    react_id <- NULL

    reg <- filter(res$reg, !is.na(fold_reg))

    gene_map <- filter(res$gene_map, react_id %in% reg$react_id)

    sub_map <- extract_subsystems(database)

    reg <- left_join(reg,
                     sub_map,
                     by = 'react_id')

    if(!save_memory) {

      model <- buildSBMLFromReactionIDs(reaction.ids = reg$react_id,
                                        database = database)

      mc <- res$mc

    } else {

      model <- NULL
      mc <- NULL

    }

    geneSBML(model = model,
             reg = reg,
             gene_map = gene_map,
             mc = mc)

  }

# END -----
