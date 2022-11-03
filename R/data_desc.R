# Description of the external data sets.

# Example regulation data ------

#' Gene regulation estimates between collagen-subtypes of prostate carcinoma.
#'
#' @description A data frame with gene expression regulation estimates between
#' Collagen^hi^ and Collagen^low^ tumors. The estimates are log~2~
#' fold regulation estimates versus the Collagen^low^ tumor subtype computed by
#' linear regression. The values were calculated for the prostate cancer TCGA
#' (The Cancer Gene Atlas) data set.
#' @format A data frame with 36118 rows and 8 variables.
#' * __gene_symbol__: HGNC gene symbol
#' * __entrez_id__: Entrez ID identifier
#' * _estimate_: log~2~ regulation estimate versus Collagen^low^ subtype
#' * _se_: standard error of the mean for the regulation estimate
#' * __lower_ci__: lower 95%CI interval limit for the regulation estimate
#' * __upper_ci__: upper 95%CI interval limit for the regulation estimate
#' * __p_adjusted__: FDR-corrected p value.
#' * __regulation__: regulation sign.
#' @source TCGA prostate adenocarcinoma data set. Publication pending.
#' @docType data
#' @name tcga_data
#' @usage data(tcga_data)
NULL

# Reaction annotation data -------

#' BiGG reaction annotation data.
#'
#' @description A data frame with BiGG database annotation of metabolic
#' reactions.
#' @format A data frame with 28301 rows and 6 variables:
#' * __bigg_id__: BiGG reaction identifier
#' * __name__: reaction name
#' * __reaction_string__: reaction equation
#' * __model_list__: list of BiGG models containing the reaction
#' * __database_links__: links to the BiGG on-line database
#' * __old_bigg_ids__: legacy BiGG ID
#' @source BiGG: http://bigg.ucsd.edu/data_access
#' @docType data
#' @name reactions
#' @usage data(reactions)
NULL

# Metabolite annotation data ------

#' BiGG metabolite annotation data.
#'
#' @description A data frame with BiGG database annotation of metabolites.
#' @format A data frame with 15724 rows and 6 variables:
#' * __bigg_id__: BiGG metabolite identifier
#' * __universal_bigg_id__: universal BiGG metabolite identifier
#' * __name__: metabolite name
#' * __model_list__: list of models with containing the metabolite
#' * __database_links__: links to the BiGG on-line database
#' * __old_bigg_ids__: legacy BiGG ID
#' @source BiGG: http://bigg.ucsd.edu/data_access
#' @docType data
#' @name metabolites
#' @usage data(metabolites)
NULL