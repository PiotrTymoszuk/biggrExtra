# Description of the external data sets.

# Example regulation data ------

#' Gene regulation estimates between collagen-subtypes of prostate carcinoma.
#'
#' @description
#' A data frame with gene expression regulation estimates between
#' Collagen^hi^ and Collagen^low^ tumors. The estimates are log~2~
#' fold regulation estimates versus the Collagen^low^ tumor subtype computed by
#' linear regression. The values were calculated for the prostate cancer TCGA
#' (The Cancer Gene Atlas) data set.
#'
#' @format A data frame with 36118 rows and 8 variables.
#' * __gene_symbol__: HGNC gene symbol
#' * __entrez_id__: Entrez ID identifier
#' * _estimate_: log~2~ regulation estimate versus Collagen^low^ subtype
#' * _se_: standard error of the mean for the regulation estimate
#' * __lower_ci__: lower 95%CI interval limit for the regulation estimate
#' * __upper_ci__: upper 95%CI interval limit for the regulation estimate
#' * __p_adjusted__: FDR-corrected p value.
#' * __regulation__: regulation sign.
#'
#' @source TCGA prostate adenocarcinoma data set. Publication pending.
#'
#' @docType data
#'
#' @name tcga_data
#'
#' @usage data(tcga_data)

  NULL

# Reaction annotation data -------

#' Reaction association rules for the Recon2 model.
#'
#' @description
#' Cleared and harmonized reaction annotation data for
#' [Recon2 model of human metabolism](https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1109130000/2/MODEL1109130000_url.xml).
#' The data were processed with the [RECON_processing](https://github.com/PiotrTymoszuk/RECON_processing)
#' R pipeline.
#'
#' @format
#' a data frame with 7440 rows and 10 columns:
#' * __id__ and __metaid__ reaction identifiers
#' * __sboTerm__ Systems Biology Ontology term
#' * __name__ reaction name
#' * __gene_association__ gene - reaction association rules as character strings
#' * __subsystem__ Recon metabolic subsystem
#' * __ec_number__ enzyme's EC identifier
#' * __confidence_level__ confidence level for the gene - reaction association rule
#' * __reference__ literature references
#' * __notes__ additional comments
#'
#' @source
#' [biomodels](https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1109130000/2/MODEL1109130000_url.xml) and
#' [RECON_processing GitHub respository](https://github.com/PiotrTymoszuk/RECON_processing) .
#'
#' @docType data
#'
#' @name Recon2D
#'
#' @md
#' @usage data(Recon2D)

  NULL

#' Reaction association rules for the Recon2.2 model.
#'
#' @description
#' Cleared and harmonized reaction annotation data for
#' [Recon2.2 model of human metabolism](https://www.biomodels.org/biomodels/services/download/get-files/MODEL1603150001/2/MODEL1603150001_url.xml).
#' The data were processed with the [RECON_processing](https://github.com/PiotrTymoszuk/RECON_processing)
#' R pipeline.
#'
#' @format
#' a data frame with 7785 rows and 7 columns:
#' * __id__ and __metaid__ reaction identifiers
#' * __sboTerm__ Systems Biology Ontology term
#' * __name__ reaction name
#' * __gene_association__ gene - reaction association rules as character strings
#' * __subsystem__ Recon metabolic subsystem
#' * __confidence_level__ confidence level for the gene - reaction association rule#'
#'
#' @source
#' [biomodels](https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1109130000/2/MODEL1109130000_url.xml) and
#' [RECON_processing GitHub respository](https://github.com/PiotrTymoszuk/RECON_processing) .
#'
#' @docType data
#'
#' @name Recon2_2D
#'
#' @md
#' @usage data(Recon2_2D)

NULL

# END ------
