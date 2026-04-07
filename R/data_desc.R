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
#' @source [Heidegger et al. 2025](https://pubmed.ncbi.nlm.nih.gov/38851995/)
#'
#' @references
#' Heidegger I, Frantzi M, Salcher S, Tymoszuk P, Martowicz A, Gomez-Gomez E,
#' Blanca A, Lendinez Cano G, Latosinska A, Mischak H, et al.
#' Prediction of Clinically Significant Prostate Cancer by a Specific Collagen-related
#' Transcriptome, Proteome, and Urinome Signature.
#' Eur Urol Oncol (2024) doi:10.1016/J.EUO.2024.05.014
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
#' [Recon2 model of human metabolism](https://www.ebi.ac.uk/biomodels/MODEL1109130000).
#' The data were processed with the [RECON_processing](https://github.com/PiotrTymoszuk/RECON_processing)
#' R pipeline.
#'
#' @format
#' a data frame with 7440 rows and 10 columns:
#' * __id__ and __metaid__ BiGG reaction identifiers
#' * __sboTerm__ Systems Biology Ontology term
#' * __name__ reaction name
#' * __gene_association__ gene - reaction association rules as character strings
#' * __subsystem__ Recon metabolic subsystem
#' * __ec_number__ EC identifier of enzymatic reaction
#' * __confidence_level__ confidence level for the gene - reaction association rule
#' * __reference__ literature references
#' * __notes__ additional comments
#'
#' @references
#' Thiele I, Swainston N, Fleming RMT, Hoppe A, Sahoo S, Aurich MK,
#' Haraldsdottir H, Mo ML, Rolfsson O, Stobbe MD, et al.
#' A community-driven global reconstruction of human metabolism.
#' Nat Biotechnol 2013 315 (2013) 31:419–425. doi:10.1038/nbt.2488
#'
#' @source
#' [biomodels](https://www.ebi.ac.uk/biomodels/MODEL1109130000) and
#' [RECON_processing GitHub repository](https://github.com/PiotrTymoszuk/RECON_processing) .
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
#' [Recon2.2 model of human metabolism](https://www.ebi.ac.uk/biomodels/MODEL1603150001).
#' The data were processed with the [RECON_processing](https://github.com/PiotrTymoszuk/RECON_processing)
#' R pipeline.
#'
#' @format
#' a data frame with 7785 rows and 7 columns:
#' * __id__ and __metaid__ BiGG reaction identifiers
#' * __sboTerm__ Systems Biology Ontology term
#' * __name__ reaction name
#' * __gene_association__ gene - reaction association rules as character strings
#' * __subsystem__ Recon metabolic subsystem
#' * __confidence_level__ confidence level for the gene - reaction association rule
#'
#' @references
#' Swainston N, Smallbone K, Hefzi H, Dobson PD, Brewer J, Hanscho M, Zielinski DC,
#' Ang KS, Gardiner NJ, Gutierrez JM, et al.
#' Recon 2.2: from reconstruction to model of human metabolism.
#' Metabolomics (2016) 12. doi:10.1007/S11306-016-1051-4
#'
#' @source
#' [biomodels](https://www.ebi.ac.uk/biomodels/MODEL1603150001) and
#' [RECON_processing GitHub repository](https://github.com/PiotrTymoszuk/RECON_processing) .
#'
#' @docType data
#'
#' @name Recon2_2D
#'
#' @md
#' @usage data(Recon2_2D)

  NULL

#' Reaction association rules for the Human-GEM 2.0.0 model.
#'
#' @description
#' Cleared and harmonized reaction annotation data for
#' [Human-GEM 2.0.0 model of metabolism](https://github.com/SysBioChalmers/Human-GEM).
#' The data were processed with the [RECON_processing](https://github.com/PiotrTymoszuk/RECON_processing)
#' R pipeline.
#'
#' @format
#' a data frame with 12931 rows and 8 columns:
#' * __id__ Metabolic Atlas reaction identifiers
#' * __name__ reaction name
#' * __ec_number__ EC identifier of enzymatic reactions
#' * __gene_association__ gene - reaction association rules as character strings
#' * __subsystem__ Metabolic Atlas subsystem
#' * __confidence_level__ confidence level for the gene - reaction association rule
#' * __reference__ literature references
#' * __miriam__ cross-references to other databases
#'
#' @references
#' Robinson JL, Kocabaş P, Wang H, Cholley PE, Cook D, Nilsson A, Anton M,
#' Ferreira R, Domenzain I, Billa V, et al.
#' An atlas of human metabolism. Sci Signal (2020) 13:eaaz1482–eaaz1482.
#' doi:10.1126/SCISIGNAL.AAZ1482
#'
#'
#' @source
#' [GitHub repository of Metabolic Atlas](https://github.com/SysBioChalmers/Human-GEM) and
#' [RECON_processing GitHub repository](https://github.com/PiotrTymoszuk/RECON_processing)
#'
#' @docType data
#'
#' @name Human_GEM_2_0_0
#'
#' @md
#' @usage data(Human_GEM_2_0_0)

  NULL

# END ------
