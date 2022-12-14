% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_regulation}
\alias{get_regulation}
\title{Calculate reaction regulation estimated based on gene expression.}
\usage{
get_regulation(
  x,
  err = NULL,
  database,
  or_fun = c("mean", "median", "min", "max"),
  and_fun = c("min", "max", "mean", "median"),
  x_default = 1,
  return_genes = FALSE,
  return_mc = FALSE,
  err_method = c("norm", "mc"),
  n_iter = 1000,
  ci_method = c("perc", "bca"),
  seed = NULL,
  .parallel = TRUE
)
}
\arguments{
\item{x}{a numeric vector with fold-regulation estimates (not log2!!!).
Elements are named with Entrez IDs.}

\item{err}{a numeric vector with an error estimate such as SD or SEM
(not log2!!!), named with Entrez IDs.}

\item{database}{an object of class SBMLDocument.}

\item{or_fun}{one of 'mean' (default), 'median', 'min' or 'max'. Specifies the name of
the function used to handle the OR operator between the genes.}

\item{and_fun}{one of 'mean', 'median', 'min' (default) or 'max'. Specifies the name of
the function used to handle the OR operator between the genes.}

\item{x_default}{the default numeric value or NA for regulation of the genes
present in the database and absent from the regulation data. If NA,
the genes missing from x won't be included in subsequent SBML model creation.}

\item{return_genes}{logical, should a data frame with assignment of the
genes to reactions be returned? Defaults to FALSE.}

\item{return_mc}{logical, should all Monte Carlo-determined regulation
estimates for reactions be returned? Defaults to FALSE.}

\item{err_method}{specifies the method of calculation of errors for reactions
as described in Details: 'norm' (from normal distribution, default) or
'mc' (Monte Carlo simulation).}

\item{n_iter}{number of iterations for the Monte carlo simulation of
the reaction errors.}

\item{ci_method}{method of calculation of confidence intervals:
percentile ('perc', default) or BCA ('bca').}

\item{seed}{a seed for the rangom number generator.}
}
\value{
a data frame with the following variables: reaction identifiers
and the fold regulation estimates. If err is not NULL, an error estimate per
reaction is included in the data frame as well along with the z statistic
(fold-regulation/error), 95\% confidence intervals and
raw p-values (normal distribution) and FDR-corrected p-values.
If return_genes is set to TRUE, a gene-mapping table is returned as well in
a list with the fold regulation data frame.
}
\description{
Computes regulation estimates for the reactions present
in the database based on the regulation estimates for gene expression.
As described in the seminal BiGGR package paper
(DOI 10.1371/journal.pone.0119016), it is assumed that the change
in gene expression corresponds directly to the magnitude of modulation of
the pathway, e.g. two-fold increase in expression of an enzyme leads to a
two-fold increase in activity of the respective pathway
(gene - protein - reaction or GPR principle). In case, more genes
are mapped to a given pathway, the fold regulation is computed accordingly
to the logical operators between the genes.
If a given reaction is not associate with any genes, its regulation estimate
is set as specified by the x_default argument.
If a vector with regulation estimate errors is provided, the cumulative error
will be returned as well - see Details.
}
\details{
Works in a comparable fashion to \code{\link[BiGGR]{gprMapping}}
but considerably faster.
When a vector of errors for gene regulation is provided, there are two ways
to compute errors for the pathway regulation:
\itemize{
\item \strong{from the normal distribution}: the cumulative reaction error is
calculated as a sum of errors for the genes mapped to the reaction,
z statistic is computed as (fold-regulation - 1)/error,
raw p values and 95\% confidence intervals are obtained from the theoretical
distribution of z statistics and corrected for multiple testing with the
FDR/Benjamini-Hochberg method.
While this method is computationally fast, the error, t statistic,
confidence intervals and p values are just raw estimates of significance.
\item \strong{by Monte Carlo method}: this option is still pending.
}
}
\references{
Gavai et al. \url{https://doi.org/10.1371/journal.pone.0119016}.
}
