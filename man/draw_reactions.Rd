% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_utils.R
\name{draw_reactions}
\alias{draw_reactions}
\title{Random reaction draws.}
\usage{
draw_reactions(
  model,
  signif_type = c("fdr", "raw"),
  regulation_level = 1,
  n_iter = 1000,
  .parallel = FALSE
)
}
\arguments{
\item{model}{a geneSBML model.}

\item{signif_type}{type of significance to be used for definition of
regulated reactions:
'raw' (default) or 'fdr' (FDR-corrected).
Ignored if the model contains no error information.}

\item{regulation_level}{fold-regulation level cutoff used for definition
of regulated reactions.}

\item{n_iter}{number of random draws from the reaction set.}

\item{.parallel}{logical, should the computation be run in parallel?}
}
\value{
a list of matrices with counts of activated, inhibited
and regulated reactions in each subsystem.
}
\description{
Draws randomly i, j and k reactions from all reactions
available in the model, where i, j and k correspond to the numbers
of activated, inhibited and regulated reactions.
}
