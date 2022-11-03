% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{draw_norm}
\alias{draw_norm}
\title{Draw values from the gene regulation estimate distribution.}
\usage{
draw_norm(x, err, n, seed = NULL, .parallel = FALSE)
}
\arguments{
\item{x}{a numeric vector with fold-regulation estimates (not log2!!!).
Elements are named with Entrez IDs.}

\item{err}{a numeric vector with an error estimate such as SD or SEM
(not log2!!!), named with Entrez IDs.}

\item{n}{number of values to be drawn.}

\item{seed}{seed of the random number generator.}

\item{.parallel}{logical, should the function work in parallel? Works only if
TRUE and a parallel backend compatible with furrr::future_map() is provided.}

\item{return}{a list with n numeric values corresponding to the estimates,
named with the Entrez IDs.}
}
\description{
Draws values from the gene regulation estimate
normal distribution given the standard deviation and expected value.
}