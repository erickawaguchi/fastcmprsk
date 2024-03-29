% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.fcrr.R
\name{summary.fcrr}
\alias{summary.fcrr}
\title{Summary method for fastCrr}
\usage{
\method{summary}{fcrr}(
  object,
  conf.int = TRUE,
  alpha = 0.05,
  digits = max(options()$digits - 5, 2),
  ...
)
}
\arguments{
\item{object}{\code{fcrr} x (output from \code{fastCrr()})}

\item{conf.int}{Logical. Whether or not to outut confidence intervals.}

\item{alpha}{Significance level of the confidence intervals.}

\item{digits}{Numer of significant difits to round to.}

\item{...}{additional arguments to \code{print()}}
}
\value{
The form of the value returned by \code{summary} depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.
}
\description{
Generate and print summaries of \code{fastCrr} output.
}
\details{
The summary method produces an ANOVA table for the coefficient estimates of the Fine-Gray model.
}
