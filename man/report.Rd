% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reports.R
\name{report}
\alias{report}
\title{report}
\usage{
report(meta_list, ..., rmarkdown = FALSE, digits = 2, transf = function(x) x)
}
\arguments{
\item{meta_list}{the meta_list object}

\item{...}{things to report. One o}

\item{rmarkdown}{return results in rmarkdown?}

\item{digits}{the number of digits to return}

\item{transf}{you can supply a function to transform baseline pooled estimates}
}
\description{
Constructs written reports about the contents of models
}
