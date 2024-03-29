% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/publication_bias.r
\name{FPP}
\alias{FPP}
\title{FPP}
\usage{
FPP(m, transf = function(x) x, round = 2, alpha = 0.05)
}
\arguments{
\item{m}{a metaSEM or meta_list model}

\item{transf}{function with which to transform results}

\item{round}{number of digits for rounding}

\item{alpha}{the alpha threshold used in PET-PEESE. Defaults to 0.1 \(10\% as per Stanley 2017\)}
}
\description{
Calculates a FAT-PET-PEESEE adjusted estimate
}
\details{
PET-PEESE uses meta-regression in order to adjust estimates for publication bias. When there is publication bias the sampling variance becomes correlated with effect size. By using the standard error (PET) or the sampling variance (PEESE) as moderators in a meta-regression, estimates (the intercept) can be made which partial out the correlation between effect size and variance. PET-PEESE first tests whether the intercept is a significant when controlling for the standard error of effects (p < 0.05) predictor of effect size. If it is, PEESE is used. Otherwise, PET is used.
}
