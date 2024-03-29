% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting-functions.R
\name{moderation_matrix}
\alias{moderation_matrix}
\title{moderation_matrix}
\usage{
moderation_matrix(
  ...,
  effect_size = "Effect size",
  moderators = NULL,
  null_value = NULL,
  transf = NULL,
  leading_zero = TRUE,
  black = NULL,
  replace = c(`_` = " ")
)
}
\arguments{
\item{...}{Named meta_list models. All common moderators should have the same name.}

\item{effect_size}{a string}

\item{moderators}{a list of moderators to include, order retained.}

\item{null_value}{a scalar indicating non-significance (where the dashed line will be drawn).}

\item{transf}{function to transform values}

\item{leading_zero}{when true, leading zeros are allowed on the x-axis}

\item{black}{list of outcome names containing their moderators}

\item{replace}{a vector with names included. gsub will be applied to the moderation column such that the vector's names are replaced with the vector's contents}
}
\description{
Plots a moderation matrix
}
