% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scIntegration.v113.allow_fixed_B.R
\name{getDiffQ.gedi}
\alias{getDiffQ.gedi}
\title{Sample-variable effects on Qi}
\usage{
getDiffQ.gedi(object, contrast)
}
\arguments{
\item{object}{GEDI object}

\item{contrast}{a vector representing the sample-variable contrast to use for extracting effect on Qi; must be the same size as aux$L}
}
\value{
Matrix R
}
\description{
Return the sample-variable effects on Qi
}
