% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scIntegration.v113.allow_fixed_B.R
\name{plot_vectorField}
\alias{plot_vectorField}
\title{Visualization of the vector field embedding}
\usage{
plot_vectorField(
  embedding_mat,
  colour = 1,
  alpha = 1,
  randomize = T,
  nbin = 50,
  minNum = 10
)
}
\arguments{
\item{embedding_mat}{Embedding from the vector field}

\item{colour}{vector of variable to plot}

\item{alpha}{alpha}

\item{randomize}{Logical. Whether to randomize data before plotting.}

\item{nbin}{number of bins}

\item{minNum}{minNum}
}
\value{
ggplot2 object
}
\description{
Visualization of the vector field embedding
}
