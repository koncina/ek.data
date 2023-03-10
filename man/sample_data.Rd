% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samples.R
\name{sample_data}
\alias{sample_data}
\title{Extract available samples as a tibble
A wrapper around colData() for SummarizedExperiment
and pData() for ExpressionSet objects
returning a tibble}
\usage{
sample_data(x, ...)
}
\arguments{
\item{x}{A SummarizedExperiment or ExpressionSet object.}

\item{...}{Arguments to be passed to methods}
}
\description{
Extract available samples as a tibble
A wrapper around colData() for SummarizedExperiment
and pData() for ExpressionSet objects
returning a tibble
}
