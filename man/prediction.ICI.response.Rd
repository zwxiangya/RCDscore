% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.ICI.response.R
\name{prediction.ICI.response}
\alias{prediction.ICI.response}
\title{Testing the ICI model in the testing data and the individual independent datasets}
\usage{
prediction.ICI.response(exp.matrix = NULL, parallel.sz = NULL)
}
\arguments{
\item{exp.matrix}{A cohort matrix with the genes as the colnames, sample names as the rownamses. The colnames inlcudes the RCD.Sig}

\item{parallel.sz}{Number of processors to use when doing the calculations in parallel. This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel.sz=0) then it will use all available core processors unless we set this argument with a smaller number. If snow is loaded then we must set this argument to a positive integer number that specifies the number of processors to employ in the parallel calculation.}
}
\value{
A dataframe that contains the predicted response to the ICI
}
\description{
#' This is some description of this function.
}
\examples{
\dontrun{
data("test")
data("RCD.Sig_selected")
dd = prediction.ICI.response(exp.matrix = test,parallel.sz=16)
}
}
