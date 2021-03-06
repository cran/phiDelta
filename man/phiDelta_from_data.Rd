% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phi_delta.R
\name{phiDelta_from_data}
\alias{phiDelta_from_data}
\title{phi delta matrix}
\usage{
phiDelta_from_data(stats, ratio_corrected = TRUE)
}
\arguments{
\item{stats}{c_statistics}

\item{ratio_corrected}{locigal, if true phi and delta will be calculated in respect to the ratio
of positive and negative samples}
}
\value{
dataframe, first column are the names of the features
second column the phi values
third column the delta values
}
\description{
calculates phi and delta directly from the stats
and generates a matrix with the names of the features, their phi
and their delta value
}
\examples{
x <- c_statistics(climate_data)
phiDelta <- phiDelta_from_data(x, ratio_corrected = FALSE)
with_ratio <- phiDelta_from_data(x)
}
\author{
rothe
}
