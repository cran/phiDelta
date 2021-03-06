% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phi_delta.R
\name{iso_precision}
\alias{iso_precision}
\title{isometric precision lines}
\usage{
iso_precision(ratio = 1, granularity = 0.25, lty = "longdash",
  col = "blue", ...)
}
\arguments{
\item{ratio}{numeric value for the ratio of positive and negative of the data}

\item{granularity}{numeric value between 0 and 1 for the granularity of the lines.
It is a value for the distance between 2 lines}

\item{lty}{the type of line, see \link{par}}

\item{col}{the color of the lines}

\item{...}{further graphical parameters, see \link{par}}
}
\description{
adds isometric lines for the precision to the plot depending on the ratio
}
\examples{
x <- c_statistics(climate_data)
ratio <- calculate_ratio(x)
phiDelta_plot_from_data(x)
iso_precision(ratio, col = "green")
}
\author{
rothe
}
