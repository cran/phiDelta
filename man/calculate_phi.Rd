% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phi_delta.R
\name{calculate_phi}
\alias{calculate_phi}
\title{calculate phi}
\usage{
calculate_phi(spec, sens, ratio = 1)
}
\arguments{
\item{spec}{is the specificity, the true negative rate}

\item{sens}{is the sensitivity, the true positive rate}

\item{ratio}{is the ratio of positive and negative of the data. The default is 1}
}
\value{
phi
}
\description{
calculates phi out of specificity and sensitivity
depending on the ratio
}
\examples{
calculate_phi(1,0)
calculate_phi(0.5,0.3)
}
\author{
rothe
}
