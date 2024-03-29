% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Biodiv_functions.R
\name{FDchao}
\alias{FDchao}
\title{FDchao(data, distance, tau, q, boot) Functional Diversity of N sites for various values of tau and q}
\usage{
FDchao(data, distance, tau, q, boot)
}
\arguments{
\item{data}{a list with N sites; each element of list is species abundances.}

\item{distance}{a matrix of species-pairwise distances.}

\item{tau}{a numeric or a vector of levels of threshold distinctiveness.}

\item{q}{a numeric or a vector of diversity orders; the suggested range for q is \link{0, 2}.}

\item{boot}{a numeric of number of bootstrap replications.}
}
\value{
two matrices of FD; the first matrix is the q-profile information, the second matrix is the tau profile information.
}
\description{
This function calculate Functional Diversity of N sites
for various values of tau and q. For more details on the values of each parameter, refer to the examples of the Div function.
}
\examples{
#minimum tau
min(Bdist[Bdist>0])
#maximum tau
max(Bdist)
#mean tau
library(magrittr)
Babund \%>\%
lapply(FUN = function(x) x/sum(x)) \%>\%
 do.call(cbind, .) \%>\%
 apply(1, mean) \%>\%
 {sum((.\%*\%t(.))*Bdist)}
data <- FDchao(Babund, Bdist, seq(min(Bdist[Bdist>0]), max(Bdist), length.out = 25), seq(from = 0, to = 2, length.out = 25), 50)
}
