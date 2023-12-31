% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Biodiv_functions.R
\name{plotDiv}
\alias{plotDiv}
\title{plotDiv(data, tog, cap) function to plot the value of diversity}
\usage{
plotDiv(data, tog = TRUE, cap = FALSE)
}
\arguments{
\item{data}{This parameter must necessarily be an object containing the outcome of the Div function to prepare the graph.}

\item{tog}{This parameter is an abbreviation of the word "together" and is a logical object. If it is set to TRUE, the function will plot the three graphs together (side by side).}

\item{cap}{This parameter is an abbreviation of the word "captions" and is a logical object. If it is set to TRUE, the function will ask you to set the labels for the graph.}
}
\value{
Three graphs, namely (functional diversity, taxonomic diversity, and functional redundancy).
}
\description{
This is a function to plot the results of the FDchao function
filtered by the Div function, generating three graphs:
Taxonomic Diversity, Functional Diversity, and Functional Redundancy.
}
