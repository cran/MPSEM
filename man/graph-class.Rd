% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-class.R
\docType{class}
\name{graph-class}
\alias{graph-class}
\alias{print.graph}
\title{Class and Method for Directed Graphs}
\format{
A \code{graph-class} object contains:
\describe{
  \item{ edge }{ A list whose first two unnamed members are the indices of
  the origin and destination vertices. Additional members must be named and
  are additional edge properties (e.g. length). }
  \item{ vertex }{ A list that optionally contain vertex properties, if any
  (or an empty list if none). }
}
}
\usage{
\method{print}{graph}(x, ...)
}
\arguments{
\item{x}{An object of \code{\link{graph-class}}.}

\item{...}{Additional parameters to be passed to the method. Currently
ignored.}
}
\description{
Class and methods to handle graphs.
}
\details{
Prints user-relevant information about the graph: number of edges
and vertices, edge and vertex labels, addition edge properties and vertex
properties.
}
\section{Functions}{
\itemize{
\item \code{print.graph}: Print method for graph-class objects
}}

\references{
Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
maps (PEM): a framework to model and predict species traits. Meth. Ecol.
Evol. 4: 1120--1131
}
\seealso{
\code{\link{PEM.build}}, \code{\link{PEM-class}}
}
\author{
\packageAuthor{MPSEM}
Maintainer: \packageMaintainer{MPSEM}
}
