\name{cnNodeSampleLoglik}
\alias{cnNodeSampleLoglik}
\alias{cnNodeSampleProb}

\title{Node Log-likelihood}

\description{For a given data sample, the function calculates the log-likelihood of a node 
with respect to a specified parent set.
}

\usage{
cnNodeSampleLoglik(node, pars, data, pert=NULL)
cnNodeSampleProb(node, pars, data, pert=NULL)
}

\arguments{
\item{node}{an \code{integer} or a \code{list} of \code{integer}s, node indices in the \code{data}}
\item{pars}{an \code{integer} or a \code{list} of \code{integer}s, vector of parent indices for the \code{node}s}
\item{data}{a \code{matrix} or \code{data.frame} of categories}
\item{pert}{an optional perturbation \code{matrix} or \code{data.frame}}
}

\value{a \code{numeric} value}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxpars=3, numcats=2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  cnNodeSampleLoglik(node=5, pars=c(1,2), data=psamples)
}

\seealso{\code{\link{cnLoglik}} }

\keyword{methods}
\keyword{graphs}
