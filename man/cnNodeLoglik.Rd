\name{cnNodeLoglik}
\alias{cnNodeLoglik}
\alias{cnNodeLoglik,catNetwork-method}

\title{Node Log-likelihood}

\description{For a given data sample, the function calculates the log-likelihood of a node 
with respect to a specified parent set.
}

\usage{
cnNodeLoglik(object, node, data, pert=NULL, softmode=FALSE, klmode=FALSE)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{node}{an \code{integer} or a \code{list} of \code{integer}s, node indices in the \code{data}}
\item{data}{a \code{matrix} or \code{data.frame} of categories}
\item{pert}{an optional perturbation \code{matrix} or \code{data.frame}}
\item{softmode}{a \code{logical}, turns on/off the soft quantization mode}
\item{klmode}{a \code{logical}, use the negative KL-distance or the log-likelihood}
}

\value{a \code{numeric} value}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxpars=3, numcats=2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  cnNodeLoglik(cnet, node=5, data=psamples)
}

\seealso{\code{\link{cnLoglik}} }

\keyword{methods}
\keyword{graphs}
