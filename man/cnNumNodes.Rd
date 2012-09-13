\name{cnNumNodes-method}
\alias{cnNumNodes}
\alias{cnNumNodes,catNetwork-method}

\title{Network Size}

\description{Returns the number of nodes of a \code{catNetwork} object.}

\usage{
cnNumNodes(object)
}

\arguments{
\item{object}{a \code{catNetwork}}
}

\value{an \code{integer}}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxpars=3, numcats=2)
  cnNumNodes(object=cnet)
}

\seealso{\code{\link{cnNodes}} }

\keyword{methods}
\keyword{graphs}
