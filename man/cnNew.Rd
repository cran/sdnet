\name{cnNew}
\alias{cnNew}

\title{New catNetwork}

\description{Creates a new \code{catNetwork} with specified nodes, categories, parent sets and probability structure.}

\usage{
cnNew(nodes, cats, pars, probs=NULL, p.delta1=0.01, p.delta2=0.01, dagonly=FALSE)
}

\arguments{
\item{nodes}{a \code{vector} of nodes names}
\item{cats}{a \code{list} of node categories}
\item{pars}{a \code{list} of node pars}
\item{probs}{a \code{list} of probabilities}
\item{p.delta1}{a \code{numeric}}
\item{p.delta2}{a \code{numeric}}
\item{dagonly}{a \code{logical}, selects between catNetwork and DAG}
}

\details{
If \code{probs} is not specified, then a random probability model is assigned with conditional probability values in the union of the intervals [p.delta1, 0.5-p.delta2] and [0.5+p.delta2, 1-p.delta1].
Because of the nested list hierarchy of the probability structure, 
specifying the probability argument explicitly can be very elaborated task for large networks. 
In the following example we create a small network with only three nodes. 
The first node has no pars and only its marginal distribution is given, \code{c(0.2,0.8)}.
Note that all inner most vectors in the \code{probs} argument, such as \code{(0.4,0.6)}, 
represent conditional distributions and thus sum to 1. 
}

\value{A \code{catNetwork} object.}

\author{N. Balov}

\examples{
cnet <- cnNew(
	nodes = c("a", "b", "c"),
	cats = list(c("1","2"), c("1","2"), c("1","2")), 
	pars = list(NULL, c(1), c(1,2)), 
	probs = list(	c(0.2,0.8), 
			list(c(0.6,0.4),c(0.4,0.6)), 
			list(list(c(0.3,0.7),c(0.7,0.3)), 
			list(c(0.9,0.1),c(0.1,0.9))))
	)
}

\seealso{\code{\link{catNetwork-class}}, \code{\link{cnRandomCatnet}} }

\keyword{classes}
\keyword{graphs}
