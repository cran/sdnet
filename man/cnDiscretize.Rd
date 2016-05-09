\name{cnDiscretize}
\alias{cnDiscretize}

\title{Data Categorization}

\description{Soft and hard  discretization of numerical data.
}

\usage{
cnDiscretize(data, numcats=3, mode = "soft", marginal = "quantile",
learnset=NULL, cover=0.95, maxiter=100, eps=1e-8, weights=NULL)
}

\arguments{
\item{data}{a numerical \code{matrix} or \code{data.frame}}
\item{numcats}{an \code{integer}, the number of categories per node}
\item{mode}{a \code{character}, the discretization method to be used, "soft" or "hard"}
\item{marginal}{a \code{character}, the marginal model, "quantile", "uniform" or "gauss"}
\item{learnset}{observation indices to be used for setting the discretization parameters}
\item{cover}{a numerical between 0 and 1, proportion of the sample to be used for setting the discretization parameters}
\item{maxiter}{an \code{integer}, maximum iterations for the mixture of
  Gaussians EM algorithm}
\item{eps}{a \code{numerical}, convergence precision number}
\item{weights}{a \code{numerical}, sample record weights}
}

\details{
The numerical \code{data} is discretized into given number of categories, \code{numcats}, using the empirical node quantiles. 
As in all functions of \code{catnet} package that accept data, if the \code{data} parameter is a \code{matrix} then it is organized in the row-node format. If it is a \code{data.frame}, the column-node format is assumed.

The \code{mode} specifies the discretization model.
Currently, two hard discretization methods are supported - "quantile" and "uniform", which is the default choice.
The quantile-based discretization method is applied as follows. 
For each node, the sample node distribution is constructed, 
which is then represented by a sum of non-intersecting classes separated by the quantile points of the sample distribution. Each node value is assigned the class index in which it falls into. 

The uniform discretization breaks the range of values of each node into \code{numcats} equal intervals or of lengths proportional to the corresponding \code{qlevels} values. 

Currently, the function assigns equal number of categories for each node of the data.
}

\value{A \code{matrix} or \code{data.frame} of indices.}

\examples{
  ps <- t(sapply(1:10, function(i) rnorm(20, i, 0.1)))
  dps1 <- cnDiscretize(data=ps, numcats=3, mode="hard", marginal="quantile")
  hist(dps1[1,])
}

\seealso{\code{\link{cnSamples}}}

\author{N. Balov}

