\name{expres}
\alias{sdnLearn}
\alias{sdnPredict}
\alias{sdnEvaluate}

\title{Naive Classifier of Gene Expression Profiles}

\description{Implements a naive classifier using soft discretization. 
}

\usage{
sdnLearn(data, cls, clslevs = NULL, ncats = 3, nodeCats = NULL, quant="uniform", std=TRUE)
sdnPredict(model, data, std=TRUE) 
sdnEvaluate(train, test, ncats = 3, nodeCats = NULL, std=FALSE)
}

\arguments{
\item{data}{a numerical \code{matrix} in row-genes format}
\item{train}{a numerical \code{matrix} in row-genes format}
\item{test}{a numerical \code{matrix} in row-genes format}
\item{cls}{a \code{factor} or \code{integer}, the sample labels}
\item{clslevs}{an optional vector of labels, should include training data's labels}
\item{ncats}{an \code{integer}, the number of categories per node}
\item{nodeCats}{a \code{list}, custom node categories}
\item{quant}{quantization method}
\item{std}{a \code{logical}, should the data rows be standardized}
\item{model}{a \code{list} of components for the training model}
}

\details{
The \code{model} contains a vector of gene names \code{geneset}, a vector of sample labels \code{clslevs}, class \code{catNetworks}: \code{nets}, a list of node categories \code{nodeCats} and a training quantization model \code{quant}. 
}

\value{\code{sdnPredict} returns the log-ratio of the class conditional probabilities for each test observation. \code{sdnEvaluate} handles 2-class problems and returns the prediction accuracy and predicted classes.}

\author{N. Balov}

