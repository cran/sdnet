\name{sdnet-package}
\alias{sdnet-package}
\alias{sdnet}

\docType{package}

\title{sdnet: Soft Discretization-based Bayesian Network Inference}

\description{
\code{sdnet} Fitting discrete Bayesian networks using marginal mixture distributions. 
A Bayesian network is defined by a graphical structure in form of directed acyclic graph and a probability model given as a set of conditional distributions, one for each node in the network. 
Considered in the package are only categorical Bayesian networks - networks which nodes represent discrete random variables.
The searching functions implemented in \code{sdnet} output sets of networks with increasing complexity that fit the data 
according to the MLE criterion. 
These optimal networks is believed to explain and represent the relations between the node-variables. The final network selection is left to the user.

Before starting to use the package, we suggest the user to take a look at some of the main objects used in \code{sdnet} such as \code{catNetwork} and \code{catNetworkEvaluate} and then familiarize with the main search functions \code{cnSearchOrder} and \code{cnSearchSA}. 
More details and examples can be found in the manual pages and the vignettes accompanying the package. 

Since \code{sdnet} does not have its own plotting abilities, the user needs to setup some external tools in order to visualize networks, or more precisely, \code{catNetwork} objects. 
There are two options: first one is to use \code{igraph} package and second, and better one, 
is to use \code{Graphviz} library. \code{Graphviz} is not a R-package but a platform independent library that the user have to install in advance on its machine in order to use this option.

If the user choose the first option, \code{igraph}, the only thing he/she has to do is to install the library in R and set the environment variable \code{R_SDNET_USE_IGRAPH} to TRUE. 
A convenient place to do this is in the R \code{.First} function

   .First <- function() \{

	......................

  	Sys.setenv(R_SDNET_USE_IGRAPH=TRUE)

   \}

In order to use \code{Graphviz}, in addition to installing the library, the user has to register a environmental variable with name \code{R_DOTVIEWER} with the path to the \code{Dot} executable file of \code{Graphviz}. The \code{Dot} routine generates a postscript or \code{pdf}-file from a text dot-file. 
Also, the user needs a postscript and \code{pdf}-viewer. The full path to it 
has to be given in another variable with name \code{R_PDFVIEWER}. 
Note that \code{R_PDFVIEWER} variable might be already setup. To check this call \code{Sys.getenv("R_PDFVIEWER")} in R.

The variables \code{R_DOTVIEWER} and eventually \code{R_PDFVIEWER} can be registered in the \code{.First} function residing in the \code{.Rprofile} initializing file.

Below we give two examples. 
On UNIX platform the user may use code like this one

   .First <- function() \{

	......................


  	Sys.setenv(R_DOTVIEWER="/usr/bin/dot")

   \}

On Windows platform the user may have the following two lines in its \code{.First} function 

  .First <- function() \{

	......................

	Sys.setenv(R_PDFVIEWER="\"C:/Program Files (x86)/Adobe/Reader 9.0/Reader/AcroRd32\"")

  	Sys.setenv(R_DOTVIEWER="\"C:/Program Files (x86)/Graphviz 2.26.3/bin/Dot\"")

 \}

  Note that all paths in Windows should be embraced by comment marks, "\"".

}

\author{N. Balov}

