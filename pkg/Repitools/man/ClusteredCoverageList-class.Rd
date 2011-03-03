\name{ClusteredCoverageList}
\docType{class}
\alias{ClusteredCoverageList}
\alias{ClusteredCoverageList-class}
\alias{ClusteredCoverageList,CoverageList-method}
\alias{show,ClusteredCoverageList-method}
\alias{[,ClusteredCoverageList-method}

\title{Container for coverage matrices with clustering results.}
\description{
  Contains a list of coverage matrices, the parameters that were used to
  generate them origin, and also cluster membership and expression data.

  This allows \code{\link{clusterPlots}} to be called multiple times without
  having to redo the clustering each time. It also allows the user to take
  the \code{\linkS4class{CoverageList}} output of \code{\link{featureCoverage}},
  and do their own custom clustering on the coverage matrices, then save the
  clustering results in this container.
}
\section{Constructor}{
  \describe{
    \item{}{
    	\code{ClusteredCoverageList(x, cvgs = tables(x), expr, cluster.id,
                                    sort.name = NULL, sort.data = NULL)}
	Creates a ClusteredCoverageList object.
	\describe{
	    \item{\code{x}}{A \code{\linkS4class{CoverageList}} object.}
	    \item{\code{cvgs}}{A list of coverage matrices.}
	    \item{\code{expr}}{A numeric vector, same length as number of
                               rows of every coverage matrix.}
	    \item{\code{cluster.id}}{A vector, same length as number of
                                     rows of every coverage matrix.}
	    \item{\code{sort.data}}{Vector of data to order features
                                     within clusters by.}
	    \item{\code{sort.name}}{Human readable description of what
                                     the sorting data is of.}
	}
    }
  }
}
\section{Subsetting}{
  In the following code snippets, \code{x} is a ClusteredCoverageList object.

  \describe{
    \item{}{
	\code{x[i]}
	Creates a ClusteredCoverageList object, keeping only the \code{i} matrices.
    }
  }
}

\author{Dario Strbenac}
