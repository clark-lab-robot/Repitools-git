\name{CoverageList}
\docType{class}
\alias{CoverageList}
\alias{CoverageList-class}
\alias{names,CoverageList-method}
\alias{names<-,CoverageList-method}
\alias{tables}
\alias{tables,CoverageList-method}
\alias{length,CoverageList-method}
\alias{show,CoverageList-method}
\alias{[,CoverageList-method}

\title{Container for \code{featureCoverage()} output.}
\description{
  Contains a list of coverage tables, and the parameters that were used to
  generate them.
}
\section{Accessors}{
  In the following code snippets, \code{x} is a CoverageList object.
  
  \describe{
    \item{}{
    	\code{names(x)}, \code{names(x) <- value}
	Gets and sets the experiment type names.
    }
    \item{}{
	\code{tables(x)}
	Gets the list of coverage matrices.
    }
    \item{}{
	\code{length(x)}
	Gets the number of coverage matrices.
    }
  }
}
\section{Subsetting}{
  In the following code snippets, \code{x} is a CoverageList object.

  \describe{
    \item{}{
	\code{x[i]}
	Creates a CoverageList object, keeping only the \code{i} matrices.
    }
  }
}

\author{Dario Strbenac}
