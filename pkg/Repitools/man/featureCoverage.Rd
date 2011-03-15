\name{featureCoverage}
\alias{featureCoverage}
\alias{featureCoverage,ANY,data.frame-method}
\alias{featureCoverage,GenomeDataList,ANY-method}
\alias{featureCoverage,character,GRanges-method}
\alias{featureCoverage,GRangesList,GRanges-method}
\alias{featureCoverage,GRanges,GRanges-method}

\title{Get coverages at regular sample points around genomic features.}
\description{
  Given a \code{GRanges} / \code{GenomeDataList} / \code{GRangesList} object,
  or BAM file paths, of reads for each experimental condition, and an anntotation
  of features of interest, the smoothed coverage of reads at regularly spaced
  positions around the features is calculated.
}
\usage{
  \S4method{featureCoverage}{ANY,data.frame}(
    x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
    s.width = 1000, verbose = TRUE)
  \S4method{featureCoverage}{GenomeDataList,ANY}(
    x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
    s.width = 1000, verbose = TRUE)
  \S4method{featureCoverage}{character,GRanges}(
    x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
    s.width = 1000, verbose = TRUE)  
  \S4method{featureCoverage}{GRangesList,GRanges}(
    x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
    s.width = 1000, verbose = TRUE)
  \S4method{featureCoverage}{GRanges,GRanges}(
    x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
    s.width = 1000, verbose = TRUE)
}
\arguments{
  \item{x}{If it is a \code{character} vector, then paths to BAM files, otherwise
           a collection of mapped short reads.}
  \item{anno}{Annotation of the features to sample around.}
  \item{up}{How far to go from the features' reference points in one direction.}
  \item{down}{How far to go from the features' reference points in the opposite
              direction.}
  \item{dist}{The type of distance measure to use, in determining the boundaries
              of the sampling area.}
  \item{freq}{Coverage sampling frequency.}
  \item{s.width}{The width of smoothing to apply to the coverage. Either length 1
                , or a numeric vector of the same length as \code{x}.}
  \item{verbose}{Whether to print the progess of processing.}
}
\value{
  A \code{\linkS4class{CoverageList}} object, that holds a list of coverage
  matrices, one for each experiment type, and the parameters that were used
  to create the coverage matrices
}
\details{
  If \code{x} is a vector of paths or a \code{GenomeDataList} or \code{GRangesList}
  object, then \code{names(x)} should contain the types of the experiments.

  If \code{anno} is a \code{data.frame}, it must contan the columns \code{chr},
  \code{start}, and \code{end}. Optional columns are \code{strand} and \code{name}.
  If \code{anno} is a \code{GRanges} object, then the name can be present as a column
  called \code{name} in the element metadata of the GRanges object. If names
  are given, then the coverage matrices will use the names as their row names.

  An approximation to running mean smoothing of the coverage is used. Reads are
  extended to the smoothing width, rather than to their fragment size, and
  coverage is used directly. This method is faster than a running mean of the
  calculated coverage, and qualtatively almost identical.

  The annotation can be stranded or not. if the annotation is stranded, then
  the reference point is the start coordinate for features on the + strand,
  and the end coordinate for features on the - strand. If the annotation is
  unstranded (e.g. annotation of CpG islands), then the midpoint of the feature
  is used for the reference point.

  The \code{up} and \code{down} values give how far up and down from the
  reference point to find sample smoothed coverages. The semantics of them depend
  on if the annotation is stranded or not. If the annotation is stranded, then
  they give how far upstream and downstream will be sampled. If the annotation is
  unstranded, then \code{"up"} gives how far towards the start of a chromosome to go,
  and \code{"down"} gives how far towards the end of a chromosome to go. If
  \code{dist} is \code{"percent"}, then they give how many percent of
  each feature's width away from the reference point the sampling boundaries are.
  If \code{dist} is \code{"base"}, then the boundaries of the sampling region
  are a fixed width for every feature, and the units of \code{"up"} and \code{"down"}
  are bases. \code{"up"} and \code{"down"} must be identical if the features
  are unstranded. The units of \code{freq} are percent for \code{dist}
  being \code{"percent"}, and bases for \code{dist} being {"base"}.

  Sampling points that are off the edge of a chromosome (e.g. for a gene that
  has a TSS very near to the edge of a chromosome) are given a coverage of 0.

  There is no correction applied for the total number of reads in each
  \code{GRanges} / \code{GenomeData} / BAM file, or for the fact that different marks
  have different peak profiles (i.e. sharp, high peaks vs. low, spread out peaks),
  so different coverage matrices are not directly comparable.
}

\seealso{\code{\link{mergeReplicates}} for merging replicates of an experiment type.}
\author{Dario Strbenac}
