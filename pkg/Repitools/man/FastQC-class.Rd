\name{FastQC-class}
\docType{class}
\alias{class:FastQC}
\alias{FastQC-class}
\alias{class:SequenceQC}
\alias{SequenceQC-class}
\alias{class:SequenceQCSet}
\alias{SequenceQCSet-class}
\alias{readFastQC}
\alias{SequenceQCSet}
\alias{[,SequenceQCSet-method}
\alias{[[,SequenceQCSet-method}
\alias{names,SequenceQCSet-method}
\title{FastQC and associated classes}

\description{
	The FastQC class stores results obtained from the FastQC application (http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/),
	with a slot for each FastQC module. The SequenceQC class contains the QC results of a single lane of sequencing in three slots:
	\code{Unaligned} - FastQC results obtained from all reads (before alignment) 
	\code{Aligned} - FastQC results obtained from only reads which aligned
	\code{Mismatches} - a \code{data.frame} containing counts for the number of mismatches of each type found at each sequencing cycle   
}

\author{Aaron Statham}