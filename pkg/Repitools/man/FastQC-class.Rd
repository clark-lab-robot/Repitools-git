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

\alias{Basic_Statistics}
\alias{Per_base_sequence_quality}
\alias{Per_sequence_quality_scores}
\alias{Per_base_sequence_content}
\alias{Per_base_GC_content}
\alias{Per_sequence_GC_content}
\alias{Per_base_N_content}
\alias{Sequence_Length_Distribution}
\alias{Sequence_Duplication_Levels}
\alias{Overrepresented_sequences}
\alias{Mismatches}
\alias{MismatchTable}

\alias{show,FastQC-method}
\alias{Basic_Statistics,FastQC-method}
\alias{Per_base_sequence_quality,FastQC-method}
\alias{Per_sequence_quality_scores,FastQC-method}
\alias{Per_base_sequence_content,FastQC-method}
\alias{Per_base_GC_content,FastQC-method}
\alias{Per_sequence_GC_content,FastQC-method}
\alias{Per_base_N_content,FastQC-method}
\alias{Sequence_Length_Distribution,FastQC-method}
\alias{Sequence_Duplication_Levels,FastQC-method}
\alias{Overrepresented_sequences,FastQC-method}

\alias{show,SequenceQC-method}
\alias{Basic_Statistics,SequenceQC-method}
\alias{Per_base_sequence_quality,SequenceQC-method}
\alias{Per_sequence_quality_scores,SequenceQC-method}
\alias{Per_base_sequence_content,SequenceQC-method}
\alias{Per_base_GC_content,SequenceQC-method}
\alias{Per_sequence_GC_content,SequenceQC-method}
\alias{Per_base_N_content,SequenceQC-method}
\alias{Sequence_Length_Distribution,SequenceQC-method}
\alias{Sequence_Duplication_Levels,SequenceQC-method}
\alias{Overrepresented_sequences,SequenceQC-method}
\alias{Mismatches,SequenceQC-method}
\alias{MismatchTable,SequenceQC-method}

\alias{show,SequenceQCSet-method}
\alias{Basic_Statistics,SequenceQCSet-method}
\alias{Per_base_sequence_quality,SequenceQCSet-method}
\alias{Per_sequence_quality_scores,SequenceQCSet-method}
\alias{Per_base_sequence_content,SequenceQCSet-method}
\alias{Per_base_GC_content,SequenceQCSet-method}
\alias{Per_sequence_GC_content,SequenceQCSet-method}
\alias{Per_base_N_content,SequenceQCSet-method}
\alias{Sequence_Length_Distribution,SequenceQCSet-method}
\alias{Sequence_Duplication_Levels,SequenceQCSet-method}
\alias{Overrepresented_sequences,SequenceQCSet-method}
\alias{Mismatches,SequenceQCSet-method}
\alias{MismatchTable,SequenceQCSet-method}
\alias{[,SequenceQCSet-method}
 

\title{FastQC and associated classes}

\description{
	The FastQC class stores results obtained from the FastQC application (see \code{references}),
	with a slot for each FastQC module. The SequenceQC class contains the QC results of a single lane of sequencing in three slots:
	\code{Unaligned} - FastQC results obtained from all reads (before alignment) 
	\code{Aligned} - FastQC results obtained from only reads which aligned
	\code{Mismatches} - a \code{data.frame} containing counts for the number of mismatches of each type found at each sequencing cycle   
}

\section{Slots of a FastQC object}{
	\code{Basic_Statistics}
	
	\code{Per_base_sequence_quality}
	
	\code{Per_sequence_quality_scores}
	
	\code{Per_base_sequence_content}
	
	\code{Per_base_GC_content}
	
	\code{Per_sequence_GC_content}
	
	\code{Per_base_N_content}
	
	\code{Sequence_Length_Distribution}
	
	\code{Sequence_Duplication_Levels}
	
	\code{Overrepresented_sequences}
}

\section{Slots of a SequenceQC object}{
	\code{Unaligned} - FastQC results obtained from all reads (before alignment)
	 
	\code{Aligned} - FastQC results obtained from only reads which aligned
	
	\code{Mismatches} - a \code{data.frame} containing counts for the number of mismatches of each type found at each sequencing cycle

	\code{MismatchTable} - a \code{data.frame} containing counts of how many mismatches aligned sequences contain
}

\author{Aaron Statham}

\references{FastQC - http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/}