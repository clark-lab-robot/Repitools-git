#FastQC-class
setClass("FastQC",
    representation(
        Basic_Statistics="data.frame", 
        Per_base_sequence_quality="data.frame",
        Per_sequence_quality_scores="data.frame",
        Per_base_sequence_content="data.frame",
        Per_base_GC_content="data.frame",
        Per_sequence_GC_content="data.frame",
        Per_base_N_content="data.frame",
        Sequence_Length_Distribution="data.frame",
        Sequence_Duplication_Levels="data.frame",
        Overrepresented_sequences="data.frame"
    )
)

#Generics for accessors
setGeneric("Basic_Statistics", function(x, ...) {standardGeneric("Basic_Statistics")})
setGeneric("Per_base_sequence_quality", function(x, ...) {standardGeneric("Per_base_sequence_quality")})
setGeneric("Per_sequence_quality_scores", function(x, ...) {standardGeneric("Per_sequence_quality_scores")})
setGeneric("Per_base_sequence_content", function(x, ...) {standardGeneric("Per_base_sequence_content")})
setGeneric("Per_base_GC_content", function(x, ...) {standardGeneric("Per_base_GC_content")})
setGeneric("Per_sequence_GC_content", function(x, ...) {standardGeneric("Per_sequence_GC_content")})
setGeneric("Per_base_N_content", function(x, ...) {standardGeneric("Per_base_N_content")})
setGeneric("Sequence_Length_Distribution", function(x, ...) {standardGeneric("Sequence_Length_Distribution")})
setGeneric("Sequence_Duplication_Levels", function(x, ...) {standardGeneric("Sequence_Duplication_Levels")})
setGeneric("Overrepresented_sequences", function(x, ...) {standardGeneric("Overrepresented_sequences")})
setGeneric("Mismatches", function(x, ...) {standardGeneric("Mismatches")})
setGeneric("MismatchTable", function(x, ...) {standardGeneric("MismatchTable")})

setMethod("show", "FastQC", 
    function(object) {
        cat("FastQC object of ", object@Basic_Statistics$Value[3], " ", object@Basic_Statistics$Value[4], "bp sequences, containing:\n", sep = "")
        cat(slotNames(object), sep="\n")
    }
)

#Accessors
setMethod("Basic_Statistics", "FastQC", function(x) x@Basic_Statistics)
setMethod("Per_base_sequence_quality", "FastQC", function(x) x@Per_base_sequence_quality)
setMethod("Per_sequence_quality_scores", "FastQC", function(x) x@Per_sequence_quality_scores)
setMethod("Per_base_sequence_content", "FastQC", function(x) x@Per_base_sequence_content)
setMethod("Per_base_GC_content", "FastQC", function(x) x@Per_base_GC_content)
setMethod("Per_sequence_GC_content", "FastQC", function(x) x@Per_sequence_GC_content)
setMethod("Per_base_N_content", "FastQC", function(x) x@Per_base_N_content)
setMethod("Sequence_Length_Distribution", "FastQC", function(x) x@Sequence_Length_Distribution)
setMethod("Sequence_Duplication_Levels", "FastQC", function(x) x@Sequence_Duplication_Levels)
setMethod("Overrepresented_sequences", "FastQC", function(x) x@Overrepresented_sequences)

readFastQC <- function(filename) {
    if (length(filename)!=1) stop("Length of filename must be 1")
    if (!file.exists(filename)) stop(filename, "does not exist")
    temp <- readLines(filename)
    #remove #'s
    temp <- gsub("#", "", temp)
    #remove END_MODULE lines
    temp <- temp[!grepl(">>END_MODULE", temp)]
    #split into modules
    temp <- split(temp, cumsum(grepl("^>>", temp)))[-1]
    names(temp) <- sapply(temp, function(x) {
                gsub("^>>", "", gsub("\t.*", "", gsub(" ", "_", x[1])))
            })
    #Parse each module
    temp <- lapply(temp, function(x) {
                if (length(x)==1) return(data.frame())
                x <- strsplit(x[-1], split="\t")
                tab <- as.data.frame(do.call(rbind, x[-1]), stringsAsFactors=FALSE)
                for (i in 1:ncol(tab)) if (!any(is.na(suppressWarnings(as.numeric(tab[,i]))))) tab[,i] <- as.numeric(tab[,i]) 
                colnames(tab) <- x[[1]]
                tab
            })
    do.call(new, c(list(Class="FastQC"), temp))
}

##SequenceQC-class
setClass("SequenceQC",
    representation(
        Unaligned="FastQC",
        Aligned="FastQC",
        Mismatches="data.frame",
        MismatchTable="data.frame"
    )
)

setMethod("show", "SequenceQC", function(object) {
    cat("SequenceQC object of ", object@Unaligned@Basic_Statistics$Value[4], "bp sequences\n", sep="")
    cat(object@Unaligned@Basic_Statistics$Value[3], " before alignment\n", sep="")
    cat(object@Aligned@Basic_Statistics$Value[3], " after alignment\n", sep="")
})

#Accessors
setMethod("Basic_Statistics", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Basic_Statistics(x@Aligned) else Basic_Statistics(x@Unaligned)
})
setMethod("Per_base_sequence_quality", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Per_base_sequence_quality(x@Aligned) else Per_base_sequence_quality(x@Unaligned)
})
setMethod("Per_sequence_quality_scores", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Per_sequence_quality_scores(x@Aligned) else Per_sequence_quality_scores(x@Unaligned)
})
setMethod("Per_base_sequence_content", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Per_base_sequence_content(x@Aligned) else Per_base_sequence_content(x@Unaligned)
})
setMethod("Per_base_GC_content", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Per_base_GC_content(x@Aligned) else Per_base_GC_content(x@Unaligned)
})
setMethod("Per_sequence_GC_content", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Per_sequence_GC_content(x@Aligned) else Per_sequence_GC_content(x@Unaligned)
})
setMethod("Per_base_N_content", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Per_base_N_content(x@Aligned) else Per_base_N_content(x@Unaligned)
})
setMethod("Sequence_Duplication_Levels", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Sequence_Duplication_Levels(x@Aligned) else Sequence_Duplication_Levels(x@Unaligned)
})
setMethod("Overrepresented_sequences", "SequenceQC", function(x, which=c("Unaligned", "Aligned")) {
    which <- match.arg(which)
    if (which=="Aligned") Overrepresented_sequences(x@Aligned) else Overrepresented_sequences(x@Unaligned)
})
setMethod("Mismatches", "SequenceQC", function(x) x@Mismatches)
setMethod("MismatchTable", "SequenceQC", function(x) x@MismatchTable)

##SequenceQCSet-class
setClass("SequenceQCSet",
    contains="list",
    prototype=list()
)

setValidity("SequenceQCSet", function(object) {
    if (!all(sapply(object@.Data, class)=="SequenceQC"))
        return("All elements must be 'SequenceQC' objects")
    TRUE
})

setMethod("show", "SequenceQCSet", function(object) {
    cat("SequenceQCSet object containing ", length(object@.Data), " element(s)\n", sep="")
    cat(paste(names(object), collapse=" "), "\n")
})

#Accessors
setMethod("Basic_Statistics", "SequenceQCSet", function(x, ...) lapply(x, Basic_Statistics, ...))
setMethod("Per_base_sequence_quality", "SequenceQCSet", function(x, ...) lapply(x, Per_base_sequence_quality, ...))
setMethod("Per_sequence_quality_scores", "SequenceQCSet", function(x, ...) lapply(x, Per_sequence_quality_scores, ...))
setMethod("Per_base_sequence_content", "SequenceQCSet", function(x, ...) lapply(x, Per_base_sequence_content, ...))
setMethod("Per_base_GC_content", "SequenceQCSet", function(x, ...) lapply(x, Per_base_GC_content, ...))
setMethod("Per_sequence_GC_content", "SequenceQCSet", function(x, ...) lapply(x, Per_sequence_GC_content, ...))
setMethod("Per_base_N_content", "SequenceQCSet", function(x, ...) lapply(x, Per_base_N_content, ...))
setMethod("Sequence_Length_Distribution", "SequenceQCSet", function(x, ...) lapply(x, Sequence_Length_Distribution, ...))
setMethod("Sequence_Duplication_Levels", "SequenceQCSet", function(x, ...) lapply(x, Sequence_Duplication_Levels, ...))
setMethod("Overrepresented_sequences", "SequenceQCSet", function(x, ...) lapply(x, Overrepresented_sequences, ...))
setMethod("Mismatches", "SequenceQCSet", function(x) lapply(x, Mismatches))
setMethod("MismatchTable", "SequenceQCSet", function(x) lapply(x, MismatchTable))

setMethod("[", "SequenceQCSet", function(x, i, j, ..., drop){
    newSet <- x@.Data[i]
    names(newSet) <- names(x)[i]
    new("SequenceQCSet", newSet)
})
