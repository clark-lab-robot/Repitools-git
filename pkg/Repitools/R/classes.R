setClass("CoverageList", representation(
                                        marks = "character",
					cvgs = "list",
					anno = "GRanges",
					up = "numeric",
					down = "numeric",
					dist = "character",
					freq = "numeric",
					s.width = "numeric"
					))

setMethod("names", "CoverageList", function(x) x@marks)
setGeneric("tables", function(x) {standardGeneric("tables")})
setMethod("tables", "CoverageList", function(x) x@cvgs)
setMethod("length", "CoverageList", function(x) length(x@marks))

setMethod("show", "CoverageList",
    function(object)
    {
	distLabel <- ifelse(object@dist == "percent", '%', "bases")
	cat("An object of class 'CoverageList'.\n")
	cat("Tables: ", paste(object@marks, collapse = ", "), ".\n", sep = '')
	cat("Features:\n")
	print(object@anno)
	cat("Region:",  paste(object@up, distLabel, "up to", object@down,
                         distLabel, "down.\n"))
	cat("Smoothing:", paste(object@s.width, collapse = ", "), "bases.\n")
	cat("Sampling : ", object@freq, ' ', distLabel, ".\n\n",  sep = '')
    })

setMethod("[", "CoverageList",
    function(x, i)
    {
	new("CoverageList", marks = x@marks[i], anno = x@anno, cvgs = x@cvgs[i],
	                    up = x@up, down = x@down, dist = x@dist,
			    freq = x@freq, s.width = x@s.width)
	})

setReplaceMethod("names", "CoverageList",
    function(x, value)
    {
	if(length(value) != length(x@marks))
	    stop("New mark name(s) are a different length to previous mark
                  name(s).\n")
	x@marks <- value
	x
    }
)

setClass("ClusteredCoverageList", representation(
                                    cluster.id = "numeric",
				    expr = "numeric",
				    sort.data = "ANY",
				    sort.name = "ANY"),
                                 prototype(
				    sort.data = NULL,
				    sort.name = NULL),
                       contains = "CoverageList")

setMethod("show", "ClusteredCoverageList",
    function(object)
    {
	distLabel <- ifelse(object@dist == "percent", '%', "bases")
	cat("An object of class 'ClusteredCoverageList'.\n")
	cat("Tables: ", paste(object@marks, collapse = ", "), ".\n", sep = '')
	cat("Region: ",  paste(object@up, distLabel, "up to", object@down,
	    distLabel, "down.\n"))
	cat("Features:\n")
	head(object@anno)
	cat("Smoothing:", object@s.width, "bases.\n")
	cat("Sampling: ", object@freq, ' ', distLabel, ".\n",  sep = '')
	cat("Feature Expressions:", paste(paste(head(object@expr),
	    collapse = ", "), ", ...\n", sep = ''))
	cat("Feature Clusters:", paste(paste(head(object@cluster.id),
	    collapse = ", "), ", ...\n", sep = ''))
	if(!is.null(object@sort.data))
	    cat("Within Cluster Sorting: By ", object@sort.name, ". ",
	    paste(paste(head(object@sort.data), collapse = ", "), ", ...\n", sep = ''),
	    sep = '')		
    })

# Constructor
setGeneric("ClusteredCoverageList", function(x, ...)
           {standardGeneric("ClusteredCoverageList")})
setMethod("ClusteredCoverageList", "CoverageList",
    function(x, cvgs = tables(x), expr, cluster.id, sort.data = NULL,
             sort.name = NULL)
{
	new("ClusteredCoverageList", marks = x@marks, cvgs = cvgs, anno = x@anno,
	    up = x@up, down = x@down, dist = x@dist,
	    freq = x@freq, s.width = x@s.width, cluster.id = cluster.id,
	    expr = expr, sort.name = sort.name, sort.data = sort.data)
})

setMethod("[", "ClusteredCoverageList",
    function(x, i)
    {
	new("ClusteredCoverageList", marks = x@marks[i], cvgs = x@cvgs[i],
	    anno = x@anno, up = x@up, down = x@down, dist = x@dist,
	    freq = x@freq, s.width = x@s.width, cluster.id = x@cluster.id,
	    expr = x@expr, sort.data = x@sort.data, sort.name = x@sort.name)
    })

