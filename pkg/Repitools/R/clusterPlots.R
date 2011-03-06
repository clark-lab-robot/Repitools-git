setGeneric("clusterPlots", function(c.list, ...){standardGeneric("clusterPlots")})

setMethod("clusterPlots", "ClusteredCoverageList",
	function(c.list, plot.ord = 1:length(c.list), plot.type = c("line", "heatmap"),
                 cols = NULL, t.name = "Clustered Enrichment")
{
    plot.type <- match.arg(plot.type)
    nMarks <- length(c.list)

    if(is.null(cols) == TRUE)
    {
	require(gplots)
	if(plot.type == "line")
	{
	    cols <- colorpanel(nMarks, "blue", "green", "red")
	} else {
	    cols <- colorpanel(100, "blue", "white", "red")
	}
    }

    cvgs <- tables(c.list)

    # Get median expression for each cluster. Find ascending order.
    expr <- c.list@expr
    clID <- c.list@cluster.id
    n.clusters <- length(unique(clID))
    clExpr <- tapply(expr, factor(clID), median, na.rm = TRUE)
    clOrd <- order(clExpr)

    # Get x-axis pos and labels.
    posLabels <- colnames(cvgs[[1]])
    pos <- as.integer(gsub('%', '', posLabels)) # Get raw position if labels have
                                                # percentage signs.

    if(plot.type == "line")
    {
	# Group each cluster from all epigenetic marks.

	profiles <- list() 
	for(i in 1:n.clusters)
	    profiles[[i]] <- sapply(cvgs, function(x)
                                    colMeans(x[clID == clOrd[i], , drop = FALSE]))

	# Plot the lineplots.

	invisible(mapply(function(x, y)
	{
		layout(rbind(1:2), widths=c(3, 1))
		par(mai=c(1, 1, 0.8, 0.1))
			
		yMax = max(unlist(x)) * 1.1
		matplot(pos, x, type = 'l', lty = 1, col = cols, ylim = c(0, yMax),
                        xlab = "Position", ylab = "Smoothed Coverage", yaxt = 'n',
                        xaxs = "i", yaxs = "i",
                        main = paste("Within Cluster Coverage (Median Expression : ",
                        round(y, 2), ")", sep = ''))
		axis(2, at = c(0, yMax / 2, yMax), labels = c("None", "Medium", "High"))

		plot.new()
		par(mai=c(1,0.1,0.2,0.1))
		legend("topleft", legend = names(c.list), title = "Mark", col = cols, lwd = 2)
	}, profiles, clExpr[clOrd]))

    } else { # Plot a heatmap
	
	sort.data <- c.list@sort.data
	# Get order of all features next.
	if(length(sort.data) == 0)
		ord <- order(factor(clID, levels = clOrd))
	else
		ord <- order(factor(clID, levels = clOrd), sort.data)
	
	# Re-arrange the ChIP and expression data and vector of sort data.
	cvgs <- lapply(cvgs, function(x) x[ord, ])
	expr <- expr[ord]
	if(length(sort.data) > 0) sort.data <- sort.data[ord]

	# Plot heatmap.
	if(length(sort.data) > 0)
		layout(rbind(1:(nMarks + 3)), widths=c(1, rep(3, nMarks), 2, 1))
	else
		layout(rbind(1:(nMarks + 2)), widths=c(1, rep(3, nMarks), 2))
	par(mai=c(1.02,0.50,0.82,0.05))
  	
	nBins = length(cols)
	image(y=seq(1/nBins/2, 1-(1/nBins/2), 1/nBins), z=rbind(1:nBins),
              col = cols, axes = FALSE, xlab = "Read Coverage", ylab = NA)
	axis(2, at=c(0, 0.5, 1), labels=c("None", "Medium", "High"))

	par(mai=c(1.02,0.05,0.82,0.05))
	mapply(function(x, y)
	{
	    image(pos, 1:nrow(x), t(x), xlab="Position relative to feature",
                  xaxt = "n", yaxt = "n", ylab = "Feature", col = cols, main = y)
	    axis(1, pos, labels = posLabels)		

	    # Add lines delimiting the cluster boundaries.
	    bounds <- cumsum(table(clID)[clOrd])[-n.clusters]
	    abline(h = bounds, lwd = 2)
	}, cvgs, names(c.list))

	par(mai=c(1.02,0.05,0.82,0.50))
	plot(expr, y = 1:length(expr), yaxs = 'i', pch = 19, xlab = "log2 expression",
             ylab = NA, yaxt = "n", cex = 0.5)
	if(!is.null(sort.data)) plot(sort.data, y = 1:length(sort.data), yaxs = 'i',
                     pch = 19, xlab = sort.name, ylab = NA, yaxt = "n", cex = 0.5)	

	mtext(t.name, line = 0, outer = TRUE)
    }
})

setMethod("clusterPlots", "CoverageList", function(c.list, scale = function(x) x,
          capQ = 0.95, cap.type = c("sep", "all"), n.clusters = 5,
          plot.ord = 1:length(c.list), expr, sort.data = NULL, sort.name = NULL,
          plot.type = c("line", "heatmap"), cols = NULL, t.name = "Clustered Enrichment")
{
	plot.type <- match.arg(plot.type)
	cap.type <- match.arg(cap.type)
	cvgs <- tables(c.list)

	if(nrow(cvgs[[1]]) != length(expr))
		stop("Number of features does not match length of expression vector.\n")
	if(!is.null(sort.data) && is.null(sort.name))
		stop("'sort.data' provided, 'sort.name' is not.")
	if(!is.null(sort.data) && length(expr) != length(sort.data))
		stop("'sort.data' length not the same as number of features in coverage
                      matrices or expression data.\n")

	# Precision sometimes means 0 is represented as very small negative numbers.
	cvgs <- lapply(cvgs, function(x) {x[x < 0] = 0; x})

	cvgs <- lapply(cvgs, scale)
	if(cap.type == "all") maxCvg <- quantile(do.call(cbind, cvgs), capQ)

	# Find the maximum score allowable, then make any scores bigger be the maximum score.
	cvgs <- lapply(cvgs, function(x)
	{
		if(cap.type == "sep") maxCvg <- quantile(x, capQ)
		x[x > maxCvg] = maxCvg
		return(x)
	})

	# Do the k-means clustering for all marks together.
	set.seed(100)
	all <- do.call(cbind, cvgs)
	clID <- kmeans(all, n.clusters, iter.max = 100)$cluster

	ccl <- ClusteredCoverageList(c.list, cvgs = cvgs, cluster.id = clID, expr = expr,
                                     sort.data = sort.data, sort.name = sort.name)
	clusterPlots(ccl, plot.ord = plot.ord, plot.type = plot.type, cols = cols, t.name = t.name)
	invisible(ccl)
})
