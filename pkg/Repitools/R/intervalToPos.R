intervalToPos <- function(regions)
{
	if(!all(c("start", "end")  %in% colnames(regions)))
		stop("Incorrect column headings for regions. Check documentation for details.")
	regions$position <- round((regions$start + regions$end)/2)
	regions
}
