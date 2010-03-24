getRMAdata <- function(name, chipType, units = NULL, doCore = FALSE, doNorm = FALSE, doExon = FALSE, doMerge = FALSE, doResiduals = FALSE, doWeights = FALSE, force = FALSE)
{
	if(doCore == TRUE)
		chipType <- paste(chipType, ",core", sep = "")
	cdf <- AffymetrixCdfFile$byChipType(chipType)
	celSet <- AffymetrixCelSet$byName(name, cdf = cdf)
	if(doNorm == TRUE)
	{
		bcProcess <- RmaBackgroundCorrection(celSet)
		bcCelSet <- process(bcProcess, verbose = verbose, force = force)
		qnProcess <- QuantileNormalization(bcCelSet, typesToUpdate="pm")
		qnCelSet <- process(qnProcess, verbose = verbose, force = force)
	}
	if(doExon == TRUE)
		plm <- ExonRmaPlm(qnCelSet, mergeGroups = doMerge) # Exons into transcripts.
	else
		plm <- RmaPlm(qnCelSet)
	fit(plm, verbose = verbose, units = units)
	qaModel <- QualityAssessmentModel(plm)
	chipEffectsSet <- getChipEffectSet(plm)
	chipEffectsDataFrame <- extractDataFrame(chipEffectsSet, addNames = TRUE, verbose = verbose, units = units)
	if(doResiduals == TRUE)
		res <- calculateResiduals(plm, verbose = verbose, force = force)
	if(doWeights == TRUE)
		wts <- calculateWeights(plm, verbose = verbose, force = force)

	list(name = name, chipType = chipType, exprs = chipEffectsDataFrame, cdf = cdf, plm = plm, qam = qaModel, residuals = res, weights = wts)
}
