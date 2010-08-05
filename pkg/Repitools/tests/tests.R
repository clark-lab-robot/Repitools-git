require("Repitools")
require("BSgenome.Hsapiens.UCSC.hg18")
require("chipseq")
options(warn = -1)

pathToData <- system.file("exampleData", package = "RepitoolsExamples")

if(pathToData == "")
{
	cat("Running minimal tests on functions that don't need datafiles.\n")
}

probes <- data.frame(chr = c("chr1", "chr9", "chrY", "chr1", "chr21", "chr6", "chr6", "chr2", "chrX", "chr11"), position = c(10000, 5500, 100000, 11000, 20000000, 500100, 499900, 700000, 9900, 90000), strand = c('+', '+', '+', '+', '-', '-', '-', '+', '-', '+'), stringsAsFactors = FALSE)
genes <- data.frame(chr = c("chr1", "chr9", "chr11", "chr1", "chr11", "chr6", "chr6", "chr22", "chrY", "chr21"), start = c(10000, 7900, 950000, 10500, 74000000, 450000, 5000000, 44000000, 1500, 9800000), end = c(12500, 9500, 1000000, 14500, 75000000, 500000, 9000000, 45000000, 3000, 10000000), strand = c('+', '-', '-', '+', '-', '-', '-', '+', '+', '-'), name = paste("Gene", 1:10), stringsAsFactors = FALSE)

crossMatch <- annotationLookup(probes, genes, 5000, 5000)
correctCrossMatch <- list(indexes = list(`Gene 1` = as.integer(c(1, 4)), `Gene 2` = as.integer(c(2)), `Gene 3` = NULL, `Gene 4` = as.integer(c(1, 4)), `Gene 5` = NULL, `Gene 6` = as.integer(c(6, 7)), `Gene 7` = NULL, `Gene 8` = NULL, `Gene 9` = NULL, `Gene 10` = NULL), offsets = list(`Gene 1` = as.integer(c(0, 1000)), `Gene 2` = as.integer(c(4000)), `Gene 3` = numeric(), `Gene 4` = as.integer(c(-500, 500)), `Gene 5` = numeric(), `Gene 6` = as.integer(c(-100, 100)), `Gene 7` = numeric(), `Gene 8` = numeric(), `Gene 9` = numeric(), `Gene 10` = numeric()))
names(correctCrossMatch$offsets$`Gene 1`) <- c(1, 4)
names(correctCrossMatch$offsets$`Gene 2`) <- c(2)
names(correctCrossMatch$offsets$`Gene 4`) <- c(1, 4)
names(correctCrossMatch$offsets$`Gene 6`) <- c(6, 7)

if(!isTRUE(all.equal(crossMatch, correctCrossMatch))) 
	stop("Error in annotationLookup function.")
cat("anontationLookup tested fine.\n")

lookupTable <- makeWindowLookupTable(crossMatch$indexes, crossMatch$offsets, starts = seq(-5000, 4900, 100), ends = seq(-4900, 5000, 100))
correctLookupTable <- matrix(NA, nrow = 10, ncol = 100, dimnames = list(genes$names, seq(-4950, 4950, 100)))
correctLookupTable[1, c(50, 51)] <- 1
correctLookupTable[1, c(60, 61)] <- 4
correctLookupTable[2, c(90, 91)] <- 2
correctLookupTable[4, c(45, 46)] <- 1
correctLookupTable[4, c(55, 56)] <- 4
correctLookupTable[6, c(49, 50)] <- 6
correctLookupTable[6, c(51, 52)] <- 7

if(!all(lookupTable == correctLookupTable, na.rm = TRUE))
	stop("Error in makeWindowLookupTable function")
cat("makeWindowLookupTable tested fine.\n")

cpgDensity <- cpgDensityCalc(genes, organism = Hsapiens)
if(!isTRUE(all.equal(cpgDensity, c(5.784, 7.620, 5.828, 2.928, 2.080, 1.252, 0.000, 7.404, 3.928, 0.000))))
	stop("cpgDensityCalc not working for window = 500, scaling = linear")
cpgDensity <- cpgDensityCalc(genes, window = 100, wFunction = "log", organism = Hsapiens)
if(!isTRUE(all.equal(round(cpgDensity, 3), c(2.424, 1.882, 1.436, 0.084, 0.379, 0.000, 0.000, 0.263, 1.392, 0.000))))
	stop("cpgDensityCalc not working for window = 100, scaling = log")
cpgDensity <- cpgDensityCalc(genes, window = 1000, wFunction = "exp", organism = Hsapiens)
if(!isTRUE(all.equal(round(cpgDensity, 3), c(4.874, 5.828, 4.999, 2.239, 1.567, 0.851, 0.054, 5.589, 3.229,0.062))))
	stop("cpgDensityCalc not working for window = 1000, scaling = exp")
cpgDensity <- cpgDensityCalc(genes, window = 500, wFunction = "none", organism = Hsapiens)
if(!isTRUE(all.equal(cpgDensity, c(11, 14, 16, 6, 4, 2, 0, 15, 9, 0))))
	stop("cpgDensityCalc not working for window = 500, scaling = none")
cat("cpgDensityCalc tested fine.\n")

GCpercent <- gcContentCalc(genes, organism = Hsapiens)
if(!isTRUE(all.equal(GCpercent, c(0.504, 0.586, 0.558, 0.470, 0.538, 0.304, 0.356, 0.636, 0.442, 0.388))))
	stop("Error in gcContentCalc function")
cat("gcContentCalc tested fine.\n")

findsCount <- sequenceCalc(genes, organism = Hsapiens, pattern = "AATT")
if(!isTRUE(all.equal(findsCount, c(1, 1, 0, 2, 1, 8, 2, 0, 4, 10))))
		stop("Error in sequenceCalc function counting task")

findsPlaces <- sequenceCalc(genes, organism = Hsapiens, pattern = "AATT", positions = TRUE)
correctPlaces <- list(-62, 181, NULL, c(-140, -98), 231, c(-219, -146, -88, -12, 12, 182, 209, 214), c(-61, 60), NULL, c(-115, -30, 11, 80), c(-238, -228, -202, -189, -177, -106, -21, 148, 158, 238))
if(!isTRUE(all.equal(findsPlaces, correctPlaces)))
		stop("Error in sequenceCalc function positions task")
cat("sequenceCalc tested fine.\n")

if(pathToData != "") # Do tests involving CDFs, CELs, sequences.
{
	require("aroma.affymetrix")
	require("RepitoolsExamples")
	require("edgeR")
	userWd <- getwd()
	setwd(pathToData)

	cdfFile <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02")
	cdfTableUniquePositions <- getUniqueCdf(cdfFile)
	celSetTables <- AffymetrixCelSet$byName("Tiling", cdf = cdfFile)
	MATNormalise <- MatNormalization(celSetTables, numChunks = 15)
	celSetMATNormalised <- process(MATNormalise)
	celSetMATNormalisedUniquePositions <- convertToUnique(celSetMATNormalised)
	experimentNames <- getNames(celSetMATNormalisedUniquePositions)
	annotationTable <- read.csv(paste("annotationData", "humanGenomeAnnotation.csv", sep = .Platform$file.sep))
	designMatrix <- matrix(c(-1, 1, 1, -1), ncol = 1, dimnames = list(experimentNames, c("K9Ac I.P. - Input")))
	results <- blocksStats(celSetMATNormalisedUniquePositions, annotationTable, design = designMatrix, upStream = 2000, downStream = 0)

	GSTP1row <- grep("GSTP1", results$symbol)
	KLK2row <- grep("KLK2", results$symbol)
	if(results[GSTP1row, "start"] != 67107856 || results[GSTP1row, "end"] != 67110699 || results[GSTP1row, "df.2000.0"] != 11 || round(results[GSTP1row, "pvals.K9Ac.I.P....Input"], 3) != 0.744 || results[KLK2row, "start"] != 56068500 || results[KLK2row, "end"] != 56075635 || results[KLK2row, "df.2000.0"] != 35 || round(results[KLK2row, 11], 3) != 0)
		stop("blocksStats giving unexpected results for acetylation arrays.")
	
	load(paste("rawData", "sequencing", "seq_data.Rdata", sep = .Platform$file.sep))
	annotationTable$position <- ifelse(annotationTable$strand == '+', annotationTable$start, annotationTable$end)
	results <- annotationCounts(rs, annotationTable, 1000, 1000, seqLen=300)
    if(!all(results[100,]==c(8,6,1,4)) || !all(results[5000,]==c(0,0,50,52))) stop("annotationCounts giving unexpected results.")

    results <- blocksStats(rs, annotationTable, cbind(test=c(1,1,-1,-1)), 1000, 1000, seqLen=300, verbose=FALSE)
    if (!all(round(results$FDR_test[1:10],3)==c(0.825,0.000,1.000,0.239,0.405,0.752,0.001,0.037,0.118,0.217))||!all(round(results$logFC_test[5000:5010],3)==c(-32.974,-0.701,-1.444,0.304,0.112,1.181,0.786,-1.025,-2.327,-0.018,-0.086)))
        stop("blocksStats giving unexpected results for sequence data.")

	units <- indexOf(cdfFile, "chrY")
	indices <- getCellIndices(cdfTableUniquePositions, units = units, stratifyBy = "pm", unlist = TRUE, useNames = FALSE)
    #set random number seed so results will be the same in every test
    set.seed(27) #my lucky number
    results <- regionStats(celSetMATNormalisedUniquePositions, designMatrix, ind = indices, probeWindow = 500, nPermutations = 10)$regions$`K9Ac I.P. - Input`
	bestRegion <- which(abs(results$score) == max(abs(results$score)))

	if(results[bestRegion, "start"] != 11829591 || results[bestRegion, "end"] != 11830892)
		stop("Error in regionStats function")	

	cat("regionStats tested fine.\n")
}

cat("All tests passed.\n")
setwd(userWd)
