#' @import data.table
#' @import GenomicRanges
#' @import S4Vectors

callalleleSpecificCN <- function(sample) {
  sample$nTotal <- sample$nMajor + sample$nMinor
  sample <- handlePloidy(sample)
  sample$SVType <- ifelse(sample$nTotal < 2, "loss", ifelse(sample$nTotal == 2 & sample$nMinor == 1, "norm", ifelse(sample$nTotal == 2 & sample$nMinor == 0, "cnloh", ifelse(sample$nTotal > 4, "amp", "gain"))))
  sample <- sample[!"norm", on = "SVType"]
  return(sample)
}

getHitCN <- function(sample1, sample2, pair, position=NULL, cnType = c("alleleSpecific", "VCF"), maxgap){
  cnType <- match.arg(cnType)
  
  if(cnType == "alleleSpecific"){
    sample1 <- callalleleSpecificCN(sample1)
    sample2 <- callalleleSpecificCN(sample2)
    
    sample1_lists <- c(list(sample1["loss", on = "SVType"]), list(sample1["cnloh", on = "SVType"]), list(sample1["gain", on = "SVType"]), list(sample1["amp", on = "SVType"]), list(sample1["gain", on = "SVType"]), list(sample1["loss", on = "SVType"]))
    sample2_lists <- c(list(sample2["loss", on = "SVType"]), list(sample2["cnloh", on = "SVType"]), list(sample2["gain", on = "SVType"]), list(sample2["amp", on = "SVType"]), list(sample2["amp", on = "SVType"]), list(sample2["cnloh", on = "SVType"]))
    
  } else if(cnType == "VCF"){
    presentStates <- unique(c(sample1$SVType, sample2$SVType))
    
    sample1_lists <- sapply(presentStates, function(x){c(list(sample1[x, on = "SVType"]))})
    sample2_lists <- sapply(presentStates, function(x){c(list(sample2[x, on = "SVType"]))})
  }
  
  #tryCatch creates an empty GRanges object if the list is empty - would error out otherwise
  sample1_granges <- lapply(sample1_lists, function(x){tryCatch({makeGRangesFromDataFrame(x)}, error = function(e){GRanges()})})
  sample2_granges <- lapply(sample2_lists, function(x){tryCatch({makeGRangesFromDataFrame(x)}, error = function(e){GRanges()})})
  
  hits <- mapply(function(x,y){x[queryHits(suppressWarnings(findOverlaps(x,y, type = position, maxgap = maxgap)))]}, sample1_granges, sample2_granges)
  return(hits)
}

getScoreCN <- function(segmentTable, pair, populationBreakpoints, cnType, maxgap){
  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]
  
  if(nrow(sample1) == 0 | nrow(sample2) == 0){
    if(nrow(sample1) == 0){
      warning("Sample ", pair[1], " seems to have no aberrations, check your data if that isn't expected")
    }
    if(nrow(sample2) == 0){
      warning("Sample ", pair[2], " seems to have no aberrations, check your data if that isn't expected")
    }
    return(0)
  }
  
  hits_start <- getHitCN(sample1, sample2, pair, position='start', cnType, maxgap)
  hits_end <- getHitCN(sample1, sample2, pair, position='end', cnType, maxgap)

  score_from_hits_start <- sum(unlist(lapply(hits_start, function(x){1 - countOverlaps(query = x, subject = populationBreakpoints$Starts, type = "start", maxgap = maxgap) / length(unique(segmentTable$SampleID))})))
  if(score_from_hits_start < 0){
    warning("Hit next probe! Consider lowering the maxgap")
    score_from_hits_start <- 0
  }
  score_from_hits_end <- sum(unlist(lapply(hits_end, function(x){1 - countOverlaps(query = x, subject = populationBreakpoints$Ends, type = "end", maxgap = maxgap) / length(unique(segmentTable$SampleID))})))
  if(score_from_hits_end < 0){
    warning("Hit next probe! Consider lowering the maxgap")
    score_from_hits_end <- 0
  }
  
  nconcordant_adj <- score_from_hits_start + score_from_hits_end
  
  ndiscordant <- 2 * (nrow(sample1) + nrow(sample2) - nconcordant_adj)
  
  score <- nconcordant_adj/(nconcordant_adj + 0.5 * ndiscordant)
  return(score)
}

calculatePloidy <- function(sample){
  sample$segLen <- sample$End - sample$Start
  sample_states <- aggregate(segLen ~ nTotal, sample, sum)
  sample_states$adjLen <- as.numeric(sample_states$nTotal) * sample_states$segLen
  sample_ploidy <- sum(sample_states$adjLen) / sum(sample$segLen)
  return(sample_ploidy)
}

handlePloidy <- function(sample){
  sample_ploidy <- calculatePloidy(sample)
  
  if(sample_ploidy >= 3.5){
    warning("Sample ", sample$SampleID[1], " looks WGD, inferring ancestral Cn state")
    sample$nTotal <- ifelse(sample$nMinor == (sample$nTotal / 2), 2, sample$nTotal)
    sample$nTotal <- ifelse(sample$nMinor == 0, 1, sample$nTotal)
    sample$nTotal <- ifelse(sample$nTotal >= 3 & sample$nTotal <= 4 & sample$nMinor == 1, 2, sample$nTotal)
    sample$nTotal <- ifelse(sample$nTotal >= 5 & sample$nTotal <= 10, 3, sample$nTotal)
  }
  return(sample)
}

collatePopulationBreakpoints <- function(segmentTable, cnType){
  if(cnType == "alleleSpecific"){
    segmentTable <- segmentTable[!(segmentTable$nMajor == 1 & segmentTable$nMinor == 1),]
  }
  
  populationStarts <- makeGRangesFromDataFrame(segmentTable[,c("Chr", "Start")], start.field = "Start", end.field = "Start")
  populationEnds <- makeGRangesFromDataFrame(segmentTable[,c("Chr", "End")], start.field = "End", end.field = "End")
  
  populationBreakpoints <- c(list(populationStarts), list(populationEnds))
  names(populationBreakpoints) <- c("Starts", "Ends")
  return(populationBreakpoints)
}

calculateMaxGap <- function(segmentTable, cnType){
  if(cnType == "alleleSpecific"){
    avgProbeDistance <- mean((segmentTable$End - segmentTable$Start) / segmentTable$nProbes)
    maxgap <- 5 * avgProbeDistance
  } else if(cnType == "VCF"){
    possBinSizes <- unique(segmentTable$Length / segmentTable$Bins)
    binSize <- possBinSizes[which.max(tabulate(match(segmentTable$Length / segmentTable$Bins, possBinSizes)))]
    maxgap <- 5 * binSize
  }
  return(maxgap)
}

#' Calculate relatedness scores for paired tumours
#'
#' Calculates the relatedness scores and (optionally) p-values for paired tumours from copy number data
#' @param segmentTable A segment table generated by the readAlleleSpecific or readVCFCn functions.
#' @param pairs A table of paired samples from the dataset, to test for relatedness.
#' @param reference A numeric vector of pair scores comprising the reference distribution, generated from the \code{makeReferenceCN} function. If omitted, p-value calculation will be skipped.
#' @param cnType The type of copy number data provided. Currently supported options are a custom allele specific data format, and standard copy number VCF files.
#' @param excludeChromosomes The name(s) of any chromosomes to be excluded.
#' @param maxgap The maximum gap between two breakpoints for them to be considered concordant. If unspecified, it is automatically set to 5 times the average interprobe distance of the assay.
#' @return A data frame listing the tumour pairs contained in \code{pairs}, their relatedness scores and p-values for relatedness.
#' @export
calculateRelatednessCn <- function(segmentTable, pairs, reference = NULL, cnType = c("alleleSpecific", "VCF"), excludeChromosomes = "Y", maxgap = NULL){
  cnType <- match.arg(cnType)
  segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
  populationBreakpoints <- collatePopulationBreakpoints(segmentTable, cnType)
  if(is.null(maxgap)){maxgap <- calculateMaxGap(segmentTable, cnType)}
  
  pair_scores <- apply(pairs, 1, function(x){getScoreCN(segmentTable, as.character(x), populationBreakpoints, cnType, maxgap)})
  
  if(is.null(reference)){warning("No reference supplied, p-values not calculated", immediate. = TRUE)}
  pair_ps <- unlist(lapply(pair_scores, function(x){mean(x <= reference)}))
  results <- cbind.data.frame(pairs, pair_scores, pair_ps)
  
  return(results)
}

#' Generate reference distribution from copy number data
#'
#' Generates the reference distribution of concordance scores from unpaired tumours for a given dataset.
#'
#' The default is to assume one tumour pair per patient. If the sample is found in more than one pair, indicating more than one pair per patient, an error will be thrown.
#' In that case you must specify either \code{patients} or \code{delimiter}.
#'
#' @param segmentTable A segment table generated by the readAlleleSpecific or readVCFCn functions.
#' @param pairs A table of paired samples from the dataset. All tumours present in this table will be paired with all tumours from other patients.
#' @param patients A character vector of patient IDs, parallel to the pairs table, used to prevent tumours originating from the same patient from being used in the reference distribution (optional).
#' @param delimiter A character separating patient IDs from tumour-specific identifiers in the sample IDs. Ignored if \code{patients} is provided (optional).
#' @param cnType The type of copy number data provided. Currently supported options are a custom allele specific data format, and standard copy number VCF files.
#' @param excludeChromosomes The name(s) of any chromosomes to be excluded.
#' @param maxgap The maximum gap between two breakpoints for them to be considered concordant. If unspecified, it is automatically set to 5 times the average interprobe distance of the assay.
#' @return A numeric vector of pair scores comprising the reference distribution.
#' @export
makeReferenceCN <- function(segmentTable, pairs, patients = NULL, delimiter = NULL, cnType = c("alleleSpecific", "VCF"), excludeChromosomes = "Y", maxgap = NULL){
  if(is.null(patients) & is.null(delimiter)){
    patients <- as.character(seq(1, nrow(pairs)))
  } else if(is.null(patients)){
    p1 <- sapply(strsplit(pairs$Sample1, delimiter), "[", 1)
    p2 <- sapply(strsplit(pairs$Sample2, delimiter), "[", 1)
    if(all(p1 == p2)){
      patients <- p1
    } else {
      stop("Autodetecting patient IDs failed!")
    }
  }
  patients <- rbind(as.data.table(cbind(patients, pairs[[1]])), as.data.table(cbind(patients, pairs[[2]])))
  colnames(patients) <- c("patient", "sample")
  patients <- unique(patients)
  setkey(patients, "sample")
  
  refPairs <- expand.grid(list(Sample1 = unique(pairs[[1]]), Sample2 = unique(pairs[[2]])), stringsAsFactors = FALSE)
  refPairs <- refPairs[patients[refPairs$Sample1]$patient != patients[refPairs$Sample2]$patient,]
  
  message("Making reference based on ", nrow(refPairs), " possible pairs, this might take a while")
  
  cnType <- match.arg(cnType)
  segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
  populationBreakpoints <- collatePopulationBreakpoints(segmentTable, cnType)
  if(is.null(maxgap)){maxgap <- calculateMaxGap(segmentTable, cnType)}
  
  reference <- apply(refPairs, 1, function(x){getScoreCN(segmentTable, as.character(x), populationBreakpoints, cnType, maxgap)})
  return(reference)
}

exportSharedBreaks <- function(pair, segmentTable, outdir = '.', cnType = c("alleleSpecific", "VCF"), maxgap, save){
  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]
  
  hits_start <- getHitCN(sample1, sample2, pair, position='start', cnType, maxgap)
  hits_end <- getHitCN(sample1, sample2, pair, position='end', cnType, maxgap)
  
  start <- do.call("rbind", lapply(hits_start, function(x) as.data.frame(x)))[,-c(4,5)]
  end <- do.call("rbind", lapply(hits_end, function(x) as.data.frame(x)))[,-c(4,5)]

  if (nrow(start) == 0 & nrow(end) == 0){
    shared_breaks <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("seqnames", "start", "end", "break_shared"))
  } else if (nrow(start) == 0) {
    end$break_shared <- 'end'
    shared_breaks <- end
  } else if (nrow(end)==0) {
    start$break_shared <- 'start'
    shared_breaks <- start
  } else {
    start$break_shared <- 'start'
    end$break_shared <- 'end'
    shared_breaks <- rbind(start, end)
  }
  
  if(isTRUE(save)){
    print(paste0(pair[1], "-", pair[2],  " pair doesn't contain shared breakpoints"))
    write.table(shared_breaks, paste0(outdir, "/", pair[1],"-",pair[2], ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  return(shared_breaks)
}


#' Generate the shared breakpoints per pair calculated by breakclone and used to assess relatedness.
#' 
#' @param segmentTable A segment table generated by the readAlleleSpecific or readVCFCn functions.
#' @param pairs A table of paired samples from the dataset. All tumours present in this table will be paired with all tumours from other patients.
#' @param outdir A path to save the shared breakpoints. If unspecified, it is automatically set to the working directory.
#' @param cnType The type of copy number data provided. Currently supported options are a custom allele specific data format, and standard copy number VCF files.
#' @param maxgap The maximum gap between two breakpoints for them to be considered concordant. If unspecified, it is automatically set to 5 times the average interprobe distance of the assay.
#' @param save TRUE if want to export the shared mutations.
#' @return List of shared breakpoints per pair It will generate a tsv file per pair in outdir with shared breakpoints. 
#' @export
getSharedBreaks <- function(segmentTable, pairs, outdir = '.', cnType = c("alleleSpecific", "VCF"), maxgap=NULL, save = FALSE){
  cnType <- match.arg(cnType)
  if(is.null(maxgap)){maxgap <- calculateMaxGap(segmentTable, cnType)}
  res <- apply(pairs, 1, function(x){exportSharedBreaks(as.character(x), segmentTable, outdir, cnType, maxgap, save)})
  names(res) <- paste0(pairs[[1]], '-', pairs[[2]])
  return(res)
}

#' Summarize clonality results of copy number aberration data.
#' 
#' @param clonalityResults A data frame listing the tumour pairs contained in \code{pairs}, their relatedness scores and p-values for relatedness.
#' @param segmentTable A segment table generated by the readAlleleSpecific or readVCFCn functions.
#' @param thres_ambiguous P-value to define ambiguous. thres_ambiguous=0.05 by default.
#' @param thres_related P-value to define related thres_related=0.05 by default.
#' @param cnType The type of copy number data provided. Currently supported options are a custom allele specific data format, and standard copy number VCF files.
#' @param maxgap The maximum gap between two breakpoints for them to be considered concordant. If unspecified, it is automatically set to 5 times the average interprobe distance of the assay.
#' @return A data frame listing the clonality results \code{clonalityResults}, clonality verdict based on the thresholds and number and fraction of shared and private breakpoints per sample.
#' @export
summarizeClonalityCN <- function(clonalityResults, segmentTable, thres_ambiguous = 0.05, thres_related = 0.01, cnType = c('alleleSpecific', 'VCF'), maxgap = NULL){
  cnType <- match.arg(cnType)
  if(is.null(maxgap)){maxgap <- calculateMaxGap(segmentTable, cnType)}
  
  clonalityResults$verdict <- factor(ifelse(clonalityResults$pair_ps<=thres_related, 'Related', ifelse(clonalityResults$pair_ps<=thres_ambiguous, 'Ambiguous', 'Unrelated')), levels = c("Related", "Ambiguous", "Unrelated"))
  clonalityResults$shared <- as.numeric(lapply(apply(clonalityResults[,1:2], 1, function(x){exportSharedBreaks(pair=as.character(x), segmentTable=segmentTable, cnType=cnType, maxgap=maxgap, save=FALSE)}), nrow))
  clonalityResults$private_sample1 <- apply(clonalityResults, 1, function(x) nrow(segmentTable[segmentTable$SampleID==x[[1]],])*2) - clonalityResults[['shared']] 
  clonalityResults$private_sample2 <- apply(clonalityResults, 1, function(x) nrow(segmentTable[segmentTable$SampleID==x[[2]],])*2) - clonalityResults[['shared']]  
  clonalityResults$total <- clonalityResults$private_sample1+clonalityResults$shared+clonalityResults$private_sample2
  clonalityResults$fraction_private_sample1 <- clonalityResults$private_sample1/clonalityResults$total
  clonalityResults$fraction_private_sample2 <- clonalityResults$private_sample2/clonalityResults$total
  clonalityResults$fraction_shared <- clonalityResults$shared/clonalityResults$total
  return(clonalityResults)
}
