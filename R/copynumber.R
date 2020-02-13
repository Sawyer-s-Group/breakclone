#' @import data.table
#' @import GenomicRanges
#' @import S4Vectors

getScoreCN <- function(segmentTable, pair, populationBreakpoints, cnType, maxgap){

  sample1 <- segmentTable[segmentTable$SampleID == pair[1],]
  sample2 <- segmentTable[segmentTable$SampleID == pair[2],]

  if(cnType == "alleleSpecific"){
    sample1$nTotal <- sample1$nMajor + sample1$nMinor
    sample2$nTotal <- sample2$nMajor + sample2$nMinor

    sample1 <- handlePloidy(sample1)
    sample2 <- handlePloidy(sample2)

    sample1$status <- ifelse(sample1$nTotal < 2, "loss", ifelse(sample1$nTotal == 2 & sample1$nMinor == 1, "norm", ifelse(sample1$nTotal == 2 & sample1$nMinor == 0, "cnloh", ifelse(sample1$nTotal > 4, "amp", "gain"))))
    sample2$status <- ifelse(sample2$nTotal < 2, "loss", ifelse(sample2$nTotal == 2 & sample2$nMinor == 1, "norm", ifelse(sample2$nTotal == 2 & sample2$nMinor == 0, "cnloh", ifelse(sample2$nTotal > 4, "amp", "gain"))))

    sample1 <- sample1[!"norm", on = "status"]
    sample2 <- sample2[!"norm", on = "status"]

    sample1_lists <- c(list(sample1["loss", on = "status"]), list(sample1["cnloh", on = "status"]), list(sample1["gain", on = "status"]), list(sample1["amp", on = "status"]), list(sample1["gain", on = "status"]), list(sample1["loss", on = "status"]))
    sample2_lists <- c(list(sample2["loss", on = "status"]), list(sample2["cnloh", on = "status"]), list(sample2["gain", on = "status"]), list(sample2["amp", on = "status"]), list(sample2["amp", on = "status"]), list(sample2["cnloh", on = "status"]))

  } else if(cnType == "VCF"){
    presentStates <- unique(c(sample1$SVType, sample2$SVType))

    sample1_lists <- sapply(presentStates, function(x){c(list(sample1[x, on = "SVType"]))})
    sample2_lists <- sapply(presentStates, function(x){c(list(sample2[x, on = "SVType"]))})
  }

  # sample1_lists <- c(list(), list(sample1["loss", on = "status"]), list(sample1["gain", on = "status"]))
  # sample2_lists <- c(list(), list(sample2["loss", on = "status"]), list(sample2["gain", on = "status"]))
  #

  # sample1_granges <- lapply(sample1_lists, makeGRangesFromDataFrame)
  # sample2_granges <- lapply(sample2_lists, makeGRangesFromDataFrame)

  #tryCatch creates an empty GRanges object if the list is empty - would error out otherwise
  sample1_granges <- lapply(sample1_lists, function(x){tryCatch({makeGRangesFromDataFrame(x)}, error = function(e){GRanges()})})
  sample2_granges <- lapply(sample2_lists, function(x){tryCatch({makeGRangesFromDataFrame(x)}, error = function(e){GRanges()})})

  hits_start <- mapply(function(x,y){x[queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = maxgap)))]}, sample1_granges, sample2_granges)
  hits_end <- mapply(function(x,y){x[queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = maxgap)))]}, sample1_granges, sample2_granges)

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

  # nonhits_start_1 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = 11449)))]}, sample1_granges, sample2_granges)
  # nonhits_end_1 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = 11449)))]}, sample1_granges, sample2_granges)
  #
  # nonhits_start_2 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "start", maxgap = 11449)))]}, sample2_granges, sample1_granges)
  # nonhits_end_2 <- mapply(function(x,y){x[-queryHits(suppressWarnings(findOverlaps(x,y, type = "end", maxgap = 11449)))]}, sample2_granges, sample1_granges)

  nconcordant_adj <- score_from_hits_start + score_from_hits_end
  #ndiscordant <- 2 * (nrow(sample1) + nrow(sample2)) - sum(unlist(lapply(c(hits_start, hits_end), length)))
  # Old calculation
  # ndiscordant <- nrow(sample1) + nrow(sample2) - 2 * nconcordant_adj
  # Pretty sure the below calculation makes more sense - 2 breakpoints per row! Small difference but still.
  ndiscordant <- 2 * (nrow(sample1) + nrow(sample2) - nconcordant_adj)

  score <- nconcordant_adj/(nconcordant_adj + 0.5 * ndiscordant)
  return(score)
}

handlePloidy <- function(sample){
  sample$segLen <- sample$End - sample$Start
  sample_states <- aggregate(segLen ~ nTotal, sample, sum)
  sample_states$adjLen <- as.numeric(sample_states$nTotal) * sample_states$segLen
  sample_ploidy <- sum(sample_states$adjLen) / sum(sample$segLen)

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

#' @export
calculateRelatednessCN <- function(segmentTable, pairs, reference = NULL, cnType = c("alleleSpecific", "VCF"), excludeChromosomes = "Y", maxgap = NULL){
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

#' #' @export
#' makeReference <- function(segmentTable, nperm = 10, cnType = c("alleleSpecific", "VCF"), excludeChromosomes = "Y", maxgap = NULL){
#'
#'   cnType <- match.arg(cnType)
#'   segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
#'   populationBreakpoints <- collatePopulationBreakpoints(segmentTable, cnType)
#'   if(is.null(maxgap)){maxgap <- calculateMaxGap(segmentTable, cnType)}
#'
#'   reference <- numeric()
#'   for(i in 1:nperm){
#'     message("Constructing reference: Iteration #", i)
#'     randomise <- sample(unique(segmentTable$SampleID))
#'     random_pairs <- cbind.data.frame(randomise[1:(length(randomise)/2)], randomise[(length(randomise)/2 + 1):length(randomise)])
#'     apply(random_pairs, 1, function(x){if(x[1] == x[2]){stop("yes, it's possible: ", x[1])}})
#'     pair_scores <- apply(random_pairs, 1, function(x){getScoreCN(segmentTable, as.character(x), populationBreakpoints, cnType, maxgap)})
#'     reference <- c(reference, pair_scores)
#'
#'
#'   }
#'   return(reference)
#' }
#'
#' #' @export
#' makeReferenceMixingPairs <- function(segmentTable, pairs, nperm = 10, cnType = c("alleleSpecific", "VCF"), excludeChromosomes = "Y", maxgap = NULL){
#'
#'   cnType <- match.arg(cnType)
#'   segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
#'   populationBreakpoints <- collatePopulationBreakpoints(segmentTable, cnType)
#'   if(is.null(maxgap)){maxgap <- calculateMaxGap(segmentTable, cnType)}
#'
#'   reference <- numeric()
#'   for(i in 1:nperm){
#'     message("Constructing reference: Iteration #", i)
#'
#'     random_pairs <- as.data.table(cbind(sample(pairs[[1]]), sample(pairs[[2]])))
#'     random_pairs <- random_pairs[!apply(random_pairs, 1, function(y){any(apply(pairs, 1, function(x){all(x == y)}))})]
#'
#'     pair_scores <- apply(random_pairs, 1, function(x){getScoreCN(segmentTable, as.character(x), populationBreakpoints, cnType, maxgap)})
#'     reference <- c(reference, pair_scores)
#'
#'
#'   }
#'   return(reference)
#' }

#' @export
makeReferenceCN <- function(segmentTable, pairs, patients = NULL, delimiter = "_", cnType = c("alleleSpecific", "VCF"), excludeChromosomes = "Y", maxgap = NULL){

  if(is.null(patients)){
    p1 <- sapply(strsplit(pairs$Sample1, delimiter), "[", 1)
    p2 <- sapply(strsplit(pairs$Sample2, delimiter), "[", 1)
    if(all(p1 == p2)){
      patients <- p1
    } else {
      stop("Autodetecting patient IDs failed!")
    }
  }
  patients <- rbind(as.data.table(cbind(patients, pairs$Sample1)), as.data.table(cbind(patients, pairs$Sample2)))
  colnames(patients) <- c("patient", "sample")
  patients <- unique(patients)
  setkey(patients, "sample")

  refPairs <- expand.grid(list(Sample1 = unique(pairs$Sample1), Sample2 = unique(pairs$Sample2)), stringsAsFactors = FALSE)
  refPairs <- refPairs[patients[refPairs$Sample1]$patient != patients[refPairs$Sample2]$patient,]

  message("Making reference based on ", nrow(refPairs), " possible pairs, this might take a while")

  cnType <- match.arg(cnType)
  segmentTable <- segmentTable[!excludeChromosomes, on = "Chr"]
  populationBreakpoints <- collatePopulationBreakpoints(segmentTable, cnType)
  if(is.null(maxgap)){maxgap <- calculateMaxGap(segmentTable, cnType)}

  reference <- apply(refPairs, 1, function(x){getScoreCN(segmentTable, as.character(x), populationBreakpoints, cnType, maxgap)})
  return(reference)
}