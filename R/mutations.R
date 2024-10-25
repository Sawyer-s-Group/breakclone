getHitMut <- function(sample1, sample2, pair, scaleAFs) {
  sample1_granges <- makeGRangesFromDataFrame(sample1, start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)
  sample2_granges <- makeGRangesFromDataFrame(sample2, start.field = "Pos", end.field = "Pos", keep.extra.columns = TRUE)

  if (scaleAFs) {
    sample1_granges$AF <- sample1_granges$AF / max(sample1_granges$AF)
    sample2_granges$AF <- sample2_granges$AF / max(sample2_granges$AF)
  }

  overlaps <- suppressWarnings(findOverlaps(sample1_granges, sample2_granges))

  hits_sample1 <- sample1_granges[queryHits(overlaps)]
  hits_sample2 <- sample2_granges[subjectHits(overlaps)]
  if (length(overlaps) > 0) {
    nonhits_sample1 <- sample1_granges[-queryHits(overlaps)]
    nonhits_sample2 <- sample2_granges[-subjectHits(overlaps)]
  } else {
    nonhits_sample1 <- sample1_granges
    nonhits_sample2 <- sample2_granges
  }

  return(list(hits_sample1, hits_sample2, nonhits_sample1, nonhits_sample2))
}

getScoreMutations <- function(mutationTable, pair, populationMutations, nAdditionalSamples = 0, scaleAFs) {
  sample1 <- mutationTable[mutationTable$SampleID == pair[1], ]
  sample2 <- mutationTable[mutationTable$SampleID == pair[2], ]

  if (nrow(sample1) == 0 | nrow(sample2) == 0) {
    if (nrow(sample1) == 0) {
      warning("Sample ", pair[1], " seems to have no mutations, check your data if that isn't expected")
    }
    if (nrow(sample2) == 0) {
      warning("Sample ", pair[2], " seems to have no mutations, check your data if that isn't expected")
    }
    return(0)
  }

  hits <- getHitMut(sample1, sample2, pair, scaleAFs)
  hits_sample1 <- hits[[1]]
  hits_sample2 <- hits[[2]]
  nonhits_sample1 <- hits[[3]]
  nonhits_sample2 <- hits[[4]]

  nSamples <- length(unique(mutationTable$SampleID)) + nAdditionalSamples

  score <- sum(
    (hits_sample1$AF + hits_sample2$AF) / sqrt(
      countOverlaps(hits_sample1, populationMutations) / nSamples
    )
  ) / (sum(
    (hits_sample1$AF + hits_sample2$AF) / sqrt(
      countOverlaps(hits_sample1, populationMutations) / nSamples
    )
  ) +
    (0.5 * sum(
      (nonhits_sample1$AF / sqrt(
        countOverlaps(nonhits_sample1, populationMutations) / nSamples
      )
      ),
      (nonhits_sample2$AF / sqrt(
        countOverlaps(nonhits_sample2, populationMutations) / nSamples
      )
      )
    )
    )
  )


  return(score)
}

#' Calculate relatedness scores for paired tumours
#'
#' Calculates the relatedness scores and (optionally) p-values for paired tumours from mutation data
#' @param mutationTable A table of mutations in each sample and their allele frequencies.
#' @param pairs A table of paired samples from the dataset, to test for relatedness.
#' @param additionalMutations A table of mutations to be taken into account when calculating population frequencies. At a minimum, a table of the mutations in the population being studied. This is more informative when tumour type-specific mutations are included from external sources (e.g. TCGA).
#' @param nAdditionalSamples The number of samples used to derive the additional mutations table.
#' @param reference A numeric vector of pair scores comprising the reference distribution, generated from the \code{makeReferenceMutations} function. If omitted, p-value calculation will be skipped.
#' @param excludeChromosomes The name(s) of any chromosomes to be excluded.
#' @param scaleAFs Scale AFs per-sample by the highest AF within each sample. Only recommended for data with significant normal contamination that you are confident contains at least one clonal mutation per sample.
#' @author argymeg
#' @return A data frame listing the tumour pairs contained in \code{pairs}, their relatedness scores and p-values for relatedness.
#' @export
calculateRelatednessMutations <- function(mutationTable, pairs, additionalMutations = NULL, nAdditionalSamples = 0, reference = NULL, excludeChromosomes = "chrY", scaleAFs = FALSE) {
  mutationTable <- mutationTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(mutationTable)
  if (!is.null(additionalMutations)) {
    populationMutations <- c(populationMutations, additionalMutations)
  }
  pair_scores <- apply(pairs, 1, function(x) {
    getScoreMutations(mutationTable, as.character(x), populationMutations, nAdditionalSamples, scaleAFs)
  })

  if (is.null(reference)) {
    warning("No reference supplied, p-values not calculated", immediate. = TRUE)
  }
  pair_ps <- unlist(lapply(pair_scores, function(x) {
    mean(x <= reference)
  }))
  results <- cbind.data.frame(pairs, pair_scores, pair_ps)

  return(results)
}

#' Generate reference distribution from mutation data
#'
#' Generates the reference distribution of concordance scores from unpaired tumours for a given dataset.
#' @param mutationTable A table of mutations in each sample and their allele frequencies.
#' @param pairs A table of paired samples from the dataset. All tumours present in this table will be paired with all tumours from other patients.
#' @param patients A character vector of patient IDs, parallel to the pairs table, used to prevent tumours originating from the same patient from being used in the reference distribution (optional)
#' @param delimiter A character separating patient IDs from tumour-specific identifiers in the sample IDs. Ignored if \code{patients} is provided.
#' @param additionalMutations A table of mutations to be taken into account when calculating population frequencies. At a minimum, a table of the mutations in the population being studied. This is more informative when tumour type-specific mutations are included from external sources (e.g. TCGA).
#' @param nAdditionalSamples The number of samples used to derive the additional mutations table.
#' @param excludeChromosomes The name(s) of any chromosomes to be excluded.
#' @param scaleAFs Scale AFs per-sample by the highest AF within each sample. Only recommended for data with significant normal contamination that you are confident contains at least one clonal mutation per sample.
#' @author argymeg
#' @return A numeric vector of pair scores comprising the reference distribution.
#' @export
makeReferenceMutations <- function(mutationTable, pairs, patients = NULL, delimiter = NULL, additionalMutations = NULL, nAdditionalSamples = 0, excludeChromosomes = "Y", scaleAFs = FALSE) {
  if (is.null(patients) & is.null(delimiter)) {
    patients <- as.character(seq(1, nrow(pairs)))
  } else if (is.null(patients)) {
    p1 <- sapply(strsplit(pairs$Sample1, delimiter), "[", 1)
    p2 <- sapply(strsplit(pairs$Sample2, delimiter), "[", 1)
    if (all(p1 == p2)) {
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
  refPairs <- refPairs[patients[refPairs$Sample1]$patient != patients[refPairs$Sample2]$patient, ]

  message("Making reference based on ", nrow(refPairs), " possible pairs, this might take a while")

  mutationTable <- mutationTable[!excludeChromosomes, on = "Chr"]
  populationMutations <- collatePopulationMutations(mutationTable)
  if (!is.null(additionalMutations)) {
    populationMutations <- c(populationMutations, additionalMutations)
  }

  reference <- apply(refPairs, 1, function(x) {
    getScoreMutations(mutationTable, as.character(x), populationMutations, nAdditionalSamples, scaleAFs)
  })


  return(reference)
}


collatePopulationMutations <- function(mutationTable) {
  populationMutations <- makeGRangesFromDataFrame(mutationTable[, c("Chr", "Pos")], start.field = "Pos", end.field = "Pos")
  return(populationMutations)
}

exportSharedMuts <- function(pair, mutationTable, outdir = ".", scaleAFs, save) {
  sample1 <- mutationTable[mutationTable$SampleID == pair[1], ]
  sample2 <- mutationTable[mutationTable$SampleID == pair[2], ]

  if (nrow(sample1) == 0 | nrow(sample2) == 0) {
    if (nrow(sample1) == 0) {
      warning("Sample ", pair[1], " seems to have no mutations, check your data if that isn't expected")
    }
    if (nrow(sample2) == 0) {
      warning("Sample ", pair[2], " seems to have no mutations, check your data if that isn't expected")
    }
    shared_muts <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("chr", "pos", "AF_sample1", "AF_sample2"))
  } else {
    hits <- getHitMut(sample1, sample2, pair, scaleAFs)
    shared_muts <- data.frame(chr = as.data.frame(hits[[1]])$seqnames, pos = as.data.frame(hits[[1]])$start, AF_sample1 = as.data.frame(hits[[1]])$AF, AF_sample2 = as.data.frame(hits[[2]])$AF)
  }

  if (isTRUE(save)) {
    write.table(shared_muts, paste0(outdir, "/", pair[1], "-", pair[2], ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    print(paste0("Shared mutations ", pair[1], "-", pair[2], " saved in ", outdir))
  }

  return(shared_muts)
}


#' Get common mutations
#'
#' @param mutationTable A segment table generated by the readAlleleSpecific or readVCFCn functions.
#' @param pairs A table of paired samples from the dataset. All tumours present in this table will be paired with all tumours from other patients.
#' @param outdir A path to save the shared mutations If unspecified, it is automatically set to the working directory.
#' @param scaleAFs Scale AFs per-sample by the highest AF within each sample. Only recommended for data with significant normal contamination that you are confident contains at least one clonal mutation per sample.
#' @param save TRUE if want to export the shared mutations.
#' @return List of shared mutations per pair It will generate a tsv file per pair in outdir with shared mutations. ]
#' @author Maria Roman Escorza
#' \email{maria.roman-escorza@@kcl.ac.uk}
#' @export
getSharedMuts <- function(mutationTable, pairs, outdir = ".", scaleAFs = FALSE, save = FALSE) {
  res <- apply(pairs, 1, function(x) {
    exportSharedMuts(as.character(x), mutationTable, outdir, scaleAFs, save)
  })
  names(res) <- paste0(pairs[[1]], "-", pairs[[2]])
  return(res)
}

#' Summarize clonality results of mutation data.
#'
#' @param clonalityResults A data frame listing the tumour pairs contained in \code{pairs}, their relatedness scores and p-values for relatedness.
#' @param mutationTable A table of mutations in each sample and their allele frequencies.
#' @param thres_ambiguous P-value to define ambiguous. thres_ambiguous=0.05 by default.
#' @param thres_related P-value to define related thres_related=0.05 by default.
#' @param scaleAFs Scale AFs per-sample by the highest AF within each sample. Only recommended for data with significant normal contamination that you are confident contains at least one clonal mutation per sample.
#' @author Maria Roman Escorza
#' \email{maria.roman-escorza@@kcl.ac.uk}
#' @return A data frame listing the clonality results \code{clonalityResults}, clonality verdict based on the thresholds and number and fraction of shared and private mutations per sample.
#' @export
summarizeClonalityMuts <- function(clonalityResults, mutationTable, thres_ambiguous = 0.05, thres_related = 0.01, scaleAFs = FALSE) {
  clonalityResults$verdict <- factor(ifelse(clonalityResults$pair_ps <= thres_related, "Related", ifelse(clonalityResults$pair_ps <= thres_ambiguous, "Ambiguous", "Unrelated")), levels = c("Related", "Ambiguous", "Unrelated"))
  clonalityResults$shared <- as.numeric(lapply(apply(clonalityResults[, 1:2], 1, function(x) {
    exportSharedMuts(as.character(x), mutationTable, outdir, scaleAFs, save = FALSE)
  }), nrow))
  clonalityResults$private_sample1 <- apply(clonalityResults[, 1:2], 1, function(x) nrow(mutationTable[mutationTable$SampleID == x[[1]], ])) - clonalityResults[["shared"]]
  clonalityResults$private_sample2 <- apply(clonalityResults[, 1:2], 1, function(x) nrow(mutationTable[mutationTable$SampleID == x[[2]], ])) - clonalityResults[["shared"]]
  clonalityResults$total <- clonalityResults$private_sample1 + clonalityResults$shared + clonalityResults$private_sample2
  clonalityResults$fraction_private_sample1 <- clonalityResults$private_sample1 / clonalityResults$total
  clonalityResults$fraction_private_sample2 <- clonalityResults$private_sample2 / clonalityResults$total
  clonalityResults$fraction_shared <- clonalityResults$shared / clonalityResults$total
  return(clonalityResults)
}
