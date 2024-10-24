## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  if (!require("devtools"))
#      install.packages("devtools")
#  
#  devtools::install_github("Sawyer-s-Group/breakclone")

## -----------------------------------------------------------------------------
library(breakclone)
data(segmentTable_AS)
data(segmentTable)
data(mutationTable)

## ----eval=FALSE---------------------------------------------------------------
#  segmentTable <- readVCFCn("/path/to/data")

## -----------------------------------------------------------------------------
class(segmentTable)
head(segmentTable)

## ----eval=FALSE---------------------------------------------------------------
#  segmentTable_AS <- readAlleleSpecific("/path/to/data")

## -----------------------------------------------------------------------------
class(segmentTable_AS)
head(segmentTable_AS)

## ----eval=FALSE---------------------------------------------------------------
#  mutationTable <- readVCFMutations("/path/to/data")

## -----------------------------------------------------------------------------
class(mutationTable)
head(mutationTable)

## -----------------------------------------------------------------------------
pairs_cn <- inferPairs(segmentTable)
pairs_cnAS <- inferPairs(segmentTable_AS)
pairs_muts <- inferPairs(mutationTable)
head(pairs_muts)

## -----------------------------------------------------------------------------
set.seed(111)
reference <- makeReferenceCN(segmentTable=segmentTable, 
                             pairs=pairs_cn, 
                             cnType='VCF')

## -----------------------------------------------------------------------------
results <- calculateRelatednessCn(segmentTable=segmentTable, 
                                  pairs=pairs_cn, 
                                  reference=reference, 
                                  cnType='VCF')

summary <- summarizeClonalityCN(clonalityResults=results, 
                                segmentTable, 
                                thres_ambiguous=0.05, 
                                thres_related=0.01, 
                                cnType='VCF')

head(summary)

## -----------------------------------------------------------------------------
breaks <- getSharedBreaks(segmentTable=segmentTable, 
                          pairs=pairs_cn, 
                          cnType="VCF", 
                          save=FALSE)
breaks[[1]]

## -----------------------------------------------------------------------------
data(binnedTable)
data(cnTable)

plotCNpairVCF(binnedTable, 
              cnTable, 
              pair=c('Patient53_Primary', 'Patient53_Recurrence'), 
              segmentTable=segmentTable, 
              breaks=breaks, 
              build='hg38')

## ----warning=FALSE, message=FALSE---------------------------------------------
reference <- makeReferenceCN(segmentTable=segmentTable_AS, 
                             pairs=pairs_cnAS)

## -----------------------------------------------------------------------------
results <- calculateRelatednessCn(segmentTable=segmentTable_AS, 
                                  pairs=pairs_cnAS, 
                                  reference=reference)

summary <- summarizeClonalityCN(clonalityResults=results, 
                                segmentTable=segmentTable_AS)
head(summary)

## -----------------------------------------------------------------------------
breaks <- getSharedBreaks(segmentTable=segmentTable_AS, 
                          pairs=pairs_cnAS, 
                          save = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  plotCNpairalleleSpecific(ASCATobj=ascat.bc,
#                           segmentTable=segmentTable_AS,
#                           pair=c('Patient5_Primary', 'Patient5_Recurrence'),
#                           breaks,
#                           build = 'hg19')

## -----------------------------------------------------------------------------
data(brca)
head(brca)

## -----------------------------------------------------------------------------
reference <- makeReferenceMutations(mutationTable=mutationTable, 
                                    pairs=pairs_muts, 
                                    additionalMutations=brca, 
                                    nAdditionalSamples=346)

## -----------------------------------------------------------------------------
results <- calculateRelatednessMutations(mutationTable=mutationTable, 
                                         pairs=pairs_muts, 
                                         reference=reference, 
                                         additionalMutations=brca,
                                         nAdditionalSamples = 346)

summary <- summarizeClonalityMuts(clonalityResults=results,
                                  mutationTable=mutationTable)
head(summary)

## -----------------------------------------------------------------------------
shared_muts <- getSharedMuts(mutationTable=mutationTable, 
                             pairs=pairs_muts, 
                             save=FALSE)
shared_muts[[1]]

plotScatterVAF(mutationTable=mutationTable, 
               pair=c("Patient43_Primary", "Patient43_Recurrence"), 
               annotGenes=TRUE)

## -----------------------------------------------------------------------------
plotScoresDensity(reference, results)

## -----------------------------------------------------------------------------
plotScoresHistogram(reference, results)

## ----message=FALSE------------------------------------------------------------
plotSummary(summary, 
            extraAnno = rep('Invasive', nrow(summary)), 
            colorsExtraAnno = c("DCIS" = "#e6cde3", 
                                "Invasive" = "#af67a7")
            )

## -----------------------------------------------------------------------------
library(TCGAretriever)

brca <- get_mutation_data(case_list_id='brca_tcga_pub_complete', 
                          gprofile_id='brca_tcga_pub_mutations', 
                          glist=unique(mutationTable$annotation))
brca$chr <- sub("23", "X", brca$chr) # 23 to X
brca$chr <- paste0('chr', brca$chr)
brca <- brca[brca$startPosition == brca$endPosition,] # take only SNPs mutations
brca <- brca[brca$mutationType %in% c("Missense_Mutation","Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region", "Splice_Site", "Translation_Start_Site"),] # filter specific mutations

## ----eval=FALSE---------------------------------------------------------------
#  brca <- makeGRangesFromDataFrame(brca,
#                                   start.field = "startPosition",
#                                   end.field = "endPosition",
#                                   seqnames.field = "chr")

## ----eval=FALSE---------------------------------------------------------------
#  library(liftOver)
#  
#  download.file("https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz",
#                "hg19ToHg38.over.chain.gz")
#  R.utils::gunzip("hg19ToHg38.over.chain.gz", overwrite = TRUE)
#  
#  ch <- import.chain("hg19ToHg38.over.chain")
#  brca <- unlist(liftOver(brca, ch))

## -----------------------------------------------------------------------------
sessionInfo()

