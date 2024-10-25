#' mutationTable
#'
#' List of mutation in a set of DCIS primaries and recurrences tumours, produced by sequencing. 
#'
#' @format A data table with 5 columns:
#' \describe{
#'   \item{SampleID}{sample ID column}
#'   \item{Chr}{chromosome name column}
#'   \item{Pos}{coordinate column}
#'   \item{AF}{allele frequency column}
#'   \item{annotation}{mutation details column}
#' }
#' @author Maria Roman Escorza
#' \email{maria.roman-escorza@@kcl.ac.uk}
"mutationTable"

#' segmentTable
#'
#' Copy number aberrations in a set of DCIS primaries and recurrences tumours, produced by low-coverage Whole Genome Sequencing.
#'
#' @format A data table with 7 columns:
#' \describe{
#'   \item{SampleID}{sample ID column}
#'   \item{Chr}{chromosome name column}
#'   \item{Start}{start coordinate column}
#'   \item{End}{end coordinate column}
#'   \item{Bins}{per-segment number of probes column}
#'   \item{SVType}{type of CNA column. Either DUP for region of elevated copy number relative to the reference or DEL for deletion relative to the reference}
#'   \item{Length}{lenght of the segment column}
#' }
#' @author Maria Roman Escorza
#' \email{maria.roman-escorza@@kcl.ac.uk}
"segmentTable"


#' segmentTable_AS
#'
#' Allele-specific copy number aberrations in a set of DCIS primaries and recurrences tumours, produced by SNParray.
#'
#' @format A data table with 7 columns:
#' \describe{
#'   \item{SampleID}{sample ID column}
#'   \item{Chr}{chromosome name column}
#'   \item{Start}{start coordinate column}
#'   \item{End}{end coordinate column}
#'   \item{nProbes}{per-segment number of probes column}
#'   \item{nMajor}{major copy number }
#'   \item{nMinor}{minor copy number }
#' }
#' @author Maria Roman Escorza
#' \email{maria.roman-escorza@@kcl.ac.uk}
"segmentTable_AS"

#' binnedTable
#'
#' A matrix of segment mean copy number values per bin or probe x sample
#'
#' @format A matrix
#' @author Maria Roman Escorza
#' \email{maria.roman-escorza@@kcl.ac.uk}
"binnedTable"


#' cnTable
#'
#' A matrix of relative copy number values per bin or probe x sample
#'
#' @format A matrix
#' @author Maria Roman Escorza
#' \email{maria.roman-escorza@@kcl.ac.uk}
"cnTable"
