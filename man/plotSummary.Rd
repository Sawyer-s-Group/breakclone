% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotSummary}
\alias{plotSummary}
\title{Plot stacked barplot of total, shared and non-shared breakpoints/mutations and breakclone score per sample.}
\usage{
plotSummary(
  summary,
  sortBy = c("verdict", "fraction_shared"),
  colors = c("#3987bb", "#aee7ea", "#0c0d0c"),
  colorClonality = c("#f5c61a", "#fb4f14", "#660c21"),
  colorScore = "#b48c9c",
  colorN = "grey50",
  patients = NULL,
  delimiter = NULL,
  extraAnno = NULL,
  colorsExtraAnno = NULL,
  sample1 = "Primary",
  sample2 = "Recurrence"
)
}
\arguments{
\item{summary}{A data frame listing the clonality results \code{clonalityResults}, clonality verdict based on the thresholds and number and fraction of shared and private breakpoints per sample.}

\item{sortBy}{Stacked barplot will be sorted by these columns in summary. By default: verdict and fraction_shared.}

\item{colors}{A character vector with colors for private in sample 1, in sample 2 and shared breakpoints/mutations.}

\item{colorClonality}{A character vector with colors for related, ambiguous and unrelated samples.}

\item{colorScore}{Color of breakclone score bars.}

\item{colorN}{Color for number of total breakpoints/mutations bars.}

\item{patients}{A character vector of patient IDs, parallel to the pairs table, used to prevent tumours originating from the same patient from being used in the reference distribution (optional).}

\item{delimiter}{A character separating patient IDs from tumour-specific identifiers in the sample IDs. Ignored if \code{patients} is provided (optional).}

\item{extraAnno}{A character vector with annotation of patients. e.g. Invasive or DCIS.}

\item{colorsExtraAnno}{A named color vector for the extra pheno bar.}

\item{sample1}{Label for samples 1.}

\item{sample2}{Label for samples 1.}
}
\value{
Summary plot.
}
\description{
Plot stacked barplot of total, shared and non-shared breakpoints/mutations and breakclone score per sample.
}
