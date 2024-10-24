
# breakclone

<!-- badges: start -->
<!-- badges: end -->

Breakclone is a package for assessing clonal relationships of tumour pairs based on copy number or mutation data. Further details in [the main vignette](doc/breakclone.html), 


## Example

This is a basic example of usage with copy number data:

``` r
library(breakclone)
table <- readAlleleSpecific("/path/to/data")
pairs <- data.table::fread("/path/to/sample/sheet")
reference <- makeReferenceCN(table, pairs)
results <- calculateRelatednessCn(table, pairs, reference)
plotScoresDensity(reference, results)
```

This is a very similar example with mutation data:

``` r
library(breakclone)
table <- readVCFMutations("/path/to/data")
pairs <- data.table::fread("/path/to/sample/sheet")
reference <- makeReferenceMutations(table, pairs)
results <- calculateRelatednessMutations(table, pairs, reference)
plotScoresDensity(reference, results)
```
