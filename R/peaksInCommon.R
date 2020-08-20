#' SPECLUST Peaks in Common Workflow
#' 
#' Run the Peaks in Common workflow based on the SPECLUST online tool and
#' algorithms.
#' 
#' @param peakLists \code{list} of single column \code{dataframes} with the
#' column name 'peaks'.
#' @param sigma \code{double} value to be used for m/z tolerance when merging
#' peak lists.
#' @param pairwiseCutoff \code{double} value between 0 and 1; peaks with scores
#' greater than the cutoff are determined to be common peaks.
#' @param multipleCutoff \code{double} value between 0 and 1; peaks with scores
#' greater than the cutoff are determined to be common peaks.
#' @param consensusCutoff \code{numeric} value; peaks that appear in \code{n}
#' peak lists are kept if \code{n >= cutoff}.
#' @return \code{dataframe} with each row representing shared peaks from
#' multiple peak lists based on m/z value along with a column for number of
#' peak lists it appears in.
#' @example 
#' 
#' peakList1 <- data.frame('peaks'=c(615.3456, 489.6651, 375.1968))
#' peakList2 <- data.frame('peaks'=c(615.3589, 453.3596, 357.9618))
#' peakList3 <- data.frame('peaks'=c(615.3358, 861.3456, 198.3557))
#' 
#' peakLists <- list(peakList1, peakList2, peakList3)
#' 
#' consensus <- consensusPeakList(peakLists, tol=0.2, cutoff=0.5)
#' 
#' @export
peaksInCommon <- function(peakLists, sigma, pairwiseCutoff, multipleCutoff,
                          consensusCutoff) {
  # Get all possible combinations of peak lists.
  peakListCombos <- combinations(n=length(peakLists), r=2,
                                 v=seq(1, length(peakLists)))
  
  # Get dataframes with pairwise similarity scores for each combination of peak
  # lists.
  pairwiseResults <- list()
  for (i in 1:nrow(peakListCombos)) {
    i <- peakListCombos[i,]
    index <- length(pairwiseResults) + 1
    pairwiseResults[[index]] <- pairwiseCommonPeaks(peakLists[[i[1]]],
                                                    peakLists[[i[2]]],
                                                    tol=sigma,
                                                    cutoff=pairwiseCutoff)
  }
  
  # Get dataframe with multiple peak similarity scores for each peak.
  multipleResults <- multipleCommonPeaks(peakLists, tol=sigma,
                                         cutoff=multipleCutoff)
  
  # Get dataframe with consensus peak list.
  consensusResults <- consensus(peakLists, tol=sigma, cutoff=consensusCutoff)
  
  return(list('pairwise'=pairwiseResults,
              'multiple'=multipleResults,
              'consensus'=consensusResults))
}