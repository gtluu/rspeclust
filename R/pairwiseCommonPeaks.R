#' Pairwise Peaks in Common
#' 
#' Compare two peak lists to identify peaks in common based on similarity
#' score.
#' 
#' @param x first peak list; single column \code{dataframes} with the
#' column name 'mz'.
#' @param y second peak list; single column \code{dataframes} with the
#' column name 'mz'.
#' @param tol \code{double} value to be used for m/z tolerance when merging
#' peak lists.
#' @param cutoff \code{double} value between 0 and 1; peaks with scores
#' greater than the cutoff are determined to be common peaks.
#' @return \code{dataframe} with each row representing shared peaks from
#' multiple peak lists based on m/z value along with a column for similarity
#' scores
#' @example 
#' 
#' peakList1 <- data.frame('mz'=c(615.3456, 489.6651, 375.1968))
#' peakList2 <- data.frame('mz'=c(615.3589, 453.3596, 357.9618))
#' 
#' commonPeaks <- pairwiseCommonPeaks(peakList1, peakList2, tol=0.2,
#'                                    cutoff=0.7)
#' 
#' @export
pairwiseCommonPeaks <- function(x, y, tol, cutoff) {
  # Merge peak lists.
  peaks <- fuzzyjoin::difference_full_join(x, y, by='mz', max_dist=tol)
  peaks[is.na(peaks)] <- 0
  
  # Calculate peak similarity scores.
  peaks$score <- 1 - VGAM::erf(abs(peaks$mz.x - peaks$mz.y) / (2 * tol))
  peaks <- peaks[which(peaks$score >= cutoff),]
  
  return(peaks)
}
