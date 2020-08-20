#' Multiple Peaks in Common
#' 
#' Compare multiple peak lists to identify peaks in common based on similarity
#' score.
#' 
#' @param peakLists \code{list} of single column \code{dataframes} with the
#' column name 'peaks'.
#' @param tol \code{double} value to be used for m/z tolerance when merging
#' peak lists.
#' @param cutoff \code{double} value between 0 and 1; peaks with scores
#' greater than the cutoff are determined to be common peaks.
#' @return \code{dataframe} with each row representing shared peaks from
#' multiple peak lists based on m/z value along with a column for similarity
#' scores
#' @example 
#' 
#' peakList1 <- data.frame('peaks'=c(615.3456, 489.6651, 375.1968))
#' peakList2 <- data.frame('peaks'=c(615.3589, 453.3596, 357.9618))
#' peakList3 <- data.frame('peaks'=c(615.3358, 861.3456, 198.3557))
#' 
#' peakLists <- list(peakList1, peakList2, peakList3)
#' 
#' commonPeaks <- multipleCommonPeaks(peakLists, tol=0.2, cutoff=0.5)
#' 
#' @export
multipleCommonPeaks <- function(peakLists, tol, cutoff) {
  # Combine peak lists into single large dataframe.
  bigDf <- combinePeakLists(peakLists, tol=tol)
  
  # Calculate pairwise similarity scores for each combination of peak lists.
  listOfScores <- list()
  for (i in 1:ncol(bigDf)) {
    scoresDf <- data.frame(matrix(NA, nrow=nrow(bigDf), ncol=0))
    for (j in 1:ncol(bigDf)) {
      if (i != j) {
        scoresDf[, as.character(j)] <- 1 - erf(abs(bigDf[, i] - bigDf[, j]) / (2 * tol))
        scoresDf[is.na(scoresDf)] <- 0
        scoresDf$average <- rowMeans(scoresDf)
      }
    }
    listOfScores[[i]] <- scoresDf
  }
  
  # Average pairwise similarity scores to obtain multiple similarity scores.
  scoresDf <- data.frame(matrix(NA, nrow=nrow(bigDf), ncol=length(listOfScores)))
  colnames(scoresDf) <- seq(1, length(listOfScores))
  for (i in 1:length(listOfScores)) {
    scoresDf[, i] <- listOfScores[[i]]$average
  }
  
  bigDf$score <- rowMeans(scoresDf)
  bigDf <- bigDf[which(bigDf$score >= cutoff),]
  
  return(bigDf)
}
