#' Get Consensus Peak List
#' 
#' Generate a consensus peak list with average m/z values from multiple peak
#' lists.
#' 
#' @param peakLists \code{list} of single column \code{dataframes} with the
#' column name 'mz'.
#' @param tol \code{double} value to be used for m/z tolerance when merging
#' peak lists.
#' @param cutoff \code{numeric} value; peaks that appear in \code{n} peak lists
#' are kept if \code{n >= cutoff}.
#' @return \code{dataframe} with each row representing shared peaks from
#' multiple peak lists based on m/z value along with a column for number of
#' peak lists it appears in.
#' @example 
#' 
#' peakList1 <- data.frame('mz'=c(615.3456, 489.6651, 375.1968))
#' peakList2 <- data.frame('mz'=c(615.3589, 453.3596, 357.9618))
#' peakList3 <- data.frame('mz'=c(615.3358, 861.3456, 198.3557))
#' 
#' peakLists <- list(peakList1, peakList2, peakList3)
#' 
#' consensus <- consensusPeakList(peakLists, tol=0.2, cutoff=0.5)
#' 
#' @export
consensusPeakList <- function(peakLists, tol, cutoff) {
  # Combine peak lists into single large dataframe.
  bigDf <- combinePeakLists(peakLists, tol=tol)
  
  # Find how many peak lists each peak was found in and remove if lower than
  # the cutoff.
  bigDf$N <- rowSums(!is.na(bigDf))
  bigDf <- bigDf[which(bigDf$N >= cutoff),]
  
  # Get average m/z values.
  n <- bigDf$N
  bigDf <- subset(bigDf, select=-c(N))
  bigDf$average <- rowMeans(bigDf, na.rm=TRUE)
  bigDf$N <- n
  
  return(bigDf)
}
