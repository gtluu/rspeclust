#' Combine Peak Lists
#' 
#' Combine multiple peak lists into a single \code{dataframe} by iteratively
#' merging each list.
#' 
#' @param peakLists \code{list} of single column \code{dataframes} with the
#' column name 'peaks'.
#' @param tol \code{double} value to be used for m/z tolerance when merging
#' peak lists.
#' @return \code{dataframe} with each row representing shared peaks from
#' multiple peak lists based on m/z value
#' @example 
#' 
#' peakList1 <- data.frame('peaks'=c(615.3456, 489.6651, 375.1968))
#' peakList2 <- data.frame('peaks'=c(615.3589, 453.3596, 357.9618))
#' 
#' listOfPeaklists <- list(peakList1, peakList2)
#' 
#' combinedPeakList <- combinePeakLists(listOfPeakLists, tol=0.2)
#' 
#' @export
combinePeakLists <- function(peakLists, tol) {
  # Merge peak lists based on first peak list.
  for (i in 2:length(peakLists)) {
    if (i == 2) {
      bigDf <- difference_full_join(peakLists[[i-1]], peakLists[[i]],
                                    by='peaks', max_dist=tol)
    } else {
      bigDf <- difference_full_join(big_df, peakLists[[i]], by='peaks',
                                    max_dist=tol)
    }
    colnames(bigDf)[1] <- 'peaks'
  }
  colnames(bigDf) <- as.character(seq(1, ncol(bigDf)))
  
  refDf <- bigDf
  
  finalDfs <- list()
  for (m in seq(1, ncol(bigDf)-2)) {
    # Remove first peak list and get peaks from second peak list that were not
    # merged to first list for merging with subsequent peak lists.
    tmpPeakLists <- list()
    for (i in 2:ncol(refDf)) {
      tmpDf <- as.data.frame(refDf[which(is.na(refDf[, 1])),][, i])
      colnames(tmpDf) <- c('peaks')
      tmpDf <- as.data.frame(tmpDf[which(!is.na(tmpDf$peaks)),])
      colnames(tmpDf) <- c('peaks')
      tmpPeakLists[[length(tmpPeakLists) + 1]] <- tmpDf
    }
    
    # Iteratively merge peak lists based on subsequent peak lists, similar to
    # first merge above.
    for (i in 2:length(tmpPeakLists)) {
      if (i == 2) {
        tmpDf <- difference_full_join(tmpPeakLists[[i-1]], tmpPeakLists[[i]],
                                      by='peaks', max_dist=tol)
      } else {
        tmpDf <- difference_full_join(tmpDf, tmpPeakLists[[i]], by='peaks',
                                      max_dist=tol)
      }
      colnames(tmpDf)[1] <- 'peaks'
    }
    colnames(tmpDf) <- as.character(seq(1 + m, ncol(bigDf)))
    
    refDf <- tmpDf
    
    # Add columns of NA to match full number of columns in full dataframe of
    # peak lists.
    for (j in 1:m) {
      tmpDf[, as.character(j)] <- rep(NA, nrow(tmp))
    }
    
    # Append peak lists to finalDfs list.
    if (m != ncol(bigDf) - 2) {
      finalDfs[[length(finalDfs) + 1]] <- tmpDf[which(!is.na(tmpDf[, 1])),]
    } else {
      finalDfs[[length(finalDfs) + 1]] <- tmpDf
    }
  }
  
  # Concatenate each datafram in list of finalDfs to make final dataframe of
  # combined peak lists.
  for (i in 1:length(finalDfs)) {
    if (i == 1) {
      finalDf <- rbind(bigDf[which(!is.na(bigDf[, 1])),], finalDfs)
    } else {
      finalDf <- rbind(finalDf, finalDfs[[i]])
    }
  }
  
  return(finalDf)
}
