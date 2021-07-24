#' inner function to split dataframe
#'

#' @return
#'
#' @examples
split_dataframe <- function(dataframe1, sep = '\t'){
  out1 <- strsplit(dataframe1[1,],sep)[[1]]
  out1 <- matrix(out1,ncol = length(out1))
  for (i in 2:nrow(dataframe1)) {
    out2 <- strsplit(dataframe1[i,],sep)[[1]]
    out1 <- rbind(out1,out2)
  }
  out1 <- as.data.frame(out1)
  return(out1)
}
