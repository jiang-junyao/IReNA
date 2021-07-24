#' Extract 'chr', 'start', 'end' columns of peak file
#'
#' @param peak peak file which should contain 'chr', 'start', 'end' columns
#'
#' @return return bed format data frame
#' @export
#'
#' @examples
get_bed <- function(peak) {
  col1 <- peak[, c("chr", "start", "end")]
  return(col1)
}
