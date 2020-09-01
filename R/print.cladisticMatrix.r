#' Compact display of a cladistic matrix
#'
#' @description
#'
#' Displays a compact summary of the dimensions and nature of a cladistic matrix object.
#'
#' @param cladistic_matrix An object of class \code{"cladisticMatrix"}.
#'
#' @details
#'
#' TO ADD.
#'
#' @return
#'
#' A summary describing the dimensions and nature of an object of class \code{"cladisticMatrix"}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' TO ADD
#'
#' @examples
#'
#' # Show print.cladisticMatrix version of each included data sets:
#' print.cladisticMatrix(cladistic_matrix = day_2016)
#' print.cladisticMatrix(cladistic_matrix = gauthier_1986)
#' print.cladisticMatrix(cladistic_matrix = michaux_1989)
#'
#' @export print.cladisticMatrix
#' @exportClass cladisticMatrix
print.cladisticMatrix <- function(cladistic_matrix) {
  
  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
  n_blocks <- length(cladistic_matrix) - 1
  
  block_sizes <- unname(obj = unlist(x = lapply(X = cladistic_matrix[2:length(cladistic_matrix)], FUN = function(x) ncol(x = x$matrix))))
  
  cat(paste0("Object of type \"cladisticMatrix\", containing ", n_blocks, " matrix blocks.", "\n"))
  
  cat(paste(block_sizes, collapse = " and "))
  
}
