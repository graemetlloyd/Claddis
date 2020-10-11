#' Cladistic matrix class
#'
#' @description
#'
#' Functions to deal with the cladistic matrix class.
#'
#' @param x A cladisticMatrix object.
#'
#' @details
#'
#' Claddis uses various classes to define specific types of data, here the use of cladistic data (the main input of the package) is assigned the class "cladisticMatrix".
#'
#' \code{is.cladisticMatrix} checks whether an object is or is not a valid cladisticMatrix object.
#'
#' @return \code{is.cladisticMatrix} returns either TRUE or FALSE.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Check that this is a valid cladisticMatrix object (will succeed as format and
#' # class are correct):
#' is.cladisticMatrix(x = day_2016)
#'
#' @export is.cladisticMatrix
is.cladisticMatrix <- function(x) {
  
  # Get any error messages for cladistic_matrix:
  messages <- check_cladisticMatrix(cladistic_matrix = x)
  
  # Return logical indicating whether object is a valid cladisticMatrix object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}
