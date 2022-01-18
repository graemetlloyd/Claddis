#' Compact display of a costmatrix
#'
#' @description
#'
#' Displays a compact summary of a costMatrix object.
#'
#' @param x An object of class \code{"costMatrix"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#'
#' Displays some basic summary information on a costmatrix object.
#'
#' @return
#'
#' Nothing is directly returned, instead a text summary describing a \code{"costMatrix"} object is printed to the console.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make an unordered costmatrix:
#' example_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Show print.costMatrix version:
#' print.costMatrix(x = example_costmatrix)
#'
#' @export print.costMatrix
print.costMatrix <- function(x, ...) {
  
  # ANOTHER USEFUL THING TO STATE IS WHETHER MATRIX CAN BE REPRESENTED AS AN ADJACENCY MATRIX?
  # NEED FUNCTION TO CONVERT COSTMATRIX TO Q-MATRIX PARAMETERS? NOT ALL DOABLE, BUT SEEMS USEFUL IF APPLYING LIKELIHOOD ELSEWHERE IN CLADDIS
  
  # Check x has class costMatrix and stop and warn user if not:
  if (!inherits(x = x, what = "costMatrix")) stop("x must be an object of class \"costMatrix\".")
  
  # If not a valid costMatrix object then stop and provide feedback to user on what is wrong:
  if (!is.costMatrix(x = x)) stop(check_costMatrix(costmatrix = x)[1])
  
  # Return summary information about object:
  cat(paste0(x$symmetry, " ", x$type, " costMatrix object containing ", x$size, " unique states", ifelse(test = x$includes_polymorphisms, yes = " (including polymorphic states)", no = ""), "."))
  
}