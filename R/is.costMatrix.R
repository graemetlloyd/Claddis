#' Costmatrix class
#'
#' @description
#'
#' Functions to deal with the costmatrix class.
#'
#' @param x An object of class \code{costMatrix}.
#'
#' @details
#'
#' Claddis uses various classes to define specific types of data, here the use of costmatrices (to specify the parsimony costs of transitions between character states) are assigned the class "costMatrix".
#'
#' \code{is.costMatrix} checks whether an object is or is not a valid costMatrix object.
#'
#' @return \code{is.costMatrix} returns either TRUE or FALSE.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make an unordered costmatrix:
#' costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered"
#' )
#'
#' # Check that this is a valid costMatrix object:
#' is.costMatrix(x = costmatrix)
#'
#' @export is.costMatrix
is.costMatrix <- function(x) {
  
  # Get any error messages for costmatrix:
  messages <- check_costMatrix(costmatrix = x)
  
  # Return logical indicating whether object is a valid costmatrix object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}
