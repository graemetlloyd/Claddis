#' Stepmatrix class
#'
#' @description
#'
#' Functions to deal with the stepmatrix class.
#'
#' @param x An object of class \code{stepMatrix}.
#'
#' @details
#'
#' Claddis uses various classes to define specific types of data, here the use of stepmatrices (to specify the parsimony costs of transitions between character states) are assigned the class "stepMatrix".
#'
#' \code{is.stepMatrix} checks whether an object is or is not a valid stepMatrix object.
#'
#' @return \code{is.stepMatrix} returns either TRUE or FALSE.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make an unordered stepmatrix:
#' stepmatrix <- make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Check that this is a valid stepMatrix object:
#' is.stepMatrix(x = stepmatrix)
#'
#' @export is.stepMatrix
is.stepMatrix <- function(x) {
  
  # Get any error messages for stepmatrix:
  messages <- check_stepMatrix(stepmatrix = x)
  
  # Return logical indicating whether object is a valid stepmatrix object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}
