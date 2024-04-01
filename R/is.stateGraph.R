#' Stategraph class
#'
#' @description
#'
#' Functions to deal with the stategraph class.
#'
#' @param x An object of class \code{stateGraph}.
#'
#' @details
#'
#' Claddis uses various classes to define specific types of data, here the use of a character stategraph (to specify the parsimony costs of transitions between character states) are assigned the class "stateGraph".
#'
#' \code{is.stateGraph} checks whether an object is or is not a valid stateGraph object.
#'
#' @return \code{is.stateGraph} returns either TRUE or FALSE.
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
#' # Convert costmatrix to stategraph:
#' stategraph <- convert_costmatrix_to_stategraph(costmatrix = costmatrix)
#'
#' # Check that this is a valid costMatrix object:
#' is.stateGraph(x = stategraph)
#'
#' @export is.stateGraph
is.stateGraph <- function(x) {
  
  # Get any error messages for stategraph:
  messages <- check_stateGraph(stategraph = x)
  
  # Return logical indicating whether object is a valid stategraph object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}
