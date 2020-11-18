#' Check cladisticMatrix object for errors
#'
#' @description
#'
#' Internal function to check cladisticMatrix object for errors.
#'
#' @param cladistic_matrix An object of class \code{cladisticMatrix}.
#'
#' @details
#'
#' Internal Claddis function. Nothing to see here. Carry on.
#'
#' @return An error message or empty vector if no errors found.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Check that this is a valid cladisticMatrix object (will return error message as class
#' # is not set):
#' check_cladisticMatrix(cladistic_matrix = day_2016)
#'
#' @export check_cladisticMatrix
check_cladisticMatrix <- function(cladistic_matrix) {
  
  # TO DO:
  #
  # Make this class hierarchical? E.g., add matrixBlock class. Also, maybe add just a taxon names part to output.
  
  # Check cladistic_matrix has class cladisticMatrix and add error message to output if true:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) return("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
  # Check cladistic_matrix are in form of list add error message to output if false:
  if (!is.list(x = cladistic_matrix)) return("cladistic_matrix must be in the form of a list.")
  
  # Check there are at least two elements to cladistic_matrix:
  if (length(x = cladistic_matrix) < 2) return("cladistic_matrix must have at least two elements.")

  # Check first element is called topper:
  if (names(cladistic_matrix)[1] != "topper") return("cladistic_matrix must have at first element called \"topper\".")

# WAY MORE STUFF HERE

  # Return empty vector:
  vector(mode = "character")
}
