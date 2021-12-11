#' Check costMatrix object for errors
#'
#' @description
#'
#' Internal function to check costMatrix object for errors.
#'
#' @param costMatrix A costMatrix object.
#'
#' @details
#'
#' Internal Claddis function. Nothing to see here. Carry on.
#'
#' Solves issue raised in Maddison and Maddison (2003) by checking custom costmatrices are internally consistent such that each transition cost represents the cost of the shortest possible path.
#'
#' @return An error message or empty vector if no errors found.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Maddison, D. R. and Maddison, W. P., 2003. \emph{MacClade 4: Analysis of phylogeny and character evolution}. Version 4.06. Sinauer Associates, Sunderland, Massachusetts.
#'
#' @examples
#'
#' # Make an unordered costMatrix:
#' costMatrix <- make_costMatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Check that this is a valid costMatrix object (should return empty vector):
#' check_costMatrix(costMatrix = costMatrix)
#'
#' @export check_costMatrix
check_costMatrix <- function(costMatrix) {
  
  # Check costMatrix has class costMatrix and add error message to output if true:
  if (!inherits(x = costMatrix, what = "costMatrix")) return("costMatrix must be an object of class \"costMatrix\".")
  
  # Check costMatrix is in form of list and add error message to output if false:
  if (!is.list(x = costMatrix)) return("costMatrix must be in the form of a list.")
  
  # Check length of list is exactly five and add error message to output if false:
  if (length(x = costMatrix) != 5) return("costMatrix must be a list with five items (size, type, costMatrix, symmetry, and include_polymorphisms).")
  
  # Check names are correct and in order add error message to output if false:
  if (!all(x = names(x = costMatrix) == c("size", "type", "costMatrix", "symmetry", "includes_polymorphisms"))) return("Elements of costMatrix must be \"size\", \"type\", \"costMatrix\", \"symmetry\", and \"include_polymorphisms\", in that order.")
  
  # Check size is formatted correctly and add error message to output if false:
  if (!is.numeric(x = costMatrix$size) || length(x = costMatrix$size) != 1) stop("costMatrix$size should be a single numeric value indicating the size (number of rows or columns) of the costMatrix.")
  
  # Check type is formatted correctly and add error message to output if false:
  if (!is.character(x = costMatrix$type) || length(x = costMatrix$type) != 1) stop("costMatrix$type should be a single character value indicating the type of costMatrix.")
  
  # Check type is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = costMatrix$type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy", "custom"))) > 0) stop("costMatrix$type must be one of: \"ordered\", \"unordered\", \"dollo\", \"irreversible\", \"stratigraphy\", \"custom\".")
  
  # Check costMatrix$costMatrix is an actual matrix and add error message to output if false:
  if (!is.matrix(x = costMatrix$costMatrix)) stop("costMatrix$costMatrix must be a matrix.")
  
  # Check costMatrix$costMatrix size matches costMatrix$size and add error message to output if false
  if (length(x = costMatrix$costMatrix) != costMatrix$size ^ 2) stop("The size of costMatrix$costMatrix must match costMatrix$size.")
  
  # Check costMatrix$costMatrix is only numeric values and add error message to output if false:
  if (!is.numeric(costMatrix$costMatrix)) stop("All elements of costMatrix$costMatrix must be numeric.")
  
  # Check costMatrix$costMatrix is only positive or zero values and add error message to output if false:
  if (!all(x = costMatrix$costMatrix >= 0)) stop("All elements of costMatrix$costMatrix must be either zero or positive.")
  
  # Check costMatrix$costMatrix diagonal is only zero values and add error message to output if false:
  if (!all(x = diag(x = costMatrix$costMatrix) == 0)) stop("All diagonal elements of costMatrix$costMatrix must be zero ")
  
  # Check costMatrix$costMatrix off-diagonal is only positive values and add error message to output if false:
  if (!all(x = c(costMatrix$costMatrix[lower.tri(x = costMatrix$costMatrix)], costMatrix$costMatrix[upper.tri(x = costMatrix$costMatrix)]) > 0)) stop("All non-diagonal elements of costMatrix$costMatrix must be positive ")
  
  # Check symmetry is formatted correctly and add error message to output if false:
  if (!is.character(x = costMatrix$symmetry) || length(x = costMatrix$symmetry) != 1) stop("costMatrix$symmetry should be a single character value indicating whether the costMatrix is symmetric or asymmetric.")
  
  # Check symmetry is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = costMatrix$symmetry, y = c("Asymmetric", "Symmetric"))) > 0) stop("costMatrix$type must be one of: \"Asymmetric\", \"Symmetric\".")
  
  # Check symmetry is correct and add error message to output if false:
  if (isSymmetric(object = costMatrix$costMatrix) != (costMatrix$symmetry == "Symmetric")) stop("costMatrix$costMatrix and costMatrix$symmetry must match (i.e., if symmetric must be Symmetric, if asymmetric then Asymmetric).")
  
  # Check includes_polymorphisms is formatted correctly and add error message to output if false:
  if (!is.logical(x = costMatrix$includes_polymorphisms) || length(x = costMatrix$includes_polymorphisms) != 1) stop("costMatrix$includes_polymorphisms should be a single logical value indicating whether polymorphisms are included or not.")
  
  # Check matrix makes sense (but only if custom):
  if (costMatrix$type == "custom") {
    
    # Build costMatrix where are costs are shortest path costs:
    shortest_path_costMatrix <- fix_costMatrix(costMatrix = costMatrix, message = FALSE)
    
    # If costMatrix is not all shortest paths stop and warn user:
    if (!all(x = shortest_path_costMatrix$costMatrix == costMatrix$costMatrix)) stop("costMatrix is not self-consistent (at least one path is shorter - lower cost - than stated). Fix using fix_costMatrix and try again.")
  }

  # Return empty vector:
  vector(mode = "character")
}
