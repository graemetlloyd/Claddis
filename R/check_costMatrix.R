#' Check costMatrix object for errors
#'
#' @description
#'
#' Internal function to check costMatrix object for errors.
#'
#' @param costmatrix A costMatrix object.
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
#' # Make an unordered costmatrix:
#' costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Check that this is a valid costMatrix object (should return empty vector):
#' check_costMatrix(costmatrix = costmatrix)
#'
#' @export check_costMatrix
check_costMatrix <- function(costmatrix) {
  
  # Check costmatrix has class costMatrix and add error message to output if true:
  if (!inherits(x = costmatrix, what = "costMatrix")) return("costmatrix must be an object of class \"costMatrix\".")
  
  # Check costmatrix is in form of list and add error message to output if false:
  if (!is.list(x = costmatrix)) return("costmatrix must be in the form of a list.")
  
  # Check length of list is exactly five and add error message to output if false:
  if (length(x = costmatrix) != 5) return("costmatrix must be a list with five items (size, type, costmatrix, symmetry, and include_polymorphisms).")
  
  # Check names are correct and in order add error message to output if false:
  if (!all(x = names(x = costmatrix) == c("size", "type", "costmatrix", "symmetry", "includes_polymorphisms"))) return("Elements of costmatrix must be \"size\", \"type\", \"costmatrix\", \"symmetry\", and \"include_polymorphisms\", in that order.")
  
  # Check size is formatted correctly and add error message to output if false:
  if (!is.numeric(x = costmatrix$size) || length(x = costmatrix$size) != 1) stop("costmatrix$size should be a single numeric value indicating the size (number of rows or columns) of the costmatrix.")
  
  # Check type is formatted correctly and add error message to output if false:
  if (!is.character(x = costmatrix$type) || length(x = costmatrix$type) != 1) stop("costmatrix$type should be a single character value indicating the type of costmatrix.")
  
  # Check type is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = costmatrix$type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy", "custom"))) > 0) stop("costmatrix$type must be one of: \"ordered\", \"unordered\", \"dollo\", \"irreversible\", \"stratigraphy\", \"custom\".")
  
  # Check costmatrix$costmatrix is an actual matrix and add error message to output if false:
  if (!is.matrix(x = costmatrix$costmatrix)) stop("costmatrix$costmatrix must be a matrix.")
  
  # Check costmatrix$costmatrix size matches costMatrix$size and add error message to output if false
  if (length(x = costmatrix$costmatrix) != costmatrix$size ^ 2) stop("The size of costmatrix$costmatrix must match costmatrix$size.")
  
  # Check costmatrix$costmatrix is only numeric values and add error message to output if false:
  if (!is.numeric(costmatrix$costmatrix)) stop("All elements of costmatrix$costmatrix must be numeric.")
  
  # Check costmatrix$costmatrix is only positive or zero values and add error message to output if false:
  if (!all(x = costmatrix$costmatrix >= 0)) stop("All elements of costmatrix$costmatrix must be either zero or positive.")
  
  # Check costmatrix$costmatrix diagonal is only zero values and add error message to output if false:
  if (!all(x = diag(x = costmatrix$costmatrix) == 0)) stop("All diagonal elements of costmatrix$costmatrix must be zero.")
  
  # Check costmatrix$costmatrix off-diagonal is only positive values and add error message to output if false:
  if (!all(x = c(costmatrix$costmatrix[lower.tri(x = costmatrix$costmatrix)], costmatrix$costmatrix[upper.tri(x = costmatrix$costmatrix)]) > 0)) stop("All non-diagonal elements of costmatrix$costmatrix must be positive ")
  
  # Check symmetry is formatted correctly and add error message to output if false:
  if (!is.character(x = costmatrix$symmetry) || length(x = costmatrix$symmetry) != 1) stop("costmatrix$symmetry should be a single character value indicating whether the costmatrix is symmetric or asymmetric.")
  
  # Check symmetry is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = costmatrix$symmetry, y = c("Asymmetric", "Symmetric"))) > 0) stop("costmatrix$symmetry must be one of: \"Asymmetric\", \"Symmetric\".")
  
  # Check symmetry is correct and add error message to output if false:
  if (isSymmetric(object = costmatrix$costmatrix) != (costmatrix$symmetry == "Symmetric")) stop("costmatrix$costmatrix and costmatrix$symmetry must match (i.e., if symmetric must be Symmetric, if asymmetric then Asymmetric).")
  
  # Check includes_polymorphisms is formatted correctly and add error message to output if false:
  if (!is.logical(x = costmatrix$includes_polymorphisms) || length(x = costmatrix$includes_polymorphisms) != 1) stop("costmatrix$includes_polymorphisms should be a single logical value indicating whether polymorphisms are included or not.")
  
  # Check matrix makes sense (but only if custom):
  if (costmatrix$type == "custom") {
    
    # Build costmatrix where costs are shortest path costs:
    shortest_path_costmatrix <- fix_costmatrix(costmatrix = costmatrix, message = FALSE)
    
    # If costMatrix is not all shortest paths stop and warn user:
    if (!all(x = shortest_path_costmatrix$costmatrix == costmatrix$costmatrix)) stop("costmatrix is not self-consistent (at least one path is shorter - lower cost - than stated). Fix using fix_costmatrix and try again.")
  }

  # Return empty vector:
  vector(mode = "character")
}
