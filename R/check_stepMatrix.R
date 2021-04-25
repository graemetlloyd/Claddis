#' Check stepMatrix object for errors
#'
#' @description
#'
#' Internal function to check stepMatrix object for errors.
#'
#' @param stepmatrix A stepMatrix object.
#'
#' @details
#'
#' Internal Claddis function. Nothing to see here. Carry on.
#'
#' Solves issue raised in Maddison and Maddison (2003) by checking custom stepmatrices are internally consistent such that each transition cost represents the cost of the shortest possible path.
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
#' # Check that this is a valid stepMatrix object (should return empty vector):
#' check_stepMatrix(stepmatrix = stepmatrix)
#'
#' @export check_stepMatrix
check_stepMatrix <- function(stepmatrix) {
  
  # Check stepmatrix has class stepMatrix and add error message to output if true:
  if (!inherits(x = stepmatrix, what = "stepMatrix")) return("stepmatrix must be an object of class \"stepMatrix\".")
  
  # Check stepmatrix is in form of list an add error message to output if false:
  if (!is.list(x = stepmatrix)) return("stepmatrix must be in the form of a list.")
  
  # Check length of list is exactly five and add error message to output if false:
  if (length(x = stepmatrix) != 5) return("stepmatrix must be a list with five items (size, type, stepmatrix, symmetry, and include_polymorphisms).")
  
  # Check names are correct and in order add error message to output if false:
  if (!all(x = names(x = stepmatrix) == c("size", "type", "stepmatrix", "symmetry", "includes_polymorphisms"))) return("Elements of stepmatrix must be \"size\", \"type\", \"stepmatrix\", \"symmetry\", and \"include_polymorphisms\", in that order.")
  
  # Check size is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stepmatrix$size) || length(x = stepmatrix$size) != 1) stop("stepmatrix$size should be a single numeric value indicating the size (number of rows or columns) of the stepmatrix.")
  
  # Check type is formatted correctly and add error message to output if false:
  if (!is.character(x = stepmatrix$type) || length(x = stepmatrix$type) != 1) stop("stepmatrix$type should be a single character value indicating the type of stepmatrix.")
  
  # Check type is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = stepmatrix$type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy", "custom"))) > 0) stop("stepmatrix$type must be one of: \"ordered\", \"unordered\", \"dollo\", \"irreversible\", \"stratigraphy\", \"custom\".")
  
  # Check stepmatrix$stepmatrix is an actual matrix and add error message to output if false:
  if (!is.matrix(x = stepmatrix$stepmatrix)) stop("stepmatrix$stepmatrix must be a matrix.")
  
  # Check stepmatrix$stepmatrix size matches stepmatrix$size and add error message to output if false
  if (length(x = stepmatrix$stepmatrix) != stepmatrix$size ^ 2) stop("The size of stepmatrix$stepmatrix must match stepmatrix$size.")
  
  # Check stepmatrix$stepmatrix is only numeric values and add error message to output if false:
  if (!is.numeric(stepmatrix$stepmatrix)) stop("All elements of stepmatrix$stepmatrix must be numeric.")
  
  # Check stepmatrix$stepmatrix is only positive or zero values and add error message to output if false:
  if (!all(x = stepmatrix$stepmatrix >= 0)) stop("All elements of stepmatrix$stepmatrix must be either zero or positive.")
  
  # Check stepmatrix$stepmatrix diagonal is only zero values and add error message to output if false:
  if (!all(x = diag(x = stepmatrix$stepmatrix) == 0)) stop("All diagonal elements of stepmatrix$stepmatrix must be zero ")
  
  # Check stepmatrix$stepmatrix off-diagonal is only positive values and add error message to output if false:
  if (!all(x = c(stepmatrix$stepmatrix[lower.tri(x = stepmatrix$stepmatrix)], stepmatrix$stepmatrix[upper.tri(x = stepmatrix$stepmatrix)]) > 0)) stop("All non-diagonal elements of stepmatrix$stepmatrix must be positive ")
  
  # Check symmetry is formatted correctly and add error message to output if false:
  if (!is.character(x = stepmatrix$symmetry) || length(x = stepmatrix$symmetry) != 1) stop("stepmatrix$symmetry should be a single character value indicating whether the stepmatrix is symmetric or asymmetric.")
  
  # Check symmetry is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = stepmatrix$symmetry, y = c("Asymmetric", "Symmetric"))) > 0) stop("stepmatrix$type must be one of: \"Asymmetric\", \"Symmetric\".")
  
  # Check symmetry is correct and add error message to output if false:
  if (isSymmetric(object = stepmatrix$stepmatrix) != (stepmatrix$symmetry == "Symmetric")) stop("stepmatrix$stepmatrix and stepmatrix$symmetry must match (i.e., if symmetric must be Symmetric, if asymmetric then Asymmetric).")
  
  # Check includes_polymorphisms is formatted correctly and add error message to output if false:
  if (!is.logical(x = stepmatrix$includes_polymorphisms) || length(x = stepmatrix$includes_polymorphisms) != 1) stop("stepmatrix$includes_polymorphisms should be a single logical value indicating whether polymorphisms are included or not.")
  
  # Check matrix makes sense (but only if custom):
  if (stepmatrix$type == "custom") {
    
    # Subfunction to build shortest path stepmatrix:
    build_shortest_path_stepmatrix <- function(stepmatrix) {
      
      # Get stepmatrix states:
      states <- rownames(x = stepmatrix$stepmatrix)
      
      # Get length of each shortest path:
      path_lengths <- apply(
        X = expand.grid(start = states, end = states),
        MARGIN = 1,
        FUN = function(x) {
          path <- find_shortest_path(stepmatrix = stepmatrix, start = x[1], end = x[2])
          lengths <- lapply(
            X = as.list(x = 2:length(x = path)),
            function(i) stepmatrix$stepmatrix[path[(i - 1)], path[i]]
          )
          sum(x = unlist(x = lengths))
        }
      )
      
      # Return shortest path stepmatrix:
      matrix(data = path_lengths, nrow = stepmatrix$size, dimnames = list(states, states))
    }
    
    # Build stepmatrix where are costs are shortest path costs:
    shortest_path_stepmatrix <- build_shortest_path_stepmatrix(stepmatrix = stepmatrix)
    
    # If stepmatrix is not all shortest paths stop and warn user:
    if (!all(x = shortest_path_stepmatrix == stepmatrix$stepmatrix)) stop("stepmatrix is not internally consistent (at least one path is shorter - lower cost - than stated). Fix using find_shortest_path and try again.")
  }

  # Return empty vector:
  vector(mode = "character")
}
