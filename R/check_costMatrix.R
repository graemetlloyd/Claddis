#' Check a costMatrix object for errors
#'
#' @description
#'
#' Internal function to check a costMatrix object for errors.
#'
#' @param costmatrix A costMatrix object.
#'
#' @details
#'
#' Costmatrix objects are more complex than what will typically be shown to the user. This function checks this hidden structure and reports any errors it finds.
#'
#' These checks include rules 1-7 from Hoyal Cuthill and lloyd (i prep.).
#'
#' @return An error message or empty vector if no errors found.
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
#' # Check that this is a valid costMatrix object (should return empty vector):
#' check_costMatrix(costmatrix = costmatrix)
#'
#' @export check_costMatrix
check_costMatrix <- function(costmatrix) {
  
  # Check costmatrix has class costMatrix and add error message to output if true:
  if (!inherits(x = costmatrix, what = "costMatrix")) return("costmatrix must be an object of class \"costMatrix\".")
  
  # Check costmatrix is in form of list and add error message to output if false:
  if (!is.list(x = costmatrix)) return("costmatrix must be in the form of a list.")
  
  # Check length of list is at least fifteen and add error message to output if false:
  if (length(x = costmatrix) < 15) return("costmatrix must be a list with fifteen items (size, n_states, single_states, type, costmatrix, symmetry, includes_polymorphisms, polymorphism_costs, polymorphism_geometry, polymorphism_distance, includes_uncertainties, pruned, dollo_penalty, base_age, and weight).")
  
  # Check names are correct and in order add error message to output if false:
  if (!all(x = names(x = costmatrix) == c("size", "n_states", "single_states", "type", "costmatrix", "symmetry", "includes_polymorphisms", "polymorphism_costs", "polymorphism_geometry", "polymorphism_distance", "includes_uncertainties", "pruned", "dollo_penalty", "base_age", "weight"))) return("Elements of costmatrix must be \"size\", \"n_states\", \"single_states\", \"type\", \"costmatrix\", \"symmetry\", \"include_polymorphisms\", \"polymorphism_costs\", \"polymorphism_geometry\", \"polymorphism_distance\", \"includes_uncertainties\", \"pruned\", \"dollo_penalty\", \"base_age\", and \"weight\" in that order.")
  
  # Check size is formatted correctly and add error message to output if false:
  if (!is.numeric(x = costmatrix$size) || length(x = costmatrix$size) != 1) return("costmatrix$size should be a single numeric value indicating the size (number of rows or columns) of the costmatrix.")
  
  # Check n_states is formatted correctly and add error message to output if false:
  if (!is.numeric(x = costmatrix$n_states) || length(x = costmatrix$size) != 1) return("costmatrix$n_states should be a single numeric value indicating the number of singles states (i.e., excluding polymorphic or ucertain values).")

  # Check single_states are formatted correctly and match other elements and add error message to output if false:
  if (!is.character(x = costmatrix$single_states) || !is.vector(x = costmatrix$single_states) || length(x = costmatrix$single_states) != costmatrix$n_states || any(x = is.na(x = match(x = costmatrix$single_states, table = rownames(x = costmatrix$costmatrix)))) || length(x = grep(pattern = "/|&", x = costmatrix$single_states)) > 0) return("costmatrix$single_states must be a single character vector value equal in length to costmatrix$n_states, contain no polymorphic or uncertain values and match the states used in costmatrix$costmatrix.")
  
  # Check type is formatted correctly and add error message to output if false:
  if (!is.character(x = costmatrix$type) || length(x = costmatrix$type) != 1) return("costmatrix$type should be a single character value indicating the type of costmatrix.")
  
  # Check type is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = costmatrix$type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy", "custom"))) > 0) return("costmatrix$type must be one of: \"ordered\", \"unordered\", \"dollo\", \"irreversible\", \"stratigraphy\", \"custom\".")
  
  # Check costmatrix$costmatrix is an actual matrix and add error message to output if false:
  if (!is.matrix(x = costmatrix$costmatrix)) return("costmatrix$costmatrix must be a matrix.")
  
  # Check costmatrix$costmatrix size matches costMatrix$size and add error message to output if false
  if (length(x = costmatrix$costmatrix) != costmatrix$size ^ 2) return("The size of costmatrix$costmatrix must match costmatrix$size.")
  
  # RULE ONE: Check costmatrix$costmatrix is only numeric values and add error message to output if false:
  if (!is.numeric(costmatrix$costmatrix)) return("All elements of costmatrix$costmatrix must be numeric.")
  
  # Check costmatrix$costmatrix is only positive or zero values and add error message to output if false:
  if (!all(x = costmatrix$costmatrix >= 0)) return("All elements of costmatrix$costmatrix must be either zero or positive.")
  
  # RULE TWO: Check costmatrix$costmatrix diagonal is only zero values and add error message to output if false:
  if (!all(x = diag(x = costmatrix$costmatrix) == 0)) return("All diagonal elements of costmatrix$costmatrix must be zero.")
  
  # RULE THREE: Check costmatrix$costmatrix off-diagonal is only positive values and add error message to output if false:
  if (!all(x = c(costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states][lower.tri(x = costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states])], costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states][upper.tri(x = costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states])]) > 0)) return("All non-diagonal elements of costmatrix$costmatrix must be positive (except for uncertainties, if included).")
  
  # RULE FOUR: Check no state has infinite costs for every row and column (i.e., no state is disconnected from all others) and add error message if so:
  infinity_states <- do.call(
    what = rbind,
    args = lapply(
      X = as.list(costmatrix$single_states),
      FUN = function(x) {
        row_x <- costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states][x, ]
        column_x <- costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states][, x]
        row_x <- row_x[row_x > 0]
        column_x <- column_x[column_x > 0]
        is_infinity <- c(row_x == Inf, column_x == Inf)
      }
    )
  )
  rownames(x = infinity_states) <- costmatrix$single_states
  infinite_states <- names(x = costmatrix$single_states[apply(X = infinity_states, MARGIN = 1, FUN = all)])
  if (length(x = infinite_states) > 0) return("Costmatrix contains disconnected states (row ad column all infinite costs).")
  
  # RULE FIVE: Check no more than one state has all infinite "to" transition costs and add error message if not:
  infinite_columns <- costmatrix$single_states[apply(
    X = costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
    MARGIN = 2,
    FUN = function(x) {
      x <- x[x > 0] # Remove diagonal zero
      all(x == Inf)
    }
  )]
  if (length(infinite_columns) > 1) return("Multiple columns of costmatrix are all infinite cost meaning a tree must have multiple roots (impossible).")

  # RULE SIX: Check at least one state has all non-infinite "from" transition costs and add error message if not:
  infinite_rows <- costmatrix$single_states[apply(
    X = costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
    MARGIN = 1,
    FUN = function(x) {
      x <- x[x > 0] # Remove diagonal zero
      all(x != Inf)
    }
  )]
  if (length(infinite_rows) < 1) return("No row of costmatrix is free of infinite costs meaning no state can be a root state (impossible).")

  # Check symmetry is formatted correctly and add error message to output if false:
  if (!is.character(x = costmatrix$symmetry) || length(x = costmatrix$symmetry) != 1) return("costmatrix$symmetry should be a single character value indicating whether the costmatrix is symmetric or asymmetric.")
  
  # Check symmetry is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = costmatrix$symmetry, y = c("Asymmetric", "Symmetric"))) > 0) return("costmatrix$symmetry must be one of: \"Asymmetric\", \"Symmetric\".")
  
  # Check symmetry is correct and add error message to output if false:
  if (isSymmetric(object = costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states]) != (costmatrix$symmetry == "Symmetric")) return("costmatrix$costmatrix and costmatrix$symmetry must match (i.e., if symmetric must be Symmetric, if asymmetric then Asymmetric).")
  
  # Check includes_polymorphisms is formatted correctly and add error message to output if false:
  if (!is.logical(x = costmatrix$includes_polymorphisms) || length(x = costmatrix$includes_polymorphisms) != 1) return("costmatrix$includes_polymorphisms should be a single logical value indicating whether polymorphisms are included or not.")
  
  # Check that if includes_polymorphisms is TRUE there actually are polymorphisms and add error message to output if not:
  if (costmatrix$includes_polymorphisms && length(x = grep(pattern = "&", x = colnames(x = costmatrix$costmatrix))) == 0) return("costmatrix$includes_polymorphisms is set to TRUE, but there are no actual polymorphisms in costmatrix$costmatrix.")
  
  # Check that if includes_polymorphisms is FALSE there actually are no polymorphisms and add error message to output if not:
  if (!costmatrix$includes_polymorphisms && length(x = grep(pattern = "&", x = colnames(x = costmatrix$costmatrix))) > 0) return("costmatrix$includes_polymorphisms is set to FALSE, but there are actual polymorphisms in costmatrix$costmatrix.")
  
  # Check polymorphism_costs is a valid value and stop and add error message to output if not:
  if (length(x = costmatrix$polymorphism_costs) != 1 || !is.character(x = costmatrix$polymorphism_costs) || length(x = setdiff(x = costmatrix$polymorphism_costs, y = c("additive", "geometric", "maddison", "stratigraphic"))) > 0) return("costmatrix$polymorphism_costs must be a single character value and one of \"additive\", \"geometric\", \"maddison\", \"stratigraphic\".")
  
  # Check polymorphism_geometry is a valid value and stop and add error message to output if not:
  if (length(x = costmatrix$polymorphism_geometry) != 1 || !is.character(x = costmatrix$polymorphism_geometry) || length(x = setdiff(x = costmatrix$polymorphism_geometry, y = c("hypercube", "hypersphere", "simplex"))) > 0) stop("costmatrix$polymorphism_geometry must be a single character value and one of \"hypercube\", \"hypersphere\", \"simplex\".")
  
  # Check polymorphism_distance is a valid value and stop and add error message to output if not:
  if (length(x = costmatrix$polymorphism_distance) != 1 || !is.character(x = costmatrix$polymorphism_distance) || length(x = setdiff(x = costmatrix$polymorphism_distance, y = c("manhattan", "euclidean", "great_circle"))) > 0) stop("costmatrix$polymorphism_distance must be a single character value and one of \"manhattan\", \"euclidean\", \"great_circle\".")
  
  # Check includes_polymorphisms is formatted correctly and add error message to output if false:
  if (!is.logical(x = costmatrix$includes_uncertainties) || length(x = costmatrix$includes_uncertainties) != 1) return("costmatrix$includes_uncertainties should be a single logical value indicating whether uncertainties are included or not.")
  
  # Check that if includes_uncertainties is TRUE there actually are uncertainties and add error message to output if not:
  if (costmatrix$includes_uncertainties && length(x = grep(pattern = "/", x = colnames(x = costmatrix$costmatrix))) == 0) return("costmatrix$includes_uncertainties is set to TRUE, but there are no actual uncertainties in costmatrix$costmatrix.")
  
  # Check that if includes_uncertainties is FALSE there actually are no uncertainties and add error message to output if not:
  if (!costmatrix$includes_uncertainties && length(x = grep(pattern = "/", x = colnames(x = costmatrix$costmatrix))) > 0) return("costmatrix$includes_uncertainties is set to FALSE, but there are actual uncertainties in costmatrix$costmatrix.")

  # Check pruned is formatted correctly and add error message to output if false:
  if (!is.logical(x = costmatrix$pruned) || length(x = costmatrix$pruned) != 1) return("costmatrix$pruned should be a single logical value indicating whether the matrix is pruned (not all states are included) or not.")
  
  # Check dollo_penalty is formatted correctly and add error message to output if false:
  if (!is.numeric(x = costmatrix$dollo_penalty) || length(x = costmatrix$dollo_penalty) != 1 || costmatrix$dollo_penalty <= 0 || costmatrix$dollo_penalty == Inf) return("costmatrix$dollo_penalty should be a single finite positive numeric value that indicates the penalty to be used for a Dollo character.")
  
  # Check base_age is formatted correctly and add error message to output if false:
  if (!is.numeric(x = costmatrix$base_age) || length(x = costmatrix$base_age) != 1 || costmatrix$base_age <= 0 || costmatrix$base_age == Inf) return("costmatrix$base_age should be a single finite positive numeric value that indicates the oldest age for a stratigrphic character.")

  # Check weight is formatted correctly and add error message to output if false:
  if (!is.numeric(x = costmatrix$weight) || length(x = costmatrix$weight) != 1 || costmatrix$weight < 0 || costmatrix$weight == Inf) return("costmatrix$weight should be a single finite non-negative numeric value that indictaes the weight the character has.")
  
  # Check matrix makes sense (but only if custom):
  if (costmatrix$type == "custom") {
    
    # Build costmatrix where costs are shortest path costs:
    shortest_path_costmatrix <- fix_costmatrix(costmatrix = costmatrix, message = FALSE)
    
    # RULE SEVEN: If costMatrix is not all shortest paths add error message to output and warn user:
    if (!all(x = shortest_path_costmatrix$costmatrix == costmatrix$costmatrix)) return("costmatrix is not self-consistent (at least one path is shorter - lower cost - than stated). Fix using fix_costmatrix and try again.")
  }

  # Return empty vector:
  return(vector(mode = "character"))
}
