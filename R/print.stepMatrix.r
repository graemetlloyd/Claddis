#' Compact display of a step matrix
#'
#' @description
#'
#' Displays a compact summary of a stepMatrix object.
#'
#' @param x An object of class \code{"stepMatrix"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#'
#' Displays some basic summary information on a stepmatrix object.
#'
#' @return
#'
#' Nothing is directly returned, instead a text summary describing a \code{"stepMatrix"} object is printed to the console.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make an unordered stepmatrix:
#' example_stepmatrix <- make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Show print.stepMatrix version of each included data sets:
#' print.stepMatrix(x = example_stepmatrix)
#'
#' @export print.stepMatrix
print.stepMatrix <- function(x, ...) {
  
  # Check x has class stepMatrix and stop and warn user if not:
  if (!inherits(x = x, what = "stepMatrix")) stop("x must be an object of class \"stepMatrix\".")
  
  # If not a valid stepMatrix object then stop and provide feedback to user on what is wrong:
  #if (!is.stepMatrix(x = x)) stop(check_taxonGroups(taxon_groups = x)[1])
  
  
  # Return summary information about object:
  cat(paste0("taxonGroups object composed of ", length(x = x), " groups:"), "\n", unlist(x = lapply(X = as.list(x = names(x = x)), function(y) paste0(" ", y, paste0(rep(x = " ", times = max(x = nchar(x = names(x = x))) - nchar(x = y) + 1), collapse = ""), "(", length(x = x[[y]]), " taxa", format_taxon_names(x = x[[y]]), ")\n"))))
  
}
