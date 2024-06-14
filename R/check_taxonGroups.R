#' Check taxonGroups object for errors
#'
#' @description
#'
#' Internal function to check taxonGroups object for errors.
#'
#' @param taxon_groups An object of class \code{taxonGroups}.
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
#' # Create a taxon groups object:
#' taxon_groups <- list(
#'   Group_A = c("Species_1", "Species_2", "Species_3"),
#'   Group_B = c("Species_3", "Species_4"),
#'   Group_C = c("Species_5", "Species_6", "Species_7", "Species_8")
#' )
#'
#' # Check that this is a valid taxonGroups object (will return error message as class
#' # is not set):
#' check_taxonGroups(taxon_groups = taxon_groups)
#'
#' @export check_taxonGroups
check_taxonGroups <- function(taxon_groups) {
  
  # Check taxon_groups has class taxonGroups and add error message to output if true:
  if (!inherits(x = taxon_groups, what = "taxonGroups")) return("taxon_groups must be an object of class \"taxonGroups\".")
  
  # Check taxon_groups are in form of list add error message to output if false:
  if (!is.list(x = taxon_groups)) return("taxon_groups must be in the form of a list.")
  
  # Check taxon_groups has at least one group:
  if (length(x = taxon_groups) == 0) return("taxon_groups must have at least one element.")
  
  # Check taxon_groups collectively contains at least one taxon:
  if (all(x = !unlist(x = lapply(X = taxon_groups, FUN = function(x) length(x = x))) > 0)) return("taxon_groups must have at least one element containing taxa.")
  
  # Check group names are set:
  if (is.null(x = names(x = taxon_groups))) return("taxon_groups must have names set for each group.")
  
  # Check group names are all unique:
  if (any(x = duplicated(x = names(x = taxon_groups)))) return("taxon_groups must have unique group names.")
  
  # Check group names are not empty strings:
  if (any(x = nchar(x = names(x = taxon_groups)) == 0)) return("taxon_groups must have group names of positive length.")
  
  # Check taxa are in form of vectors:
  if (any(x = !unlist(x = lapply(X = taxon_groups, FUN = function(x) is.vector(x = x))))) return("taxon_groups must be composed of vectors of taxon names.")
  
  # Check taxa are vaid character formats:
  if (any(!unlist(x = lapply(X = taxon_groups, FUN = function(x) is.character(x = x))))) return("taxon_groups must be composed of character vectors.")
  
  # Check no taxa are duplicated within groups:
  if (any(x = unlist(x = lapply(X = taxon_groups, FUN = function(x) duplicated(x = x))))) return("taxon_groups must not contain duplicated taxa within groups.")
  
  # Check no taxa are empty strings:
  if (any(unlist(x = lapply(X = taxon_groups, FUN = function(x) nchar(x = x))) == 0)) return("taxon_groups must contain taxon anmes of positive length.")
  
  # Return empty vector:
  vector(mode = "character")
}
