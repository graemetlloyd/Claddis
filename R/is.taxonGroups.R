#' Taxon groups class
#'
#' @description
#'
#' Functions to deal with the taxon groups class.
#'
#' @param x An object of class \code{taxonGroups}.
#'
#' @details
#'
#' Claddis uses various classes to define specific types of data, here the use of taxon groups (to delineate different groups of taxa, e.g., clades, time bins, geographic regions etc.) ae assigned the class "taxonGroups".
#'
#' \code{is.taxonGroups} checks whether an object is or is not a valid taxonGroups object.
#'
#' @return \code{is.taxonGroups} returns either TRUE or FALSE.
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
#' # Check that this is a valid taxonGroups object (will fail as class is not set):
#' is.taxonGroups(x = taxon_groups)
#'
#' # Set class as taxonGroups:
#' class(taxon_groups) <- "taxonGroups"
#'
#' # Check that this is a valid taxonGroups object (will succeed as format and
#' # class are correct):
#' is.taxonGroups(x = taxon_groups)
#'
#' @export is.taxonGroups
is.taxonGroups <- function(x) {
  
  # Get any error messages for taxon_groups:
  messages <- check_taxonGroups(taxon_groups = x)
  
  # Return logical indicating whether object is a valid taxonGroups object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}
