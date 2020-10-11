#' Compact display of taxon groups
#'
#' @description
#'
#' Displays a compact summary of a taxonGroups object.
#'
#' @param x An object of class \code{"taxonGroups"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#'
#' Displays some basic summary information on a taxon groups object, including number of groups and their names and partial contents.
#'
#' @return
#'
#' Nothing is directly returned, instead a text summary describing a \code{"taxonGroups"} object is printed to the console.
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
#' # Set class as taxonGroups:
#' class(taxon_groups) <- "taxonGroups"
#'
#' # Show print.taxonGroups version of each included data sets:
#' print.taxonGroups(x = taxon_groups)
#' @export print.taxonGroups
print.taxonGroups <- function(x, ...) {
  
  # Check x has class taxonGroups and stop and warn user if not:
  if (!inherits(x = x, what = "taxonGroups")) stop("x must be an object of class \"taxonGroups\".")
  
  # If not a valid taxonGroups object then stop and provide feedback to user on what is wrong:
  if (!is.taxonGroups(x = x)) stop(check_taxonGroups(taxon_groups = x)[1])
  
  # Sub function to formt taxon names output:
  format_taxon_names <- function(x) {
    if (length(x) == 0) return("")
    if (length(x) == 1) return(paste0(": ", x, collapse = ""))
    if (length(x) == 2) return(paste0(": ", x[1], ", ", x[2], collapse = ""))
    if (length(x) == 3) return(paste0(": ", x[1], ", ", x[2], ", ", x[3], collapse = ""))
    if (length(x) > 3) return(paste0(": ", x[1], ", ", x[2], ", ", x[3], ", ...", collapse = ""))
  }
  
  # Return summary information about object:
  cat(paste0("taxonGroups object composed of ", length(x = x), " groups:"), "\n", unlist(x = lapply(X = as.list(x = names(x = x)), function(y) paste0(" ", y, paste0(rep(x = " ", times = max(x = nchar(x = names(x = x))) - nchar(x = y) + 1), collapse = ""), "(", length(x = x[[y]]), " taxa", format_taxon_names(x = x[[y]]), ")\n"))))
  
}
