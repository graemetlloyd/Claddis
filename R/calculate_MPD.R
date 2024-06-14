#' Calculate mean pairwise distances
#'
#' @description
#'
#' Given distanceMatrices and taxonGroups objects calculates their mean pairwise distances.
#'
#' @param distances An object of class \code{distanceMatrices}.
#' @param taxon_groups An object of class \code{taxonGroups}.
#'
#' @details
#'
#' Not all measures of disparity (morphological distance) require an ordination space. For example, the pariwise distances between taxa are themselves a disparity metric. This function takes the output from \link{calculate_morphological_distances} and a set of taxon groups and returns the mean pairwise distance for each of those groups.
#'
#' @return
#'
#' A labelled vector of weighted mean pairwise distances.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get morphological distances for the Day et al. (2016) data set:
#' distances <- calculate_morphological_distances(
#'   cladistic_matrix = day_2016,
#'   distance_metric = "mord",
#'   distance_transformation = "none"
#' )
#'
#' # Build simple taxonomic groups for Day et al. (2016) data set:
#' taxon_groups <- list(nonBurnetiamorpha = c("Biarmosuchus_tener", "Hipposaurus_boonstrai",
#'   "Bullacephalus_jacksoni", "Pachydectes_elsi", "Niuksenitia_sukhonensis", "Ictidorhinus_martinsi",
#'   "RC_20", "Herpetoskylax_hopsoni"), Burnetiamorpha = c("Lemurosaurus_pricei", "Lobalopex_mordax",
#'   "Lophorhinus_willodenensis", "Proburnetia_viatkensis", "Lende_chiweta",
#'   "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098"))
#'
#' # Set class as taxonGroups:
#' class(taxon_groups) <- "taxonGroups"
#'
#' # Calculate mean pairwise distances:
#' calculate_MPD(distances, taxon_groups)
#'
#' @export calculate_MPD
calculate_MPD <- function(distances, taxon_groups) {
  
  # If not a valid taxonGroups object then stop and provide feedback to user on what is wrong:
  if (!is.taxonGroups(x = taxon_groups)) stop(check_taxonGroups(taxon_groups = taxon_groups)[1])
  
  # Calculate and return mean pairwise distance for each taxon group:
  unlist(
    x = lapply(
      X = taxon_groups,
      FUN = function(x) {
        mean(
          x = distances$distance_matrix[x, x][lower.tri(x = distances$distance_matrix[x, x])],
          na.rm = TRUE
        )
        
      }
    )
  )
}
