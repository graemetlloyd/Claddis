#' Compact display of a cladistic matrix
#'
#' @description
#'
#' Displays a compact summary of the dimensions and nature of a cladistic matrix object.
#'
#' @param x An object of class \code{"cladisticMatrix"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#'
#' Displays some basic summary information on a cladistic matrix object, including number and type of characters, information about ordering, and whether variable weights are used.
#'
#' @return
#'
#' Nothing is directly returned, instead a text summary describing the dimensions and nature of an object of class \code{"cladisticMatrix"} is printed to the console.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{compactify_cladistic_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Show print.cladisticMatrix version of each included data sets:
#' print.cladisticMatrix(x = day_2016)
#' print.cladisticMatrix(x = gauthier_1986)
#' print.cladisticMatrix(x = michaux_1989)
#'
#' @export print.cladisticMatrix
print.cladisticMatrix <- function(x, ...) {
  
  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = x, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
  # Gather basic data about the cladistic matrix:
  n_blocks <- length(x = x) - 1
  n_taxa <- nrow(x = x$matrix_1$matrix)
  block_sizes <- unname(obj = unlist(x = lapply(X = x[2:length(x)], FUN = function(x) ncol(x = x$matrix))))
  n_characters <- sum(block_sizes)
  data_types <- sort(x = tolower(x = unique(x = unname(obj = unlist(x = lapply(X = x[2:length(x = x)], FUN = function(x) x$datatype))))))
  character_ordering <- unname(obj = unlist(x = lapply(X = x[2:length(x = x)], FUN = function(x) x$ordering)))
  character_weights <- unname(obj = unlist(x = lapply(X = x[2:length(x = x)], FUN = function(x) x$character_weights)))
  n_unordered_characters <- sum(character_ordering == "unordered")
  n_ordered_characters <- sum(character_ordering == "ordered")
  n_continuous_characters <- sum(character_ordering == "continuous")
  n_costmatrix_characters <- length(x = grep(pattern = "step", x = character_ordering))
  
  # Calculate character sizes for ordering types:
  character_sizes <- nchar(x = c(n_unordered_characters, n_ordered_characters, n_continuous_characters, n_costmatrix_characters))
  character_spacings <- unlist(x = lapply(X = as.list(x = max(x = character_sizes) - character_sizes), function(x) paste(rep(x = " ", times = x), collapse = "")))
  
  # Get number of unique weights (excluding continuous characters):
  unique_non_continuous_weights <- unique(x = character_weights[character_ordering != "continuous"])

  # Create block parentheses text:
  if (length(x = block_sizes) == 1) block_parentheses <- ""
  if (length(x = block_sizes) == 2) block_parentheses <- paste0(" (in 2 matrix blocks of ", paste0(block_sizes, collapse = " and "), " characters, respectively)")
  if (length(x = block_sizes) > 2) block_parentheses <- paste0(" (in ", n_blocks, " matrix blocks of ", paste0(paste0(block_sizes[1:(length(x = block_sizes) - 1)], collapse = ", "), ", and ", block_sizes[length(x = block_sizes)]), " characters, respectively)")
  
  # Create block parentheses text:
  if (length(x = data_types) == 1) data_type_lines <- data_types
  if (length(x = data_types) == 2) data_type_lines <- paste0(paste0(data_types, collapse = " and "))
  if (length(x = data_types) > 2) data_type_lines <- paste0(paste0(data_types[1:(length(x = data_types) - 1)], collapse = ", "), ", and ", data_types[length(x = data_types)])
  
  # Make weights line:
  if (all(data_types == "continuous")) weight_line <- paste0("No non-continuous characters are present.\n")
  if (length(x = unique_non_continuous_weights) == 1) weight_line <- paste0("All non-continuous characters are weighted ", unique_non_continuous_weights, ".\n")
  if (length(x = unique_non_continuous_weights) > 1) weight_line <- paste0("Non-continuous characters have variable weights.\n")

  # Output:
  cat(paste0("Cladistic matrix containing ", n_taxa, " taxa and ", n_characters, " ", data_type_lines, " type characters", block_parentheses, ", of which:\n",
      paste0("  ", character_spacings[1], n_unordered_characters, " are unordered,\n"),
      paste0("  ", character_spacings[2], n_ordered_characters, " are ordered,\n"),
      paste0("  ", character_spacings[3], n_continuous_characters, " are continuous, and\n"),
      paste0("  ", character_spacings[4], n_costmatrix_characters, " are costmatrix characters\n"),
      weight_line)
  )
}
