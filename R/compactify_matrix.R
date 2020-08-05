#' Collapses matrix to unique character state distributions
#'
#' @description
#'
#' Collapses a cladistic matrix to just unique character state distributions and taxon names.
#'
#' @param cladistic_matrix The cladistic matrix in the format imported by \link{read_nexus_matrix}.
#' @param message Logical indicating whether or not a message should be printed to the screen if the matrix cannot be compactified.
#'
#' @details
#'
#' Important: not recommended for general use.
#'
#' This function is intended to make a matrix with redundant character state distributions smaller by collapsing these to single characters and upweighting them accordingly. It is intended purely for use with MRP matrices, but may have some very restricted uses elsewhere.
#'
#' The function also deletes any characters weighted zero from the matrix and will merge duplicate taxon names into unique character strings.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Examine the matrix pre-compactification:
#' michaux_1989$matrix_1$matrix
#'
#' # Examine the weights pre-compactification:
#' michaux_1989$matrix_1$character_weights
#'
#' # Compactify the matrix:
#' michaux_1989compact <- compactify_matrix(michaux_1989)
#'
#' # Examine the matrix post-compactification:
#' michaux_1989compact$matrix_1$matrix
#'
#' # Examine the weights post-compactification:
#' michaux_1989compact$matrix_1$character_weights
#' @export compactify_matrix
compactify_matrix <- function(cladistic_matrix, message = TRUE) {

  # FUTURE COULD CHECK FOR UNORD AND ORD WHEN BINARY AND HENCE MEANINGLESS

  # List any zero weight characters:
  zero_weight_characters <- which(x = unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")) == 0)

  # If there are zero weight characters then prune them:
  if (length(x = zero_weight_characters) > 0) cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix, characters2prune = zero_weight_characters)

  # For each matrix block:
  for (i in 2:length(x = cladistic_matrix)) {

    # Get strings for each character distribution, including ordering:
    character_distribution_strings <- paste(apply(cladistic_matrix[[i]]$matrix, 2, paste, collapse = ""), cladistic_matrix[[i]]$ordering, sep = " ")

    # Case if matrix can be compactified:
    if (length(x = unique(x = character_distribution_strings)) < length(x = character_distribution_strings) || any(duplicated(rownames(x = cladistic_matrix[[i]]$matrix)))) {

      # If collapsing characters because they are duplicated:
      if (length(x = unique(x = character_distribution_strings)) < length(x = character_distribution_strings)) {

        # Get rle of character distribution strings:
        rle_character_distribution_strings <- rle(sort(x = character_distribution_strings, decreasing = TRUE))

        # Set ordering of newly collapsed characters:
        cladistic_matrix[[i]]$ordering <- unlist(x = lapply(X = strsplit(rle_character_distribution_strings$values, " "), "[", 2))

        # Set character_weights of newly collapsed characters by aggregating weights of source characters:
        cladistic_matrix[[i]]$character_weights <- unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = as.list(x = rle_character_distribution_strings$values), "==", character_distribution_strings), which), function(x) cladistic_matrix[[i]]$character_weights[x]), sum))

        # Build new collapsed matrix:
        cladistic_matrix[[i]]$matrix <- matrix(unlist(x = lapply(X = lapply(X = strsplit(rle_character_distribution_strings$values, " "), "[", 1), strsplit, split = "")), nrow = nrow(cladistic_matrix[[i]]$matrix), dimnames = list(rownames(x = cladistic_matrix[[i]]$Matrix), c()))

        # Get ranges of values for characters in new collapsed matrix:
        character_ranges <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = apply(cladistic_matrix[[i]]$matrix, 2, strsplit, split = "/"), unlist), strsplit, split = "&"), unlist), unique), as.numeric), range)

        # Set new minimum values for collapsed matrix:
        cladistic_matrix[[i]]$minimum_values <- unlist(x = lapply(X = character_ranges, "[", 1))

        # Set new maximum values for collapsed matrix:
        cladistic_matrix[[i]]$maximum_values <- unlist(x = lapply(X = character_ranges, "[", 2))
      }

      # If collapsing matrix because taxa are duplicated:
      if (any(duplicated(rownames(x = cladistic_matrix[[i]]$matrix)))) {

        # Find duplicated taxa:
        duplicated_taxa <- unique(x = rownames(x = cladistic_matrix[[i]]$matrix)[duplicated(rownames(x = cladistic_matrix[[i]]$matrix))])

        # For each duplicated taxon:
        for (j in duplicated_taxa) {

          # Get rows for taxon:
          duplicate_rows <- which(x = rownames(x = cladistic_matrix[[i]]$matrix) == j)

          # Only continue if the duplicated rows are actually variable:
          if (length(x = unique(x = apply(cladistic_matrix[[i]]$matrix[duplicate_rows, ], 1, paste, collapse = ""))) > 1) {

            # Build duplicated matrix from other taxa:
            temporary_matrix <- matrix(rep(cladistic_matrix[[i]]$matrix[-duplicate_rows, ], length(x = duplicate_rows)), ncol = ncol(cladistic_matrix[[i]]$matrix) * length(x = duplicate_rows), dimnames = list(rownames(x = cladistic_matrix[[i]]$matrix)[-duplicate_rows], c()))

            # Add duplicated taxon as single row:
            temporary_matrix <- rbind(temporary_matrix, as.vector(t(cladistic_matrix[[i]]$matrix[duplicate_rows, ])))

            # Add duplicated taxon name:
            rownames(x = temporary_matrix)[nrow(temporary_matrix)] <- j

            # Update stored MRP matrix:
            cladistic_matrix[[i]]$matrix <- temporary_matrix

            # Update other parts of matrix:
            cladistic_matrix[[i]]$ordering <- rep(cladistic_matrix[[i]]$ordering, length(x = duplicate_rows))
            cladistic_matrix[[i]]$character_weights <- rep(cladistic_matrix[[i]]$character_weights, length(x = duplicate_rows))
            cladistic_matrix[[i]]$minimum_values <- rep(cladistic_matrix[[i]]$minimum_values, length(x = duplicate_rows))
            cladistic_matrix[[i]]$maximum_values <- rep(cladistic_matrix[[i]]$maximum_values, length(x = duplicate_rows))

            # BELOW IS EFFECTIVELY A RECURSION OF THIS FUNCTION!

            # Get strings for each character distribution, including ordering:
            secondary_distribution_strings <- paste(apply(cladistic_matrix[[i]]$matrix, 2, paste, collapse = ""), cladistic_matrix[[i]]$ordering, sep = " ")

            # If need to collapse characters because they are duplicated:
            if (length(x = unique(x = secondary_distribution_strings)) < length(x = secondary_distribution_strings)) {

              # Get rle of character distribution strings:
              rle_secondary_distribution_strings <- rle(sort(x = secondary_distribution_strings, decreasing = TRUE))

              # Set ordering of newly collapsed characters:
              cladistic_matrix[[i]]$ordering <- unlist(x = lapply(X = strsplit(rle_secondary_distribution_strings$values, " "), "[", 2))

              # Set weights of newly collapsed characters by aggregating weights of source characters:
              cladistic_matrix[[i]]$character_weights <- unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = as.list(x = rle_secondary_distribution_strings$values), "==", secondary_distribution_strings), which), function(x) cladistic_matrix[[i]]$character_weights[x]), sum))

              # Build new collapsed matrix:
              cladistic_matrix[[i]]$matrix <- matrix(unlist(x = lapply(X = lapply(X = strsplit(rle_secondary_distribution_strings$values, " "), "[", 1), strsplit, split = "")), nrow = nrow(cladistic_matrix[[i]]$matrix), dimnames = list(rownames(x = cladistic_matrix[[i]]$matrix), c()))

              # Get ranges of values for characters in new collapsed matrix:
              character_ranges <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = apply(cladistic_matrix[[i]]$matrix, 2, strsplit, split = "/"), unlist), strsplit, split = "&"), unlist), unique), as.numeric), range)

              # Set new minimum values for collapsed matrix:
              cladistic_matrix[[i]]$minimum_values <- unlist(x = lapply(X = character_ranges, "[", 1))

              # Set new maximum values for collapsed matrix:
              cladistic_matrix[[i]]$maximum_values <- unlist(x = lapply(X = character_ranges, "[", 2))
            }

            # If duplicated rows are not variable:
          } else {

            # Remove all but one duplicated row from the matrix:
            cladistic_matrix[[i]]$matrix <- cladistic_matrix[[i]]$matrix[-duplicate_rows[2:length(x = duplicate_rows)], , drop = FALSE]
          }
        }
      }

      # Case if matrix cannot be compactified:
    } else {

      # Print message to user:
      if (message) print("Matrix cannot be compactified. All character distributions are unique and character_weights are greater than zero.")
    }
  }

  # Output unaltered matrix:
  return(cladistic_matrix)
}
