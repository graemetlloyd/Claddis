#' Prunes a character matrix of characters or taxa
#'
#' @description
#'
#' Prunes a character matrix of characters, taxa, or both.
#'
#' @param cladistic_matrix The cladistic matrix in the format imported by \link{read_nexus_matrix}.
#' @param blocks2prune A vector of number(s) of any blocks to prune.
#' @param characters2prune A vector of character numbers to prune.
#' @param taxa2prune A vector of taxon names to prune (these must be present in \code{rownames(x = cladistic_matrix$matrix}).
#' @param remove_invariant A logical for whether invariant characters should (TRUE) or should not (FALSE, default) be pruned.
#'
#' @details
#'
#' Removing characters or taxa from a matrix imported using \link{read_nexus_matrix} is not simple due to associated vectors for ordering, character weights etc. To save repetitively pruning each part this function takes the matrix as input and vector(s) of either block numbers, character numbers, taxon names, or any combination thereof and returns a matrix with these items removed. Minimum and maximum values (used by \link{calculate_morphological_distances}) are also updated and the user has the option to remove constant characters this way as well (e.g, to reduce the memory required for a DNA matrix).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{compactify_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Remove the outgroup taxon and characters 11 and 53 from gauthier_1986:
#' prunedmatrix <- prune_cladistic_matrix(
#'   cladistic_matrix =
#'     gauthier_1986, characters2prune = c(11, 53), taxa2prune =
#'     c("Outgroup")
#' )
#'
#' # Show priuned matrix:
#' prunedmatrix$matrix_1$matrix
#' @export prune_cladistic_matrix
prune_cladistic_matrix <- function(cladistic_matrix, blocks2prune = c(), characters2prune = c(), taxa2prune = c(), remove_invariant = FALSE) {

  # How do blocks and characters to prune interact? (explain to user in manual)

  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")

  # Subfunction to find length of character types for each character (i.e., unique values excluding polymorphisms but included inapplicables):
  find_length <- function(x) {

    # Convert each column of matrix to a list of numeric values:
    x <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = apply(x, 2, as.list), unlist), strsplit, split = "&|/"), unlist), unique), sort), length)

    # Return(x):
    return(x)
  }

  # Check that something to prune has been specified:
  if (is.null(blocks2prune) && is.null(characters2prune) && is.null(taxa2prune) && remove_invariant == FALSE) stop("No blocks, taxa, or characters to prune specified.")

  # Check blocks specified exist and stop and warn :
  if (length(x = setdiff(x = blocks2prune, y = 1:(length(x = cladistic_matrix) - 1))) > 0) stop("Block numbers specified that are not present in data.")

  # Check characters specified exist and stop and warn user if not:
  if (length(x = setdiff(x = characters2prune, y = 1:sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol))))) > 0) stop("characters specified that are outside the scope of the matrix. Check and retry.")

  # Check taxa specified exist and stop and warn user if not:
  if (length(x = setdiff(x = taxa2prune, y = rownames(x = cladistic_matrix$matrix_1$matrix))) > 0) stop("Taxa specified that are not found in the matrix. Check and retry.")

  # Get number of characters:
  n_characters <- sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)))

  # If there are characters to prune:
  if (!is.null(characters2prune)) {

    # Get character blocks for each character in descendng order (as want to work backwards so things match up properly):
    character_blocks <- unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = sort(x = characters2prune, decreasing = TRUE)), ">", cumsum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)))), which), length)) + 1

    # Initial build of characters in list form:
    characters_as_list <- lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), function(x) 1:ncol(x))

    # Actually form list of character numbers (i.e., renumber characters in second or higher blocks):
    if (length(x = characters_as_list) > 1) for (i in 2:length(x = characters_as_list)) characters_as_list[[i]] <- characters_as_list[[i]] + max(characters_as_list[[(i - 1)]])

    # For each unique character block:
    for (i in unique(x = character_blocks)) {

      # Find columns to delete in ith matrix:
      columns_to_delete <- match(sort(x = characters2prune, decreasing = TRUE)[character_blocks == i], characters_as_list[[i]])

      # Remove characters from matrix:
      cladistic_matrix[[(i + 1)]]$matrix <- cladistic_matrix[[(i + 1)]]$matrix[, -columns_to_delete, drop = FALSE]

      # Remove characters from ordering:
      cladistic_matrix[[(i + 1)]]$ordering <- cladistic_matrix[[(i + 1)]]$ordering[-columns_to_delete]

      # Remove characters from weights:
      cladistic_matrix[[(i + 1)]]$character_weights <- cladistic_matrix[[(i + 1)]]$character_weights[-columns_to_delete]

      # Remove characters from minimum values:
      cladistic_matrix[[(i + 1)]]$minimum_values <- cladistic_matrix[[(i + 1)]]$minimum_values[-columns_to_delete]

      # Remove characters from maximum values:
      cladistic_matrix[[(i + 1)]]$maximum_values <- cladistic_matrix[[(i + 1)]]$maximum_values[-columns_to_delete]
    }
  }

  # If there are taxa to prune:
  if (!is.null(taxa2prune)) {

    # Remove pruned taxa from each block:
    for (i in 2:length(x = cladistic_matrix)) cladistic_matrix[[i]]$matrix <- cladistic_matrix[[i]]$matrix[-match(taxa2prune, rownames(x = cladistic_matrix[[i]]$matrix)), , drop = FALSE]
  }

  # If there are blocks to prune:
  if (!is.null(blocks2prune)) {

    # Remove blocks to be rpuned:
    cladistic_matrix <- cladistic_matrix[-(blocks2prune + 1)]

    # Rename (renumber) remaining matrix blocks:
    names(cladistic_matrix[2:length(x = cladistic_matrix)]) <- paste("matrix_", 1:(length(x = cladistic_matrix) - 1), sep = "")
  }

  # If there are invariant characters:
  if (remove_invariant) {

    # Find any invariant characters:
    invariants_as_list <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), find_length), unlist), "<", 1), which)

    # If there are invariant characters:
    if (length(x = unlist(x = invariants_as_list)) > 0) {

      # For each matrix block:
      for (i in 1:length(x = invariants_as_list)) {

        # Only if there are invraints for this block:
        if (length(x = invariants_as_list[[i]]) > 0) {

          # Remove characters from matrix:
          cladistic_matrix[[(i + 1)]]$matrix <- cladistic_matrix[[(i + 1)]]$matrix[, -invariants_as_list[[i]], drop = FALSE]

          # Remove characters from ordering:
          cladistic_matrix[[(i + 1)]]$ordering <- cladistic_matrix[[(i + 1)]]$ordering[-invariants_as_list[[i]]]

          # Remove characters from weights:
          cladistic_matrix[[(i + 1)]]$character_weights <- cladistic_matrix[[(i + 1)]]$character_weights[-invariants_as_list[[i]]]

          # Remove characters from minimum values:
          cladistic_matrix[[(i + 1)]]$minimum_values <- cladistic_matrix[[(i + 1)]]$minimum_values[-invariants_as_list[[i]]]

          # Remove characters from maximum values:
          cladistic_matrix[[(i + 1)]]$maximum_values <- cladistic_matrix[[(i + 1)]]$maximum_values[-invariants_as_list[[i]]]
        }
      }
    }
  }

  # Check for empty blocks and store them as blocks to delete if found:
  new_blocks_to_delete <- which(x = unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)) == 0)

  # If there are new blocks to prune:
  if (length(x = new_blocks_to_delete) > 0) {

    # Remove blocks to be rpuned:
    cladistic_matrix <- cladistic_matrix[-(new_blocks_to_delete + 1)]
  }
  
  # Rename (renumber) matrix blocks to ensure consistent output:
  names(cladistic_matrix[2:length(x = cladistic_matrix)]) <- paste("matrix_", 1:(length(x = cladistic_matrix) - 1), sep = "")

  # Ensure class is set:
  class(cladistic_matrix) <- "cladisticMatrix"

  # Return pruned matrix:
  cladistic_matrix
}
