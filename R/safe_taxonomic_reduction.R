#' Safe Taxonomic Reduction
#'
#' @description
#'
#' Performs Safe Taxonomic Reduction (STR) on a character-taxon matrix.
#'
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#'
#' @details
#'
#' Performs Safe Taxonomic Reduction (Wilkinson 1995).
#'
#' If no taxa can be safely removed will print the text "No taxa can be safely removed", and the \code{str_taxa} and \code{removed_matrix} will have no rows.
#'
#' NB: If your data contains inapplicable characters these will be treated as missing data, but this is inappropriate. Thus the user is advised to double check that any removed taxa make sense in the light of inapplicable states. (As far as I am aware this same behaviour occurs in the TAXEQ3 software.)
#'
#' @return
#'
#' \item{str_taxa}{A matrix listing the taxa that can be removed (\code{junior}), the taxa which they are equivalent to (\code{senior}) and the rule under which they can be safely removed (\code{rule}).}
#' \item{reduced_matrix}{A character-taxon matrix excluding the taxa that can be safely removed.}
#' \item{removed_matrix}{A character-taxon matrix of the taxa that can be safely removed.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references Wilkinson, M., 1995. Coping with abundant missing entries in phylogenetic inference using parsimony. \emph{Systematic Biology}, \bold{44}, 501-514.
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{compactify_matrix}, \link{prune_cladistic_matrix}, \link{safe_taxonomic_reinsertion}, \link{read_nexus_matrix}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Performs STR on the Gauthier 1986 dataset used in Wilkinson (1995):
#' str_data <- safe_taxonomic_reduction(cladistic_matrix = gauthier_1986)
#'
#' # View deleted taxa:
#' str_data$str_taxa
#'
#' # View reduced matrix:
#' str_data$reduced_matrix
#'
#' # View removed matrix:
#' str_data$removed_matrix
#' @export safe_taxonomic_reduction
safe_taxonomic_reduction <- function(cladistic_matrix) {

  # Store unaltered version of matrix to return to later:
  clean_matrix <- cladistic_matrix

  # If there is more than one matrix block (will need to temporarily combine into single matrix):
  if (length(x = cladistic_matrix) > 2) {

    # Get list of just the matrix blocks:
    matrix_blocks <- lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix")

    # For each additional block add at end of first block:
    for (i in 2:length(x = matrix_blocks)) matrix_blocks[[1]] <- cbind(matrix_blocks[[1]], matrix_blocks[[i]][rownames(x = matrix_blocks[[1]]), ])

    # Store single matrix block as first block of cladistic_matrix:
    cladistic_matrix[[2]]$matrix <- matrix_blocks[[1]]

    # Store ordering for all blocks in first block of cladistic_matrix:
    cladistic_matrix[[2]]$ordering <- as.vector(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))

    # Store weights for all blocks in first block of cladistic_matrix:
    cladistic_matrix[[2]]$character_weights <- as.vector(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))

    # Store minimum values for all blocks in first block of cladistic_matrix:
    cladistic_matrix[[2]]$minimum_values <- as.vector(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))

    # Store maximum values for all blocks in first block of cladistic_matrix:
    cladistic_matrix[[2]]$maximum_values <- as.vector(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))

    # Isolate to just first matrix block:
    cladistic_matrix <- cladistic_matrix[1:2]
  }

  # Prune out any zero weight characters, if they exist:
  if (any(cladistic_matrix[[2]]$character_weights == 0)) cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix = cladistic_matrix, characters2prune = which(x = cladistic_matrix[[2]]$character_weights == 0))

  # Order matrix from least to most complete taxon (as least is most likely to be removed):
  cladistic_matrix[[2]]$matrix <- cladistic_matrix[[2]]$matrix[order(apply(apply(cladistic_matrix[[2]]$matrix, 1, is.na), 2, sum), decreasing = TRUE), ]

  # Subfunction to be used by lapply to check downwards (look at more complete taxa only) for safely removable taxa:
  check_downward_for_matches <- function(rownumber, cladistic_matrix) {

    # First find which characters are not scored as missing in current taxon:
    non_missing_characters <- !is.na(cladistic_matrix[[2]]$matrix[rownumber, ])

    # Build isolated matrix block from current taxon to end of matrix only for characters coded for current taxon:
    matrix_block_to_check <- cladistic_matrix[[2]]$matrix[rownumber:nrow(cladistic_matrix[[2]]$matrix), non_missing_characters, drop = FALSE]

    # Find any taxa that have missing characters inside the block (can not be true parents):
    taxa_with_missing_characters <- names(which(x = apply(apply(matrix_block_to_check, 1, is.na), 2, any)))

    # If any taxa have missing characters inside the block (can not be true parents) remove them from the block:
    if (length(x = taxa_with_missing_characters) > 0) matrix_block_to_check <- matrix_block_to_check[-match(taxa_with_missing_characters, rownames(x = matrix_block_to_check)), , drop = FALSE]

    # Set start column (first character in matrix block):
    start_column <- 1

    # Set end column (last character in matrix block):
    end_column <- ncol(matrix_block_to_check)

    # As long as there is more potential seniors (more than one row in the matrix block) and there are still characters to check:
    while (nrow(matrix_block_to_check) > 1 && start_column != (end_column + 1)) {

      # Only look at rows where taxa are coded for the ith character:
      non_missing_rows <- which(x = !is.na(matrix_block_to_check[2:nrow(matrix_block_to_check), start_column]))

      # List any insafe taxa (i.e., those with a diffrent coding to the current taxon):
      unsafe_taxa <- names(which(x = matrix_block_to_check[(non_missing_rows + 1), start_column] != matrix_block_to_check[1, start_column]))

      # Reduce the matrix block to just potential safe taxon parents (and the current taxon):
      if (length(x = unsafe_taxa) > 0) matrix_block_to_check <- matrix_block_to_check[-match(unsafe_taxa, rownames(x = matrix_block_to_check)), , drop = FALSE]

      # Iterate to next colum:
      start_column <- start_column + 1
    }

    # If safe taxonomic reduction is possible store junior and senior(s):
    if (nrow(matrix_block_to_check) > 1) {
      return(cbind(rep(x = rownames(x = matrix_block_to_check)[1], times = nrow(matrix_block_to_check) - 1), rownames(x = matrix_block_to_check)[2:nrow(matrix_block_to_check)]))
    }
  }

  # Get any safely removable taxa and their senior parents:
  safe_taxa <- lapply(X = as.list(x = 1:(nrow(cladistic_matrix[[2]]$matrix) - 1)), check_downward_for_matches, cladistic_matrix = cladistic_matrix)

  # If no taxa can be safely deleted:
  if (all(unlist(x = lapply(X = safe_taxa, is.null)))) {

    # Warn user:
    print("No taxa can be safely removed")

    # Create empty safe to remove matrix:
    safe_taxa <- matrix(nrow = 0, ncol = 3, dimnames = list(c(), c("junior", "senior", "rule")))

    # Set reduced matrix as unaltered matrix:
    removed_matrix <- reduced_matrix <- clean_matrix

    # Set removed matrix as empty matrix:
    removed_matrix <- prune_cladistic_matrix(clean_matrix, taxa2prune = rownames(x = clean_matrix$matrix_1$matrix))

    # If there are taxa that can be deleted:
  } else {

    # Collapse down to just those with data:
    safe_taxa <- safe_taxa[which(x = !unlist(x = lapply(X = safe_taxa, is.null)))]

    # Rebuild as matrix:
    safe_taxa <- cbind(unlist(x = lapply(X = safe_taxa, "[", , 1)), unlist(x = lapply(X = safe_taxa, "[", , 2)))

    # Add column names:
    colnames(x = safe_taxa) <- c("junior", "senior")

    # Create empty vector to store rule data:
    str_rule <- c()

    # For each STR pair:
    for (i in 1:nrow(safe_taxa)) {

      # Check for rule 1 (symmetrical coding):
      if (paste(cladistic_matrix$matrix_1$matrix[safe_taxa[i, 1], ], collapse = "") == paste(cladistic_matrix$matrix_1$matrix[safe_taxa[i, 2], ], collapse = "")) {

        # Check for any missing data:
        if (any(is.na(cladistic_matrix$matrix_1$matrix[safe_taxa[i, 1], ]))) {

          # Is rule 1b:
          str_rule <- c(str_rule, "rule_1b")

          # If no missing data:
        } else {

          # Is rule 1a:
          str_rule <- c(str_rule, "rule_1a")
        }

        # Must be rule 2 (asymmetric coding):
      } else {

        # Check for any missing data:
        if (any(is.na(cladistic_matrix$matrix_1$matrix[safe_taxa[i, 1], ]))) {

          # Is rule 2B:
          str_rule <- c(str_rule, "rule_2b")

          # If no missing data is rule 1a:
        } else {

          # Is rule 2A:
          str_rule <- c(str_rule, "rule_2a")
        }
      }
    }

    # Add rule to output:
    safe_taxa <- cbind(safe_taxa, str_rule)

    # Update column heading for rule:
    colnames(x = safe_taxa)[3] <- "rule"

    # Set reduced matrix as complete matrix:
    reduced_matrix <- prune_cladistic_matrix(clean_matrix, taxa2prune = unique(x = safe_taxa[, "junior"]))

    # Set removed matrix as empty matrix:
    removed_matrix <- prune_cladistic_matrix(clean_matrix, taxa2prune = setdiff(x = rownames(x = clean_matrix$matrix_1$matrix), y = unique(x = safe_taxa[, "junior"])))
  }

  # Invisibly return output compiled into a single list:
  invisible(list(str_taxa = safe_taxa, reduced_matrix = reduced_matrix, removed_matrix = removed_matrix))
}
