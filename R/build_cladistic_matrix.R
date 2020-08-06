#' Creates a morphological data file from a matrix
#'
#' @description
#'
#' Creates a morphological data file from a character-taxon matrix.
#'
#' @param character_taxon_matrix A Character-Taxon (columns-rows) matrix, with taxon names as rownames.
#' @param header A scalar indicating any header text (defaults to an empty string: "").
#' @param character_weights A vector specifying the weights used (if not specified defaults to 1).
#' @param ordering A vector indicating whether characters are ordered ("ord") or unordered ("unord") (if no specified defaults to ordered).
#' @param symbols The symbols to use if writing to a file (defaults to the numbers 0:9 then the letters A to V).
#' @param equalise.weights Optional that overrides the weights specified above make all characters truly equally weighted.
#' @param ignore_duplicate_taxa Logical indicating whether or not to ignore (allow; TRUE) duplicate taxa or not (FALSE; default).
#'
#' @details
#'
#' Claddis generally assumes that matrices will be imported into R from the #NEXUS format, but in some cases (e.g., when using simulated data) it might be desirable to build a matrix within R. This function allows the user to convert such a matrix into the format required by other Claddis functions as long as it only contains a single block.
#'
#' NB: Currently the function cannot deal directly with step matrices or continuous characters.
#'
#' @return
#'
#' \item{topper}{Contains any header text or step matrices and pertains to the entire file.}
#' \item{matrix_N}{One or more matrix blocks (numbered 1 to N) with associated information pertaining only to that matrix block. This includes the block name (if specificed, NA if not), the block datatype (one of "CONTINUOUS", "DNA", "NUCLEOTIDE", "PROTEIN", "RESTRICTION", "RNA", or "STANDARD"), the actual matrix (taxa as rows, names stored as rownames and characters as columns), the ordering type of each character ("ord" = ordered, "unord" = unordered), the character weights, the minimum and maximum values (used by Claddis' distance functions), and the original characters (symbols, missing, and gap values) used for writing out the data.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{compactify_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Create random 10-by-50 matrix:
#' character_taxon_matrix <- matrix(sample(c("0", "1", "0&1", NA, ""),
#'   500,
#'   replace = TRUE
#' ),
#' nrow = 10, dimnames =
#'   list(apply(matrix(sample(LETTERS, 40,
#'     replace = TRUE
#'   ), nrow = 10), 1, paste,
#'   collapse = ""
#'   ), c())
#' )
#'
#' # Reformat for use elsewhere in Claddis:
#' build_cladistic_matrix(character_taxon_matrix)
#' @export build_cladistic_matrix
build_cladistic_matrix <- function(character_taxon_matrix, header = "", character_weights = NULL, ordering = NULL, symbols = NULL, equalise.weights = FALSE, ignore_duplicate_taxa = FALSE) {

  # Use manual to inform user about defaults for weighting and ordering etc.

  # Check input is a matrix:
  if (!is.matrix(character_taxon_matrix)) stop("character_taxon_matrix must be a matrix.")

  # Check taxon names are supplied:
  if (is.null(rownames(x = character_taxon_matrix))) stop("character_taxon_matrix must have rownames indicating taxa.")

  # Check taxon names are unique (could cause downstream issues if not):
  if (!ignore_duplicate_taxa) if (any(duplicated(rownames(x = character_taxon_matrix)))) stop("Taxon names must be unique.")

  # Delete any column names (could cause downstream issues otherwise):
  if (!is.null(colnames(x = character_taxon_matrix))) colnames(x = character_taxon_matrix) <- NULL




  # WILL ONLY WORK FOR DISCRETE CHARACTERS!




  # List any mystery character types:
  mystery_characters <- setdiff(x = unique(x = unlist(x = strsplit(as.character(unique(x = as.vector(character_taxon_matrix))), split = "&|/"))), y = c(as.character(0:31), NA))






  # If mystery character types are present warn user:
  if (length(x = mystery_characters) > 0) stop("characters must either be the integers 0 to 31, NA for missing, & for polymorphisms, or / for uncertainties.")

  # Check supplied weights are correct length:
  if (!is.null(character_weights) && length(x = character_weights) != ncol(character_taxon_matrix)) stop("character_weights must have same length as number of characters in character_taxon_matrix.")

  # Check supplied weights are numeric:
  if (!is.null(character_weights) && !is.numeric(character_weights)) stop("character_weights must be numeric.")

  # Check supplied weights are non-negative:
  if (!is.null(character_weights) && any(character_weights < 0)) stop("character_weights must not be negative.")

  # Check supplied ordering is the correct length:
  if (!is.null(ordering) && length(x = ordering) != ncol(character_taxon_matrix)) stop("ordering must have same length as number of characters in character_taxon_matrix.")

  # Check ordering is of unord or ord type only:
  if (!is.null(ordering) && length(x = setdiff(x = ordering, y = c("unord", "ord"))) > 0) stop("ordering must be unord or ord only.")

  # Check symbols are of correct length:
  if (!is.null(symbols) && length(x = symbols) >= (diff(range(as.numeric(unique(x = sort(x = unlist(x = strsplit(as.vector(character_taxon_matrix), split = "&|/"))))))) + 1)) stop("symbols must be at least as long as the range of character values in character_taxon_matrix.")

  # Check symbols are single characters only:
  if (!is.null(symbols) && any(nchar(x = symbols) != 1)) stop("symbols must be single characters only.")

  # Check header is a single value:
  if (length(x = header) != 1) stop("header text must be a single value.")

  # If no ordering set default to ordered.
  if (is.null(ordering)) ordering <- rep("ord", ncol(character_taxon_matrix))

  # If no weights are set:
  if (is.null(character_weights)) character_weights <- rep(1, ncol(character_taxon_matrix))

  # Calculate minimum values:
  minimum_values <- apply(character_taxon_matrix, 2, function(x) sort(x = as.numeric(unlist(x = strsplit(x, split = "&|/"))), decreasing = FALSE)[1])

  # Calculate maximum values:
  maximum_values <- apply(character_taxon_matrix, 2, function(x) sort(x = as.numeric(unlist(x = strsplit(x, split = "&|/"))), decreasing = TRUE)[1])

  # If any NAs in minimum_values replace with zero:
  if (any(is.na(minimum_values))) minimum_values[is.na(minimum_values)] <- 0

  # If any NAs in maximum_values replace with zero:
  if (any(is.na(maximum_values))) maximum_values[is.na(maximum_values)] <- 0

  # Default step matrices to NULL for now (may add this option in future):
  step_matrices <- NULL

  # If symbols were not set:
  if (is.null(symbols)) symbols <- c(c(0:9), LETTERS[1:22])[c(min(minimum_values):max(maximum_values)) + 1]

  # If wanting to equalise weights:
  if (equalise.weights) {

    # Get starting weights:
    character_weights <- apply(rbind(c(maximum_values - minimum_values), rep(1, nchar)), 2, max)

    # Update weights for unordered characters:
    character_weights[ordering == "unord"] <- 1

    # Update weights for ordered characters:
    character_weights[ordering == "ord"] <- 1 / character_weights[ordering == "ord"]

    # If there are step matrices (not technically using these yet, but I guess this will have to exist eventually):
    if (!is.null(step_matrices)) {

      # Get maximum distances for each step matrix:
      step_maxima <- unlist(x = lapply(X = lapply(X = step_matrices, as.numeric), max))

      # Update weights for step matrices:
      for (i in 1:length(x = step_maxima)) character_weights[ordering == names(step_matrices)[i]] <- 1 / step_maxima[i]
    }

    # Ensure all weights are integers by multiplying by product of all reciprocals:
    character_weights <- prod(unique(x = round(1 / character_weights))) * character_weights

    # Sub function to get all factors of an integer (stolen from: "http://stackoverflow.com/questions/6424856/r-function-for-returning-all-factors"):
    get_all_factors <- function(x) {

      # Ensure input is an integer:
      x <- as.integer(x)

      # Ensure x is positive and get sequence of 1 to x:
      div <- seq_len(abs(x))

      # Get factors of x (i.e. numbers whose remainders are zero):
      factors <- div[x %% div == 0L]

      # Output answer:
      return(factors)
    }

    # Get factors of every weight currently applied:
    out <- sort(x = unlist(x = apply(matrix(unique(x = character_weights)), 1, get_all_factors)))

    # As long as the maximum possible factor is greater than 1:
    while (max(rle(out)$values[rle(out)$lengths == length(x = unique(x = character_weights))]) > 1) {

      # Divide through weights by largest common factor:
      character_weights <- character_weights / max(rle(out)$values[rle(out)$lengths == length(x = unique(x = character_weights))])

      # Update factors for new weights:
      out <- sort(x = unlist(x = apply(matrix(unique(x = character_weights)), 1, get_all_factors)))
    }
  }

  # Build matrix topper:
  topper <- list(header = header, step_matrices = step_matrices)

  # Build characters list:
  characters <- list(symbols = symbols, missing = "?", gap = "-")

  # Build matrix_1 list:
  matrix_1 <- list(block_name = NA, datatype = "STANDARD", matrix = character_taxon_matrix, ordering = ordering, character_weights = character_weights, minimum_values = minimum_values, maximum_values = maximum_values, characters = characters)

  # Return assimilated output:
  list(topper = topper, matrix_1 = matrix_1)
}
