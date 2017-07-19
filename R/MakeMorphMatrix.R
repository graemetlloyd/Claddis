#' Creates a morphological data file from a matrix
#' 
#' Creates a morphological data file from a character-taxon matrix.
#' 
#' Claddis generally assumes that matrices will be imported into R from the #NEXUS format, but in some cases (e.g., when using simulated data) it might be desirable to build a matrix within R. This function allows the user to convert such a matrix into the format required by other Claddis functions.
#'
#' NB: Currently the function cannot deal with step matrices or continuous characters.
#' 
#' @param CTmatrix A Character-Taxon (columns-rows) matrix, with taxa names as rownames.
#' @param header A scalar indicating any header text (defaults to an empty string: "").
#' @param weights A vector specifying the weights used (if not specified defaults to 1).
#' @param ordering A vector indicating whether characters are ordered ("ord") or unordered ("unord") (if no specified defaults to ordered).
#' @param symbols The symbols to use if writing to a file (defaults to the numbers 0:9 then the letters A to V).
#' @param equalise.weights Optional that overrides the weights specified above make all characters truly equally weighted.
#'
#' @return
#'
#' \item{header}{Any header text included in the file given between square brackets as character.}
#' \item{matrix}{A matrix of taxa (rows) and characters (columns). The matrix is in character format in order to deal with polymorphisms, which are separated by ampersands.}
#' \item{ordering}{A character vector of the same length as the number of morphological characters indicating whether they are ordered (\code{ord}) or unordered (\code{unord}).}
#' \item{weights}{A numeric vector of the same length as the number of morphological characters indicating their weights.}
#' \item{max.vals}{A numeric vector of the same length as the number of morphological characters indicating the maximum state values.}
#' \item{min.vals}{A numeric vector of the same length as the number of morphological characters indicating the minimum state values.}
#' \item{step.matrices}{A list of any step matrices supplied in the input file. Is \code{NULL} if none are specified.}
#' \item{symbols}{The original symbols used in the input data (these are replaced in the matrix by integers starting at zero.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{ReadMorphNexus}
#'
#' @examples
#'
#' # Create random 10-by-50 matrix:
#' CTmatrix <- matrix(sample(c("0", "1", "0&1", NA), 500,
#'   replace = TRUE), nrow = 10, dimnames =
#'   list(apply(matrix(sample(LETTERS, 40,
#'   replace = TRUE), nrow = 10), 1, paste,
#'   collapse = ""), c()))
#'
#' # Reformat for use elsewhere in Claddis:
#' MakeMorphMatrix(CTmatrix)
#'
#' @export MakeMorphMatrix
MakeMorphMatrix <- function(CTmatrix, header = "", weights = NULL, ordering = NULL, symbols = NULL, equalise.weights = FALSE) {
    
  # Check input is a matrix:
  if(!is.matrix(CTmatrix)) stop("CTmatrix must be a matrix.")
  
  # Check taxon names are supplied:
  if(is.null(rownames(CTmatrix))) stop("CTmatrix must have rownames indicating taxa.")
  
  # Check taxon names are unique (could cause downstream issues if not):
  if(any(duplicated(rownames(CTmatrix)))) stop("Taxon names must be unique.")
  
  # Delete any column names (could cause downstream issues otherwise):
  if(!is.null(colnames(CTmatrix))) colnames(CTmatrix) <- NULL
  
  # List any mystery character types:
  mystery.characters <- setdiff(unique(unlist(strsplit(as.character(unique(as.vector(CTmatrix))), "&"))), c(as.character(0:31), NA, "-"))
  
  # If mystery character types are present warn user:
  if(length(mystery.characters) > 0) stop("Characters must either be the integers 0 to 31, NA for missing, or & for polymorphisms.")

  # Check supplied weights are correct length:
  if(!is.null(weights) && length(weights) != ncol(CTmatrix)) stop("Weights must have same length as number of characters in CTmatrix.")
  
  # Check supplied weights are numeric:
  if(!is.null(weights) && !is.numeric(weights)) stop("Weights must be numeric.")
  
  # Check supplied weights are non-negative:
  if(!is.null(weights) && any(weights < 0)) stop("Weights must not be negative.")
  
  # Check supplied ordering is the correct length:
  if(!is.null(ordering) && length(ordering) != ncol(CTmatrix)) stop("Ordering must have same length as number of characters in CTmatrix.")
  
  # Check ordering is of unord or ord type only:
  if(!is.null(ordering) && length(setdiff(ordering, c("unord", "ord"))) > 0) stop("Ordering must be unord or ord only.")
  
  # Check symbols are of correct length:
  if(!is.null(symbols) && length(symbols) >= (diff(range(sort(as.numeric(unique(unlist(strsplit(as.character(unique(as.vector(CTmatrix))), "&"))))))) + 1)) stop("Symbols must be at least as long as the range of character values in CTmatrix.")

  # Check symbols are single characters only:
  if(!is.null(symbols) && any(nchar(symbols) != 1)) stop("Symbols must be single characters only.")

  # Check header is a single value:
  if(length(header) != 1)  stop("Header text must be a single value.")
  
  # If no ordering set default to ordered.
  if(is.null(ordering)) ordering <- rep("ord", ncol(CTmatrix))

  # If no weights are set:
  if(is.null(weights)) weights <- rep(1, ncol(CTmatrix))

  # Converting the inapplicables in the matrix if any
  if(any(CTmatrix == "-")) {
    CTmatrix_tmp <- gsub("-", NA, CTmatrix)
  } else {
    CTmatrix_tmp <- CTmatrix
  }

  # Calculate minimum values:
  min.vals <- unlist(lapply(lapply(lapply(lapply(apply(apply(CTmatrix_tmp, 2, as.character), 2, strsplit, split = "&"), unlist), as.numeric), sort), min))

  # Calculate maximum values:
  max.vals <- unlist(lapply(lapply(lapply(lapply(apply(apply(CTmatrix_tmp, 2, as.character), 2, strsplit, split = "&"), unlist), as.numeric), sort), max))

  # Default step matrices to NULL for now (may add this option in future):
  step.matrices <- NULL

  # If symbols were not set:
  if(is.null(symbols)) symbols <- c(c(0:9), LETTERS[1:22])[c(min(min.vals):max(max.vals)) + 1]

  # If wanting to equalise weights:
  if(equalise.weights) {
    
    # Get starting weights:
    weights <- apply(rbind(c(max.vals - min.vals), rep(1, nchar)), 2, max)
    
    # Updaye weights for unordered characters:
    weights[ordering == "unord"] <- 1
    
    # Update weights for ordered characters:
    weights[ordering == "ord"] <- 1 / weights[ordering == "ord"]
    
    # If there are step matrices:
    if(!is.null(step.matrices)) {
        
      # Get maximum distances for each step matrix:
      step.maxes <- unlist(lapply(lapply(step.matrices, as.numeric), max))
        
      # Update weights for step matrices:
      for(i in 1:length(step.maxes)) weights[ordering == names(step.matrices)[i]] <- 1 / step.maxes[i]
        
    }
    
    # Ensure all weights are integers by multiplying by product of all reciprocals:
    weights <- prod(unique(round(1 / weights))) * weights
    
    # Sub function to get all factors of an integer (stolen from: "http://stackoverflow.com/questions/6424856/r-function-for-returning-all-factors"):
    get.all.factors <- function(x) {
        
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
    out <- sort(unlist(apply(matrix(unique(weights)), 1, get.all.factors)))
    
    # As long as the maximum possible factor is greater than 1:
    while(max(rle(out)$values[rle(out)$lengths == length(unique(weights))]) > 1) {
        
      # Divide through weights by largest common factor:
      weights <- weights / max(rle(out)$values[rle(out)$lengths == length(unique(weights))])
        
      # Update factors for new weights:
      out <- sort(unlist(apply(matrix(unique(weights)), 1, get.all.factors)))
        
    }
    
  }

  # Create output formatted data:
  result <- list(header, CTmatrix, ordering, weights, max.vals, min.vals, step.matrices, symbols)

  # Add names to result:
  names(result) <- c("header", "matrix", "ordering", "weights", "max.vals", "min.vals", "step.matrices", "symbols")

  # Return result:
  return(result)

}
