#' Creates a morphological data file from a matrix
#'
#' @description
#'
#' Creates a morphological data file from a character-taxon matrix.
#' 
#' @param CharacterTaxonMatrix A Character-Taxon (columns-rows) matrix, with taxon names as rownames.
#' @param header A scalar indicating any header text (defaults to an empty string: "").
#' @param weights A vector specifying the weights used (if not specified defaults to 1).
#' @param ordering A vector indicating whether characters are ordered ("ord") or unordered ("unord") (if no specified defaults to ordered).
#' @param symbols The symbols to use if writing to a file (defaults to the numbers 0:9 then the letters A to V).
#' @param equalise.weights Optional that overrides the weights specified above make all characters truly equally weighted.
#'
#' @details
#'
#' Claddis generally assumes that matrices will be imported into R from the #NEXUS format, but in some cases (e.g., when using simulated data) it might be desirable to build a matrix within R. This function allows the user to convert such a matrix into the format required by other Claddis functions as long as it only contains a single block.
#'
#' NB: Currently the function cannot deal directly with step matrices or continuous characters.
#'
#' @return
#'
#' \item{Topper}{Contains any header text or step matrices and pertains to the entire file.}
#' \item{Matrix_N}{One or more matrix blocks (numbered 1 to N) with associated information pertaining only to that matrix block. This includes the block name (if specificed, NA if not), the block datatype (one of "CONTINUOUS", "DNA", "NUCLEOTIDE", "PROTEIN", "RESTRICTION", "RNA", or "STANDARD"), the actual matrix (taxa as rows, names stored as rownames and characters as columns), the ordering type of each character ("ord" = ordered, "unord" = unordered), the character weights, the minimum and maximum values (used by Claddis' distance functions), and the original characters (symbols, missing, and gap values) used for writing out the data.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{ReadMorphNexus}
#'
#' @examples
#'
#' # Create random 10-by-50 matrix:
#' CharacterTaxonMatrix <- matrix(sample(c("0", "1", "0&1", NA, ""),
#'   500, replace = TRUE), nrow = 10, dimnames =
#'   list(apply(matrix(sample(LETTERS, 40,
#'   replace = TRUE), nrow = 10), 1, paste,
#'   collapse = ""), c()))
#'
#' # Reformat for use elsewhere in Claddis:
#' MakeMorphMatrix(CharacterTaxonMatrix)
#'
#' @export MakeMorphMatrix
MakeMorphMatrix <- function(CharacterTaxonMatrix, header = "", weights = NULL, ordering = NULL, symbols = NULL, equalise.weights = FALSE) {
  
  # Check input is a matrix:
  if(!is.matrix(CharacterTaxonMatrix)) stop("CharacterTaxonMatrix must be a matrix.")
  
  # Check taxon names are supplied:
  if(is.null(rownames(CharacterTaxonMatrix))) stop("CharacterTaxonMatrix must have rownames indicating taxa.")
  
  # Check taxon names are unique (could cause downstream issues if not):
  if(any(duplicated(rownames(CharacterTaxonMatrix)))) stop("Taxon names must be unique.")
  
  # Delete any column names (could cause downstream issues otherwise):
  if(!is.null(colnames(CharacterTaxonMatrix))) colnames(CharacterTaxonMatrix) <- NULL
  
  
  
  
  # WILL ONLY WORK FOR DISCRETE CHARACTERS!
  
  
  
  
  # List any mystery character types:
  mystery.characters <- setdiff(unique(unlist(strsplit(as.character(unique(as.vector(CharacterTaxonMatrix))), split = "&|/"))), c(as.character(0:31), NA))
  
  
  
  
  
  
  # If mystery character types are present warn user:
  if(length(mystery.characters) > 0) stop("Characters must either be the integers 0 to 31, NA for missing, & for polymorphisms, or / for uncertainties.")

  # Check supplied weights are correct length:
  if(!is.null(weights) && length(weights) != ncol(CharacterTaxonMatrix)) stop("Weights must have same length as number of characters in CharacterTaxonMatrix.")
  
  # Check supplied weights are numeric:
  if(!is.null(weights) && !is.numeric(weights)) stop("Weights must be numeric.")
  
  # Check supplied weights are non-negative:
  if(!is.null(weights) && any(weights < 0)) stop("Weights must not be negative.")
  
  # Check supplied ordering is the correct length:
  if(!is.null(ordering) && length(ordering) != ncol(CharacterTaxonMatrix)) stop("Ordering must have same length as number of characters in CharacterTaxonMatrix.")
  
  # Check ordering is of unord or ord type only:
  if(!is.null(ordering) && length(setdiff(ordering, c("unord", "ord"))) > 0) stop("Ordering must be unord or ord only.")
  
  # Check symbols are of correct length:
  if(!is.null(symbols) && length(symbols) >= (diff(range(sort(as.numeric(unique(unlist(strsplit(as.character(unique(as.vector(CharacterTaxonMatrix))), "&"))))))) + 1)) stop("Symbols must be at least as long as the range of character values in CharacterTaxonMatrix.")

  # Check symbols are single characters only:
  if(!is.null(symbols) && any(nchar(symbols) != 1)) stop("Symbols must be single characters only.")

  # Check header is a single value:
  if(length(header) != 1)  stop("Header text must be a single value.")
  
  # If no ordering set default to ordered.
  if(is.null(ordering)) ordering <- rep("ord", ncol(CharacterTaxonMatrix))

  # If no weights are set:
  if(is.null(weights)) weights <- rep(1, ncol(CharacterTaxonMatrix))

  # Calculate minimum values:
  min.vals <- unlist(lapply(lapply(lapply(lapply(apply(apply(CharacterTaxonMatrix, 2, as.character), 2, strsplit, split = "&|/"), unlist), as.numeric), sort), min))

  # Calculate maximum values:
  max.vals <- unlist(lapply(lapply(lapply(lapply(apply(apply(CharacterTaxonMatrix, 2, as.character), 2, strsplit, split = "&|/"), unlist), as.numeric), sort), max))

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
    
    # If there are step matrices (not technically using these yet, but I guess this will have to exist eventually):
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
  
  # Build matrix topper:
  Topper <- list(header, step.matrices)
  
  # Add names to topper:
  names(Topper) <- c("Header", "StepMatrices")
  
  # Build characters list:
  Characters <- list(symbols, "?", "-")
  
  # Add names to characters list:
  names(Characters) <- c("Symbols", "Missing", "Gap")
  
  # Build Matrix_1 list:
  Matrix_1 <- list(NA, "STANDARD", CharacterTaxonMatrix, ordering, weights, min.vals, max.vals, Characters)
  
  # Add names to Matrix_1:
  names(Matrix_1) <- c("BlockName", "Datatype", "Matrix", "Ordering", "Weights", "MinVals", "MaxVals", "Characters")

  # Assimilate into output:
  result <- list(Topper, Matrix_1)

  # Add names to output:
  names(result) <- c("Topper", "Matrix_1")

  # Return output:
  return(result)

}
