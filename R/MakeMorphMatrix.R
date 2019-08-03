#' Creates a morphological data file from a matrix
#'
#' @description
#'
#' Creates a morphological data file from a character-taxon matrix.
#' 
#' @param CharacterTaxonMatrix A Character-Taxon (columns-rows) matrix, with taxon names as rownames.
#' @param header A scalar indicating any header text (defaults to an empty string: "").
#' @param Weights A vector specifying the weights used (if not specified defaults to 1).
#' @param ordering A vector indicating whether characters are ordered ("ord") or unordered ("unord") (if no specified defaults to ordered).
#' @param symbols The symbols to use if writing to a file (defaults to the numbers 0:9 then the letters A to V).
#' @param equalise.weights Optional that overrides the weights specified above make all characters truly equally weighted.
#' @param ignore.duplicate.taxa Logical indicating whether or not to ignore (allow; TRUE) duplicate taxa or not (FALSE; default).
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
MakeMorphMatrix <- function(CharacterTaxonMatrix, header = "", Weights = NULL, ordering = NULL, symbols = NULL, equalise.weights = FALSE, ignore.duplicate.taxa = FALSE) {
  
  # Check input is a matrix:
  if(!is.matrix(CharacterTaxonMatrix)) stop("CharacterTaxonMatrix must be a matrix.")
  
  # Check taxon names are supplied:
  if(is.null(rownames(CharacterTaxonMatrix))) stop("CharacterTaxonMatrix must have rownames indicating taxa.")
  
  # Check taxon names are unique (could cause downstream issues if not):
  if(!ignore.duplicate.taxa) if(any(duplicated(rownames(CharacterTaxonMatrix)))) stop("Taxon names must be unique.")
  
  # Delete any column names (could cause downstream issues otherwise):
  if(!is.null(colnames(CharacterTaxonMatrix))) colnames(CharacterTaxonMatrix) <- NULL
  
  
  
  
  # WILL ONLY WORK FOR DISCRETE CHARACTERS!
  
  
  
  
  # List any mystery character types:
  mystery.characters <- setdiff(unique(unlist(strsplit(as.character(unique(as.vector(CharacterTaxonMatrix))), split = "&|/"))), c(as.character(0:31), NA))
  
  
  
  
  
  
  # If mystery character types are present warn user:
  if(length(mystery.characters) > 0) stop("Characters must either be the integers 0 to 31, NA for missing, & for polymorphisms, or / for uncertainties.")

  # Check supplied weights are correct length:
  if(!is.null(Weights) && length(Weights) != ncol(CharacterTaxonMatrix)) stop("Weights must have same length as number of characters in CharacterTaxonMatrix.")
  
  # Check supplied weights are numeric:
  if(!is.null(Weights) && !is.numeric(Weights)) stop("Weights must be numeric.")
  
  # Check supplied weights are non-negative:
  if(!is.null(Weights) && any(Weights < 0)) stop("Weights must not be negative.")
  
  # Check supplied ordering is the correct length:
  if(!is.null(ordering) && length(ordering) != ncol(CharacterTaxonMatrix)) stop("Ordering must have same length as number of characters in CharacterTaxonMatrix.")
  
  # Check ordering is of unord or ord type only:
  if(!is.null(ordering) && length(setdiff(ordering, c("unord", "ord"))) > 0) stop("Ordering must be unord or ord only.")
  
  # Check symbols are of correct length:
  if(!is.null(symbols) && length(symbols) >= (diff(range(as.numeric(unique(sort(unlist(strsplit(as.vector(CharacterTaxonMatrix), split = "&|/"))))))) + 1)) stop("Symbols must be at least as long as the range of character values in CharacterTaxonMatrix.")

  # Check symbols are single characters only:
  if(!is.null(symbols) && any(nchar(symbols) != 1)) stop("Symbols must be single characters only.")

  # Check header is a single value:
  if(length(header) != 1)  stop("Header text must be a single value.")
  
  # If no ordering set default to ordered.
  if(is.null(ordering)) ordering <- rep("ord", ncol(CharacterTaxonMatrix))

  # If no weights are set:
  if(is.null(Weights)) Weights <- rep(1, ncol(CharacterTaxonMatrix))

  # Calculate minimum values:
  min.vals <- apply(CharacterTaxonMatrix, 2, function(x) sort(as.numeric(unlist(strsplit(x, split = "&|/"))), decreasing = FALSE)[1])
  
  # Calculate maximum values:
  max.vals <- apply(CharacterTaxonMatrix, 2, function(x) sort(as.numeric(unlist(strsplit(x, split = "&|/"))), decreasing = TRUE)[1])
  
  # If any NAs in min.vals replace with zero:
  if(any(is.na(min.vals))) min.vals[is.na(min.vals)] <- 0
  
  # If any NAs in max.vals replace with zero:
  if(any(is.na(max.vals))) max.vals[is.na(max.vals)] <- 0

  # Default step matrices to NULL for now (may add this option in future):
  step.matrices <- NULL

  # If symbols were not set:
  if(is.null(symbols)) symbols <- c(c(0:9), LETTERS[1:22])[c(min(min.vals):max(max.vals)) + 1]

  # If wanting to equalise weights:
  if(equalise.weights) {
    
    # Get starting weights:
    Weights <- apply(rbind(c(max.vals - min.vals), rep(1, nchar)), 2, max)
    
    # Updaye weights for unordered characters:
    Weights[ordering == "unord"] <- 1
    
    # Update weights for ordered characters:
    Weights[ordering == "ord"] <- 1 / Weights[ordering == "ord"]
    
    # If there are step matrices (not technically using these yet, but I guess this will have to exist eventually):
    if(!is.null(step.matrices)) {
        
      # Get maximum distances for each step matrix:
      step.maxes <- unlist(lapply(lapply(step.matrices, as.numeric), max))
        
      # Update weights for step matrices:
      for(i in 1:length(step.maxes)) Weights[ordering == names(step.matrices)[i]] <- 1 / step.maxes[i]
        
    }
    
    # Ensure all weights are integers by multiplying by product of all reciprocals:
    Weights <- prod(unique(round(1 / Weights))) * Weights
    
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
    out <- sort(unlist(apply(matrix(unique(Weights)), 1, get.all.factors)))
    
    # As long as the maximum possible factor is greater than 1:
    while(max(rle(out)$values[rle(out)$lengths == length(unique(Weights))]) > 1) {
        
      # Divide through weights by largest common factor:
      Weights <- Weights / max(rle(out)$values[rle(out)$lengths == length(unique(Weights))])
        
      # Update factors for new weights:
      out <- sort(unlist(apply(matrix(unique(Weights)), 1, get.all.factors)))
        
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
  Matrix_1 <- list(NA, "STANDARD", CharacterTaxonMatrix, ordering, Weights, min.vals, max.vals, Characters)
  
  # Add names to Matrix_1:
  names(Matrix_1) <- c("BlockName", "Datatype", "Matrix", "Ordering", "Weights", "MinVals", "MaxVals", "Characters")

  # Assimilate into output:
  result <- list(Topper, Matrix_1)

  # Add names to output:
  names(result) <- c("Topper", "Matrix_1")

  # Return output:
  return(result)

}
