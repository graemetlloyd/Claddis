#' Writes out a morphological #NEXUS data file
#'
#' @description
#'
#' Writes out a morphological data file in #NEXUS format.
#'
#' @param cladistic_matrix The cladistic matrix in the format imported by \link{read_nexus_matrix}.
#' @param file_name The file name to write to. Should end in \code{.nex}.
#'
#' @details
#'
#' Writes out a #NEXUS (Maddison et al. 1997) data file representing the distribution of characters in a set of taxa. Data must be in the format created by importing data with \link{read_nexus_matrix}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{write_tnt_matrix}
#'
#' @references
#'
#' Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. \emph{Systematic Biology}, \bold{46}, 590-621.
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{compactify_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Write out Michaux 1989 to current working directory:
#' write_nexus_matrix(cladistic_matrix = michaux_1989, file_name = "michaux_1989.nex")
#'
#' # Remove file when finished:
#' file.remove("michaux_1989.nex")
#'
#' @export write_nexus_matrix
write_nexus_matrix <- function(cladistic_matrix, file_name) {
  
  # Subfunction to convert matrices back to symbols, missing and gap characters:
  convert_matrix <- function(DataMatrix) {
    
    # If there are missing characters replace with missing symbol:
    if (any(is.na(DataMatrix$matrix))) DataMatrix$matrix[is.na(DataMatrix$matrix)] <- DataMatrix$characters$missing
    
    # If there are gap characters replace with gap symbol:
    if (sum(as.vector(DataMatrix$matrix) == "") > 0) DataMatrix$matrix[DataMatrix$matrix == ""] <- DataMatrix$characters$gap
    
    # If there are symbols (i.e., non-continuous data):
    if (length(DataMatrix$characters$symbols) > 0) {
      
      # In reverse order go through numbers:
      for(i in rev(DataMatrix$characters$symbols)) {
        
        # Replace current number with appropriate symbol:
        if (length(grep(as.character(which(DataMatrix$characters$symbols == i) - 1), DataMatrix$matrix)) > 0) DataMatrix$matrix <- gsub(as.character(which(DataMatrix$characters$symbols == i) - 1), i, DataMatrix$matrix)
        
      }
      
    }
    
    # If there are uncertainties:
    if (length(grep("/", DataMatrix$matrix)) > 0) {
      
      # Find cells that haev uncertainties:
      Uncertainties <- grep("/", DataMatrix$matrix)
      
      # Repalce with all possible values in curly braces:
      DataMatrix$matrix[Uncertainties] <- paste("{", unlist(lapply(strsplit(DataMatrix$matrix[Uncertainties], split = "/"), paste, collapse = "")), "}", sep = "")
      
    }
    
    # If there are polymorphisms:
    if (length(grep("&", DataMatrix$matrix)) > 0) {
      
      # Find cells with polymorphsims:
      Polymorphisms <- grep("&", DataMatrix$matrix)
      
      # Repale with values inside parentheses:
      DataMatrix$matrix[Polymorphisms] <- paste("(", unlist(lapply(strsplit(DataMatrix$matrix[Polymorphisms], split = "&"), paste, collapse = "")), ")", sep = "")
      
    }
    
    # Get equal length taxon names (with added spaces):
    TaxonNamesWithTrailingSpaces <- paste(rownames(DataMatrix$matrix), unlist(lapply(lapply(as.list((max(nchar(rownames(DataMatrix$matrix))) + 2) - nchar(rownames(DataMatrix$matrix))), rep, x = " "), paste, collapse = "")), sep = "")
    
    # If block is continuous:
    if (DataMatrix$datatype == "CONTINUOUS") {
      
      # Format rows with spaces between values:
      DataMatrix$matrix <- paste(TaxonNamesWithTrailingSpaces, (apply(DataMatrix$matrix, 1, paste, collapse = " ")), sep = "")
      
    # If block is non-continuous:
    } else {
      
      # Format rows without spaces between values:
      DataMatrix$matrix <- paste(TaxonNamesWithTrailingSpaces, (apply(DataMatrix$matrix, 1, paste, collapse = "")), sep = "")
      
    }
    
    # Return just the newly formatted matrix (now a vector):
    return(DataMatrix$matrix)
    
  }
  
  # Subfunction for collapsing strings of character numbers (i.e., for weight and ordering in assumptions block):
  zip_string <- function(x) {
    
    # Set up empty zipped version of strin:
    zipped_string <- vector(mode = "character")
    
    # As long as there are zip-able parts to x remaining:
    while(any(diff(x) == 1)) {
      
      # Find start of first zippable part:
      zip_start <- which(diff(x) == 1)[1]
      
      # Find length of zippable part:
      zip_length <- rle(diff(x))$lengths[which(rle(diff(x))$values == 1)[1]]
      
      # Use length and start to find end of zip:
      zip_end <- zip_start + zip_length
      
      # Add zippable part to string (and any preceding unzippable parts):
      zipped_string <- c(zipped_string, setdiff(x[1:zip_end], x[zip_start:zip_end]), paste(x[zip_start], x[zip_end], sep = "-"))
      
      # Cut x down ready to find next zippable part:
      x <- x[-(1:zip_end)]
      
    }
    
    # If there are characters left (or there are no zippable parts):
    if (length(x) > 0) zipped_string <- c(zipped_string, x)
    
    # Collapse to a single string with spaces:
    zipped_string <- paste(zipped_string, collapse = " ")
    
    # Return zipped string:
    return(zipped_string)
    
  }
  
  # Isolate just data blocks (i.e., cladistic_matrix without topper):
  DataBlocks <- cladistic_matrix[2:length(cladistic_matrix)]
  
  # Get block names:
  block_names <- unlist(lapply(DataBlocks, '[[', "block_name"))
  
  # Get datatypes:
  datatypes <- unlist(lapply(DataBlocks, '[[', "datatype"))
  
  # Get number of taxa:
  NTaxa <- unname(unlist(lapply(lapply(DataBlocks, '[[', "matrix"), nrow))[1])
  
  # Get number of characters:
  Ncharacters <- unlist(lapply(lapply(DataBlocks, '[[', "matrix"), ncol))
  
  # Get symbols strings:
  symbols <- unlist(lapply(lapply(lapply(DataBlocks, '[[', "characters"), '[[', "symbols"), paste, collapse = " "))
  
  # Get missing value:
  missing <- unlist(lapply(lapply(DataBlocks, '[[', "characters"), '[[', "missing"))
  
  # Get gap symbol:
  gap <- unlist(lapply(lapply(DataBlocks, '[[', "characters"), '[[', "gap"))
  
  # Conver matrices to vectors of text strings:
  DataBlocksAsTextStrings <- lapply(DataBlocks, convert_matrix)
  
  # Set up header block (returns empty string if nothing there):
  headerBlock <- ifelse(nchar(cladistic_matrix$topper$header) > 0, paste("[", cladistic_matrix$topper$header, "]\n\n", sep = ""), "")

  # Set up taxa block (only required if multiple matrix blocks as sets number of taxa, will be empty string otherwise):
  TaxaBlock <- ifelse(length(DataBlocks) > 1, paste("BEGIN TAXA;\n\tDIMENSIONS NTAX=", NTaxa, ";\n\tTAXLABELS\n\t\t", paste(rownames(cladistic_matrix$matrix_1$matrix), collapse = " "), "\n;\nEND;\n\n", sep = ""), "")
  
  # Set up data block (only required if a single matrix block):
  DataBlock <- ifelse(length(DataBlocks) == 1, paste("BEGIN DATA;\n\tDIMENSIONS  NTAX=", NTaxa, " NCHAR=", Ncharacters, " ;\n\tFORMAT DATATYPE=", datatypes, " SYMBOLS=\" ", symbols, "\" MISSING=", missing, " GAP=", gap, " ;\n", sep = ""), "")
  
  # Set up character block (including MATRIX that will begin data):
  CharacterBlock <- ifelse(rep(length(DataBlocks), length(DataBlocks)) > 1, paste("BEGIN CHARACTERS;\n\t", ifelse(nchar(block_names) > 0, paste("TITLE  ", block_names, ";\n", sep = ""), ""), "\tDIMENSIONS  NCHAR=", Ncharacters, ";\n\tFORMAT  DATATYPE=", ifelse(datatypes == "CONTINUOUS", "CONTINUOUS ", paste(datatypes, " SYMBOLS=\" ", symbols, "\" ", sep = "")), "MISSING=", missing, " GAP=", gap, " ;\nMATRIX\n\n", sep = ""), "MATRIX\n\n")
  
  # Take character block and meld with matri(ces) into matrix block(s):
  MatrixBlock <- paste(paste(CharacterBlock, unlist(lapply(DataBlocksAsTextStrings, paste, collapse = "\n")), "\n;\nEND;\n\n", sep = ""), collapse = "")
  
  # Make sure step matrices are a list if null:
  if (!is.list(cladistic_matrix$topper$step_matrices)) cladistic_matrix$topper$step_matrices <- list(NULL)
  
  # Create step matrix block:
  StepMatrixBlock <- paste(ifelse(!unlist(lapply(cladistic_matrix$topper$step_matrices, is.null)), paste(paste("\tUSERTYPE '", names(cladistic_matrix$topper$step_matrices), "' (STEPMATRIX) = ", unlist(lapply(cladistic_matrix$topper$step_matrices, ncol)), "\n", sep = ""), paste("\t", unlist(lapply(lapply(cladistic_matrix$topper$step_matrices, colnames), paste, collapse = " ")), "\n\t", sep = ""), unlist(lapply(lapply(lapply(cladistic_matrix$topper$step_matrices, function(x) { diag(x) <- "."; return(x) }), apply, 1, paste, collapse = " "), paste, collapse = "\n\t")), "\n\t;\n", sep = ""), ""), collapse = "")

  # Get ordering of all characters in sequence:
  ordering <- unlist(lapply(DataBlocks, '[[', "ordering"))
  
  # Get weights of all characters in sequence:
  weights <- unlist(lapply(DataBlocks, '[[', "weights"))

  # Create options block (if no block names):
  if (all(is.na(block_names))) OptionsBlock <- paste(ifelse(all(ordering == "unord"), "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n", ifelse(all(ordering == "ord"), "\tOPTIONS  DEFTYPE=ord PolyTcount=MINSTEPS ;\n", "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n")), ifelse(length(unique(ordering)) == 1 && length(setdiff(unique(ordering), c("ord", "unord"))) == 0, "", paste("\tTYPESET * UNTITLED  = ", paste(paste(sort(unique(ordering)), unlist(lapply(lapply(lapply(as.list(sort(unique(ordering))), '==', ordering), which), zip_string)), sep = ": "), collapse = ", "), ";\n", sep = "")), collapse = "")
  
  # Create options block (if there are block names):
  if (!all(is.na(unlist(block_names)))) OptionsBlock <- paste(paste("\tTYPESET * UNTITLED  (CHARACTERS = ", block_names, ")  =  ", unlist(lapply(lapply(DataBlocks, '[[', "ordering"), function(x) paste(paste(paste(sort(unique(x)), unlist(lapply(lapply(lapply(as.list(sort(unique(x))), '==', x), which), zip_string)), sep = ": "), collapse = ", "), sep = ""))), ";\n", sep = ""), collapse = "")
  
  # Replace cont with Squared if continuous characters present:
  if (length(grep(" cont: ", OptionsBlock)) > 0) OptionsBlock <- gsub(" cont: ", " Squared: ", OptionsBlock)
  
  # Convert continuosu character weights to one before making weights block:
  weights[ordering == "cont"] <- 1
  
  # Create weights block (if no block names):
  if (all(is.na(block_names))) weightsBlock <- ifelse(all(weights == 1), "", paste("\tWTSET * UNTITLED  = ", paste(paste(sort(unique(weights)), unlist(lapply(lapply(lapply(as.list(sort(unique(weights))), '==', weights), which), zip_string)), sep = ": "), collapse = ", "), ";\n", sep = ""))
  
  # Create weights block (if there are block names):
  if (!all(is.na(unlist(block_names)))) weightsBlock <- paste(paste("\tWTSET * UNTITLED  (CHARACTERS = ", block_names, ")  =  ", unlist(lapply(lapply(DataBlocks, '[[', "weights"), function(x) paste(paste(paste(sort(unique(x)), unlist(lapply(lapply(lapply(as.list(sort(unique(x))), '==', x), which), zip_string)), sep = ": "), collapse = ", "), sep = ""))), ";\n", sep = ""), collapse = "")
  
  # Build assumptions block:
  AssumptionBlock <- paste("BEGIN ASSUMPTIONS;\n", StepMatrixBlock, OptionsBlock, weightsBlock, "END;\n", sep = "")
  
  # Build full string with all blocks together:
  FullString <- paste("#NEXUS\n\n", headerBlock, TaxaBlock, DataBlock, MatrixBlock, AssumptionBlock, sep = "")

  # Write to file:
  write(FullString, file_name)

}
