#' Writes out a morphological #NEXUS data file
#'
#' Writes out a morphological data file in #NEXUS format.
#'
#' Writes out a #NEXUS (Maddison et al. 1997) data file representing the distribution of characters in a set of taxa. Data must be in the format created by importing data with \link{ReadMorphNexus}.
#'
#' Currently all empty values (missing or inapplicable) are treated as missing and will be written to file as question marks.
#'
#' @param clad.matrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#' @param filename The file name to write to. Should end in \code{.nex}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{WriteMorphTNT}
#'
#' @references Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. Systematic Biology, 46, 590-621.
#'
#' @keywords NEXUS
#'
#' @examples
#'
#' # Write out Michaux 1989 to current working directory:
#' WriteMorphNexus(clad.matrix = Michaux1989, filename = "Michaux1989.nex")
#'
#' # Remove file when finished:
#' file.remove("Michaux1989.nex")
#'
#' @export WriteMorphNexus
WriteMorphNexus <- function(clad.matrix, filename) {
  
  # Subfunction to convert matrices back to symbols, missing and gap characters:
  MatrixConversion <- function(DataMatrix) {
    
    # If there are missing characters replace with missing symbol:
    if(any(is.na(DataMatrix$Matrix))) DataMatrix$Matrix[is.na(DataMatrix$Matrix)] <- DataMatrix$Characters$Missing
    
    # If there are gap characters replace with gap symbol:
    if(sum(as.vector(DataMatrix$Matrix) == "") > 0) DataMatrix$Matrix[DataMatrix$Matrix == ""] <- DataMatrix$Characters$Gap
    
    # If there are symbols (i.e., non-continuous data):
    if(length(DataMatrix$Characters$Symbols) > 0) {
      
      # In reverse order go through numbers:
      for(i in rev(DataMatrix$Characters$Symbols)) {
        
        # Replace current number with appropriate symbol:
        if(length(grep(as.character(which(DataMatrix$Characters$Symbols == i) - 1), DataMatrix$Matrix)) > 0) DataMatrix$Matrix <- gsub(as.character(which(DataMatrix$Characters$Symbols == i) - 1), i, DataMatrix$Matrix)
        
      }
      
    }
    
    # If there are uncertainties:
    if(length(grep("/", DataMatrix$Matrix)) > 0) {
      
      # Find cells that haev uncertainties:
      Uncertainties <- grep("/", DataMatrix$Matrix)
      
      # Repalce with all possible values in curly braces:
      DataMatrix$Matrix[Uncertainties] <- paste("{", unlist(lapply(strsplit(DataMatrix$Matrix[Uncertainties], split = "/"), paste, collapse = "")), "}", sep = "")
      
    }
    
    # If there are polymorphisms:
    if(length(grep("&", DataMatrix$Matrix)) > 0) {
      
      # Find cells with polymorphsims:
      Polymorphisms <- grep("&", DataMatrix$Matrix)
      
      # Repale with values inside parentheses:
      DataMatrix$Matrix[Polymorphisms] <- paste("(", unlist(lapply(strsplit(DataMatrix$Matrix[Polymorphisms], split = "&"), paste, collapse = "")), ")", sep = "")
      
    }
    
    # Get equal length taxon names (with added spaces):
    TaxonNamesWithTrailingSpaces <- paste(rownames(DataMatrix$Matrix), unlist(lapply(lapply(as.list((max(nchar(rownames(DataMatrix$Matrix))) + 2) - nchar(rownames(DataMatrix$Matrix))), rep, x = " "), paste, collapse = "")), sep = "")
    
    # If block is continuous:
    if(DataMatrix$Datatype == "CONTINUOUS") {
      
      # Format rows with spaces between values:
      DataMatrix$Matrix <- paste(TaxonNamesWithTrailingSpaces, (apply(DataMatrix$Matrix, 1, paste, collapse = " ")), sep = "")
      
    # If block is non-continuous:
    } else {
      
      # Format rows without spaces between values:
      DataMatrix$Matrix <- paste(TaxonNamesWithTrailingSpaces, (apply(DataMatrix$Matrix, 1, paste, collapse = "")), sep = "")
      
    }
    
    # Return just the newly fromateed matrix (now a vector):
    return(DataMatrix$Matrix)
    
  }
  
  # If header text is present:
  if(nchar(paste(clad.matrix$Topper$Header, collapse = "")) > 0) {
    
    # Create header line:
    headlines <- paste("#NEXUS\n", paste(paste("[", clad.matrix$Topper$Header, "]", sep = ""), collapse= "\n"), "\n", sep = "\n")
    
  # If no text
  } else {
    
    # Create header line:
    headlines <- "#NEXUS\n\n"
    
  }
  
  # Isolate just data blocks (i.e., clad.matrix without Topper):
  DataBlocks <- clad.matrix[2:length(clad.matrix)]
  
  # Get block names:
  BlockNames <- unlist(lapply(DataBlocks, '[[', "BlockName"))
  
  # Get datatypes:
  Datatypes <- unlist(lapply(DataBlocks, '[[', "Datatype"))
  
  # Get number of taxa:
  NTaxa <- unname(unlist(lapply(lapply(DataBlocks, '[[', "Matrix"), nrow))[1])
  
  # Get number of characters:
  NCharacters <- unlist(lapply(lapply(DataBlocks, '[[', "Matrix"), ncol))
  
  # Get symbols strings:
  Symbols <- unlist(lapply(lapply(lapply(DataBlocks, '[[', "Characters"), '[[', "Symbols"), paste, collapse = " "))
  
  # Get missing value:
  Missing <- unlist(lapply(lapply(DataBlocks, '[[', "Characters"), '[[', "Missing"))
  
  # Get gap symbol:
  Gap <- unlist(lapply(lapply(DataBlocks, '[[', "Characters"), '[[', "Gap"))
  
  # Conver matrices to vectors of text strings:
  DataBlocksAsTextStrings <- lapply(DataBlocks, MatrixConversion)
  
  
  
  
  
  
  
  
  # Stet satrt of data block (data for single block, characters if multiple):
  DataBlockStart <- ifelse(length(DataBlocks) > 1, "BEGIN CHARACTERS;\n", paste("BEGIN DATA;\n\tDIMENSIONS  NTAX=", NTaxa, " ", collapse = ""))
  
  # Reformta symbols (returning empty string for continuous) eady for data or character block(s):
  Symbols <- unlist(lapply(as.list(Symbols), function(x) ifelse(nchar(x) > 0, paste(" SYMBOLS=\" ", x, "\" ", sep = ""), " ")))
  
  # Set up character and format lines:
  CharacterAndFormatLines <- paste("NCHAR=", NCharacters, ";\n\tFORMAT ", "DATATYPE=", Datatypes, Symbols, "MISSING=", Missing, " GAP=", Gap, " ;", sep = "")






  # Set up taxa block (only required if multiple matrix blocks as sets number of taxa):
  TaxaBlock <- ifelse(length(DataBlocks) > 1, paste("\tBEGIN TAXA;\n\tDIMENSIONS NTAX=", NTaxa, ";\n\tTAXLABELS\n\t\t", paste(rownames(clad.matrix$Matrix_1$Matrix), collapse = " "), "\n\t;\n\tEND;", sep = ""), "")
  
  # Make data block string:
  datablock <- paste(DataBlockStart, paste("NCHAR=", NCharacters, ";", sep = ""), paste("\t", "FORMAT SYMBOLS= \" ", Symbols, "\" MISSING=", Missing, " GAP=", Gap, " ;", sep = ""), "\n", sep = "\n")
  
  
  
  
  
    
  # TITLE if blocks have titles, check that single matrix with title can store this maybe in ReadMorphNexus?
  # Removes symbols for continuous
  # Remove ntax for multiple blocks
  # ADD TAXA BLOCK FOR MULTIPLE BLOCK CASES
  
  
  
  
  
  
  # For each state (including polymorphisms):
  for(i in unique.states) {
    
    # Store positions in matrix of ith state:
    state.positions <- which(clad.matrix$matrix == i)
    
    # Get state symbol:
    ifelse(length(grep("&", i)) > 0, state.symbol <- clad.matrix$symbols[as.numeric(strsplit(i, "&")[[1]]) + 1], state.symbol <- clad.matrix$symbols[as.numeric(i) + 1])
    
    # If state symbol is a polymorphisms modify for writing:
    if(length(state.symbol) > 1) state.symbol <- paste("(", paste(state.symbol, collapse = ""), ")", sep = "")
    
    # Write correct state symbols into the matrix block:
    matrixblock[state.positions] <- state.symbol
    
  }
  
  # Make top of matrix block:
  matrixblock.top <- c("MATRIX", "")
  
  # Get spaces required after each taxon name:
  name.spaces <- (max(nchar(rownames(clad.matrix$matrix)))+2)-nchar(rownames(clad.matrix$matrix))
  
  # Make matrix strings and store as rownames:
  for(i in 1:nrow(matrixblock)) rownames(matrixblock)[i] <- paste(rownames(matrixblock)[i], paste(rep(" ", name.spaces[i]), collapse = ""), paste(matrixblock[i, ], collapse = ""), sep = "")
  
  # Remake matrixblock as single string:
  matrixblock <- paste("MATRIX\n\n", paste(rownames(matrixblock), collapse = "\n"), "\n;\nEND;\n\n", sep = "")
  
  # Case if all characters are unordered:
  if(length(unique(clad.matrix$ordering)) == 1 && unique(clad.matrix$ordering) == "unord") ord.block <- "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;"
  
  # Case if all characters are ordered:
  if(length(unique(clad.matrix$ordering)) == 1 && unique(clad.matrix$ordering) == "ord") ord.block <- "\tOPTIONS  DEFTYPE=ord PolyTcount=MINSTEPS ;"
  
  # Case if mix of character types:
  if(length(unique(clad.matrix$ordering)) > 1) {
    
    # Empty vector to store character strings:
    character.strings <- vector(mode = "numeric")
    
    # Get character types (ord, unord, stepmatrix):
    character.types <- unique(clad.matrix$ordering)
    
    # For each character type (unord, ord, stepmatrix):
    for(i in character.types) {
      
      # Find characters of current character type:
      character.numbers <- which(clad.matrix$ordering == i)
      
      # If there are contiguous character numbers:
      if(any(diff(character.numbers) == 1)) {
        
        # Set up character string:
        char.string <- paste(i, ":", sep = "")
        
        # Set first character number:
        j <- character.numbers[1]
        
        # While we havem't reached the end:
        while(j < character.numbers[length(character.numbers)]) {
          
          # What is last character in the contiguous run:
          endrun <- character.numbers[which(character.numbers == j):length(character.numbers)][min(c(which(diff(character.numbers[which(character.numbers == j):length(character.numbers)]) > 1), length(character.numbers[which(character.numbers == j):length(character.numbers)])))]
          
          # If last character in run is after first character:
          if(endrun > j) {
            
            # Store run by collapsing with a dash:
            char.string <- c(char.string, paste(j, "-", endrun, sep = ""))
            
            # Update j to the next character (or keep as last if last character):
            j <- character.numbers[min(c(which(character.numbers == endrun) + 1, length(character.numbers)))]
            
            # If the last character is also the current character (there is no contiguous run):
          } else {
            
            # Store as single character:
            char.string <- c(char.string, j)
            
            # Update j to next character in sequence (or keep as last character if last character):
            j <- character.numbers[min(c(which(character.numbers == j) + 1, length(character.numbers)))]
            
            # Special case where last character is isolated (not part of a contiguous run):
            if(j == character.numbers[length(character.numbers)]) char.string <- c(char.string, j)
            
          }
          
        }
        
        # Complete character strong:
        char.string <- paste(char.string, collapse = " ")
        
        # If character numbers are all non-contiguous:
      } else {
        
        # Complete character string:
        char.string <- paste(paste(i, ":", sep = ""), paste(character.numbers, collapse = " "), collapse = " ")
        
      }
      
      # Combine character strings:
      character.strings <- c(character.strings, char.string)
      
    }
    
    # Create ordering block:
    ord.block <- paste("\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n\tTYPESET * UNTITLED  = ", paste(sort(character.strings, decreasing = TRUE), collapse = ", "), ";", sep = "")
    
  }
  
  # If different weights are specified:
  if(length(unique(clad.matrix$weights)) > 1) {
    
    # Empty vector to store weight strings:
    weight.strings <- vector(mode = "numeric")
    
    # Get character types (ord, unord, stepmatrix):
    weight.values <- unique(clad.matrix$weights)
    
    # For each weight value:
    for(i in weight.values) {
      
      # Find characters of current weight:
      character.numbers <- which(clad.matrix$weights == i)
      
      # If there are contiguous character numbers:
      if(any(diff(character.numbers) == 1)) {
        
        # Set up weight string:
        weight.string <- paste(i, ":", sep = "")
        
        # Set first character number:
        j <- character.numbers[1]
        
        # While we havem't reached the end:
        while(j < character.numbers[length(character.numbers)]) {
          
          # What is last character in the contiguous run:
          endrun <- character.numbers[which(character.numbers == j):length(character.numbers)][min(c(which(diff(character.numbers[which(character.numbers == j):length(character.numbers)]) > 1), length(character.numbers[which(character.numbers == j):length(character.numbers)])))]
          
          # If last character in run is after first character:
          if(endrun > j) {
            
            # Store run by collapsing with a dash:
            weight.string <- c(weight.string, paste(j, "-", endrun, sep = ""))
            
            # Update j to the next character (or keep as last if last character):
            j <- character.numbers[min(c(which(character.numbers == endrun) + 1, length(character.numbers)))]
            
            # If the last character is also the current character (there is no contiguous run):
          } else {
            
            # Store as single character:
            weight.string <- c(weight.string, j)
            
            # Update j to next character in sequence (or keep as last character if last character):
            j <- character.numbers[min(c(which(character.numbers == j) + 1, length(character.numbers)))]
            
            # Special case where last character is isolated (not part of a contiguous run):
            if(j == character.numbers[length(character.numbers)]) weight.string <- c(weight.string, j)
            
          }
          
        }
        
        # Complete character strong:
        weight.string <- paste(weight.string, collapse = " ")
        
        # If character numbers are all non-contiguous:
      } else {
        
        # Complete character string:
        weight.string <- paste(paste(i, ":", sep = ""), paste(character.numbers, collapse = " "), collapse = " ")
        
      }
      
      # Combine character strings:
      weight.strings <- c(weight.strings, weight.string)
      
    }
    
    # Create ordering block:
    weight.block <- paste("\tWTSET * UNTITLED  = ", paste(sort(weight.strings, decreasing = TRUE), collapse = ", "), ";", sep = "")
    
    # Add weights to ord block
    ord.block <- paste(ord.block, "\n", weight.block, sep = "")
    
  }
  
  # If there are step matrices:
  if(!is.null(clad.matrix$step.matrices)) {
    
    # Create empty step matrices vector:
    step_matrices <- vector(mode = "character")
    
    # For each step matrix:
    for(i in 1:length(clad.matrix$step.matrices)) {
      
      # Replace diagonal with a period:
      diag(clad.matrix$step.matrices[[i]]) <- "."
      
      # Add ith step matrix to step matrices:
      step_matrices[i] <- paste(c(paste("\tUSERTYPE '", names(clad.matrix$step.matrices)[i], "' (STEPMATRIX) = ", ncol(clad.matrix$step.matrices[[i]]), sep = ""), paste("\t", paste(colnames(clad.matrix$step.matrices[[i]]), collapse = " "), sep = ""), paste("\t", apply(clad.matrix$step.matrices[[i]], 1, paste, collapse = " "), sep = ""), "\t;\n"), collapse = "\n")
      
    }
    
  }
  
  # If there are no step matrices:
  if(is.null(clad.matrix$step.matrices)) {
    
    # Finish assumptions block:
    ass.block <- paste("BEGIN ASSUMPTIONS;\n", ord.block, "\nEND;\n", sep = "")
    
    # If there are step matrices:
  } else {
    
    # Finish assumptions block:
    ass.block <- paste("BEGIN ASSUMPTIONS;\n", paste(step_matrices, collapse = ""), ord.block, "\nEND;\n", collapse = "")
    
  }
  
  # Compile output:
  out <- paste(headlines, datablock, matrixblock, ass.block, sep = "")
  
  # Write to file:
  write(out, filename)
  
}
