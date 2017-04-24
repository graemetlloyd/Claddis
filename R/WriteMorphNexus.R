#' Writes out a morphological #NEXUS data file
#' 
#' Writes out a morphological data file in #NEXUS format.
#' 
#' Writes out a #NEXUS (Maddison et al. 1997) data file representing the distribution of discrete morphological characters in a set of taxa. Data must be in the format created by importing data with \link{ReadMorphNexus}.
#' 
#' Note that the function cannot yet deal with continuous characters, but can write out step matrices.
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

# NEED TO ADD WAY TO DEAL WITH MULTIPLE MATRICES (DISCRETE MORPHOLOGY, CONTINUOUS CHARACTERS, DNA).

    # If header text is present:
    if(nchar(paste(clad.matrix$header, collapse = "")) > 0) {
        
        # Create header line:
        headlines <- paste("#NEXUS\n", paste(paste("[", clad.matrix$header, "]", sep = ""), collapse= "\n"), "\n", sep = "\n")
        
    # If no text
    } else {
        
        # Create header line:
        headlines <- "#NEXUS\n\n"
        
    }
    
    # Make data block string:
    datablock <- paste("BEGIN DATA;", paste("\t", "DIMENSIONS  NTAX=", nrow(clad.matrix$matrix), " NCHAR=", ncol(clad.matrix$matrix), ";", sep = ""), paste("\t", "FORMAT SYMBOLS= \" ", paste(clad.matrix$symbols, collapse = " "), "\" MISSING=? GAP=- ;", sep = ""), "\n", sep = "\n")
    
    # Create matrix block:
    matrixblock <- clad.matrix$matrix
    
    # Fill with question marks (missing state):
    matrixblock[1:length(matrixblock)] <- "?"
    
    # Get unique states:
    unique.states <- unique(sort(clad.matrix$matrix))
    
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
