#' Writes out a morphological TNT data file
#' 
#' Writes out a morphological data file in Hennig86/TNT format.
#' 
#' Writes out a TNT (Goloboff et al. 2008) data file representing the distribution of discrete morphological characters in a set of taxa. Data must be in the format created by importing data with \link{ReadMorphNexus}.
#' 
#' Note that the function cannot yet deal with continuous characters or step matrices.
#' 
#' Currently all empty values (missing or inapplicable) are treated as missing and will be written to file as question marks.
#' 
#' @param clad.matrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#' @param filename The file name to write to. Should end in \code{.tnt}.
#' @param add.analysis.block Whether or not to add analysis block (i.e., tree search commands).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{WriteMorphNexus}
#'
#' @references Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. Cladistics, 24, 774-786.
#'
#' @keywords TNT
#'
#' @examples
#' 
#' # Write out Michaux 1989 to current working directory:
#' WriteMorphTNT(clad.matrix = Michaux1989, filename = "Michaux1989.tnt")
#' 
#' @export WriteMorphTNT
WriteMorphTNT <- function(clad.matrix, filename, add.analysis.block = FALSE) {

# NEED TO ADD WAY TO DEAL WITH MULTIPLE MATRICES (DISCRETE MORPHOLOGY, CONTINUOUS CHARACTERS, DNA).

    # If header text is present:
    if(nchar(paste(clad.matrix$header, collapse = "")) > 0) {
        
        # Create header line:
        headlines <- paste(paste("taxname=;\nmxram 4096;\n", paste("taxname +", max(nchar(rownames(clad.matrix$matrix))), ";\n", sep = ""), "nstates num 32;\nxread", sep = ""), paste(paste("'", clad.matrix$header, "'", sep = ""), collapse= " "), sep = "\n")
        
    # If no text
    } else {
        
        # Create header line:
        headlines <- paste("taxname=;\nmxram 4096;\n", paste("taxname +", max(nchar(rownames(clad.matrix$matrix))), ";\n", sep = ""), "nstates num 32;\nxread\n", sep = "")
        
    }
    
    # Make data block string:
    datablock <- paste(ncol(clad.matrix$matrix), nrow(clad.matrix$matrix))
    
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
        state.symbol <- clad.matrix$symbols[as.numeric(strsplit(i, "&")[[1]]) + 1]
        
        # If state symbol is a polymorphisms modify for writing:
        if(length(state.symbol) > 1) state.symbol <- paste("[", paste(state.symbol, collapse = ""), "]", sep = "")
        
        # Write correct state symbols into the matrix block:
        matrixblock[state.positions] <- state.symbol
        
    }
    
    # Make top of matrix block:
    matrixblock.top <- c("MATRIX", "")
    
    # Get spaces required after each taxon name:
    name.spaces <- (max(nchar(rownames(clad.matrix$matrix))) + 2) - nchar(rownames(clad.matrix$matrix))
    
    # Make matrix strings and store as rownames:
    for(i in 1:nrow(matrixblock)) rownames(matrixblock)[i] <- paste(rownames(matrixblock)[i], paste(rep(" ", name.spaces[i]), collapse = ""), paste(matrixblock[i, ], collapse = ""), sep = "")
    
    # Remake matrixblock as single string:
    matrixblock <- paste(paste(rownames(matrixblock), collapse = "\n"), "\n;\n", sep = "")

    # Get ordering (as plus = ordered or minus = unordered):
    ordering <- gsub("ord", "+", gsub("unord", "-", clad.matrix$ordering))

    # Set default for all characters as "on":
    onoff <- rep("[", ncol(clad.matrix$matrix))
    
    # If any characters are weighted zero then turn them "off":
    if(sum(clad.matrix$weights == 0) > 0) onoff[which(clad.matrix$weights == 0)] <- "]"

    # Set character weights:
    weights <- clad.matrix$weights
    
    # If there are zero weights convert these to 1 (as TNT does not allow them and they are dealt with by turning characters "off" instead):
    if(sum(weights == 0) > 0) weights[which(weights == 0)] <- 1

    # Round weights to two decimal places (required by TNT):
    weights <- round(weights, 2)

    # Remove leading zero from decimal weights (i.e., make 0.5 into .5) to fit TNT format:
    if(sum(weights < 1) > 0) weights[which(weights < 1)] <- gsub("0.", ".", as.character(weights[which(weights < 1)]))
    
    # Make orderin block (includes weights and on/off):
    ord.block <- paste("ccode", paste(paste(ordering, onoff, "/", weights, " ", c(0:(ncol(clad.matrix$matrix) - 1)), sep = ""), collapse = " "), ";\nproc/;\n", collapse = " ")
    
    # Compile output:
    out <- paste(paste(headlines, datablock, sep = "\n"), "\n", matrixblock, ord.block, sep = "")

    # If adding analysis block:
    if(add.analysis.block) {
        
        # If there are few enough taxa for an exact solution:
        if(nrow(clad.matrix$matrix) <= 24) {
            
            # Use implicit enumeration:
            anal <- c("collapse [;", "ienum;")
            
        # If there are too many taxa (requiring heuristics):
        } else {
            
            # First iteration with new tech (where scratch.tre is created):
            anal <- c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "mult 50 =tbr drift;", "tsave scratch.tre;", "save;", "tsave /;")
            
            # Iterations 2-20 (where trees are appended to scratch.tre):
            anal <- c(anal, rep(c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "mult 50 =tbr drift;", "tsave scratch.tre +;", "save;", "tsave /;"), 19))
            
            # Read in new technology trees (scratch.tre) and finish with a heuristic (tbr) search for all mpts:
            anal <- c(anal, c("hold 100000;", "shortread scratch.tre;", "bbreak=tbr;"))
            
        }
        
        # Get stripped file name for use in export lines:
        out.file <- strsplit(strsplit(filename, "/")[[1]][length(strsplit(filename, "/")[[1]])], "\\.")[[1]][1]
        
        # Make name for strict consensus and MPTs tree:
        strict.name <- paste("export -", out.file, "tntmpts_plus_strict.nex;", sep="")
        
        # Make MRP file name:
        mrp.name <- c("mrp;", paste("export ", out.file, "mrp.nex;", sep=""))
        
        # Add analysis block to output file:
        out <- gsub("proc/;\n", paste(paste(anal, collapse = "\n"), "nelsen*;", strict.name, paste(mrp.name, collapse = "\n"), "proc/;\n"), out)

    }


    # Write to file:
    write(out, filename)

}
