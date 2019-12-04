#' Writes out a morphological TNT data file
#'
#' @description
#'
#' Writes out a morphological data file in Hennig86/TNT format.
#'
#' @param CladisticMatrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#' @param filename The file name to write to. Should end in \code{.tnt}.
#' @param add.analysis.block Whether or not to add analysis block (i.e., tree search commands).
#'
#' @details
#'
#' Writes out a TNT (Goloboff et al. 2008) data file representing the distribution of discrete morphological characters in a set of taxa. Data must be in the format created by importing data with \link{ReadMorphNexus}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{WriteMorphNexus}
#'
#' @references Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. Cladistics, 24, 774-786.
#'
#' @examples
#'
#' # Write out Michaux 1989 to current working directory:
#' WriteMorphTNT(CladisticMatrix = Michaux1989, filename = "Michaux1989.tnt")
#'
#' # Remove file when finished:
#' file.remove("Michaux1989.tnt")
#'
#' @export WriteMorphTNT
WriteMorphTNT <- function(CladisticMatrix, filename, add.analysis.block = FALSE) {
  
  # Subfunction to convert matrices back to symbols, missing and gap characters:
  MatrixConversion <- function(DataMatrix) {
    
    # If there are missing characters replace with missing symbol:
    if(any(is.na(DataMatrix$Matrix))) DataMatrix$Matrix[is.na(DataMatrix$Matrix)] <- "?"
    
    # If there are gap characters replace with gap symbol:
    if(sum(as.vector(DataMatrix$Matrix) == "") > 0) DataMatrix$Matrix[DataMatrix$Matrix == ""] <- "-"
    
    # If there are symbols (i.e., non-continuous data):
    if(length(DataMatrix$Characters$Symbols) > 0) {
      
      # If datatype is STANDARD set TNT symbols:
      if(DataMatrix$Datatype == "STANDARD") TNTSymbols <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V")
      
      # If datatype is non-STANDARD (but still discrete):
      if(DataMatrix$Datatype != "STANDARD") TNTSymbols <- DataMatrix$Characters$Symbols
      
      # In reverse order go through numbers:
      for(i in rev(TNTSymbols)) {
        
        # Replace current number with appropriate symbol:
        if(length(grep(as.character(which(TNTSymbols == i) - 1), DataMatrix$Matrix)) > 0) DataMatrix$Matrix <- gsub(as.character(which(TNTSymbols == i) - 1), i, DataMatrix$Matrix)
        
      }
      
    }
    
    # If there are uncertainties:
    if(length(grep("/", DataMatrix$Matrix)) > 0) {
      
      # Find cells that have uncertainties:
      Uncertainties <- grep("/", DataMatrix$Matrix)
      
      # Replace with all possible values in curly braces:
      DataMatrix$Matrix[Uncertainties] <- paste("[", unlist(lapply(strsplit(DataMatrix$Matrix[Uncertainties], split = "/"), paste, collapse = "")), "]", sep = "")
      
    }
    
    # If there are polymorphisms:
    if(length(grep("&", DataMatrix$Matrix)) > 0) {
      
      # Find cells with polymorphsims:
      Polymorphisms <- grep("&", DataMatrix$Matrix)
      
      # Resplae with values inside parentheses:
      DataMatrix$Matrix[Polymorphisms] <- paste("[", unlist(lapply(strsplit(DataMatrix$Matrix[Polymorphisms], split = "&"), paste, collapse = "")), "]", sep = "")
      
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
    
    # Return just the newly formatted matrix (now a vector):
    return(DataMatrix$Matrix)
    
  }
  
  # Isolate just data blocks (i.e., CladisticMatrix without Topper):
  DataBlocks <- CladisticMatrix[2:length(CladisticMatrix)]
  
  # Get block names:
  BlockNames <- unlist(lapply(DataBlocks, '[[', "BlockName"))
  
  # Get datatypes:
  Datatypes <- unlist(lapply(DataBlocks, '[[', "Datatype"))
  
  # NEXUS versions of datatypes:
  NEXUSversion <- c("CONTINUOUS", "DNA", "NUCLEOTIDE", "PROTEIN", "RESTRICTION", "RNA", "STANDARD")
  
  # TNT versions of datatypes:
  TNTversion <- c("continuous", "dna", "dna", "proteins", "numeric", "dna", "numeric")
  
  # Replace NEXUS versions of datatypes with TNT version (as best as I can guess!):
  for(i in 1:7) Datatypes <- gsub(NEXUSversion[i], TNTversion[i], Datatypes)
  
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
  
  # Set up header block (returns empty string if nothing there):
  HeaderBlock <- ifelse(length(CladisticMatrix$Topper$Header) > 0, paste("'", CladisticMatrix$Topper$Header, "'\n", sep = ""), "")
  
  # Set up character block (including MATRIX that will begin data):
  CharacterBlock <- paste("& [", Datatypes, "]\n", sep = "")
  
  # Take character block and meld with matri(ces) into matrix block(s):
  MatrixBlock <- paste(paste(CharacterBlock, unlist(lapply(DataBlocksAsTextStrings, paste, collapse = "\n")), "\n", sep = ""), collapse = "")
  
  # Get ordering of all characters in sequence:
  Ordering <- unname(unlist(lapply(DataBlocks, '[[', "Ordering")))
  
  # Get weights of all characters in sequence:
  Weights <- unname(unlist(lapply(DataBlocks, '[[', "Weights")))

  # Make sure step matrices are a list if null:
  if(!is.list(CladisticMatrix$Topper$StepMatrices)) CladisticMatrix$Topper$StepMatrices <- list(NULL)
  
  # If there are step matrices:
  if(any(!unlist(lapply(CladisticMatrix$Topper$StepMatrices, is.null)))) {
    
    # Empty vector to store hits (characters assigned to a step matrix):
    global_hits <- vector(mode = "numeric")
    
    # Empty vector to store all step matrix lines:
    all_step_matrix_lines <- vector(mode = "character")
    
    # Check there are not too many step matrices:
    if(length(CladisticMatrix$Topper$StepMatrices) > 32) stop("Too many (>32) step matrices for TNT!")
    
    # For each step matrix:
    for(i in 1:length(CladisticMatrix$Topper$StepMatrices)) {
      
      # Set up costs vector:
      costs <- vector(mode = "character")
      
      # For each row state (from):
      for(j in rownames(CladisticMatrix$Topper$StepMatrices[[i]])) {
        
        # For each column state (to):
        for(k in colnames(CladisticMatrix$Topper$StepMatrices[[i]])) {
          
          # Add cost of j to k transition to costs vector:
          costs <- c(costs, paste(j, ">", k, " ", CladisticMatrix$Topper$StepMatrices[[i]][j, k], sep = ""))
          
        }
        
      }
      
      # Format top of step matrix code:
      step_matrix_top <- paste("smatrix = ", i - 1, " (", names(CladisticMatrix$Topper$StepMatrices)[i], ")", sep = "")
      
      # Get hits (characters assigned to ith step matrix):
      hits <- which(Ordering == names(CladisticMatrix$Topper$StepMatrices)[i])
      
      # Add huts to global hits for all step matrices:
      global_hits <- c(global_hits, hits)
      
      # Stop if no hits!:
      if(length(hits) == 0) stop(paste("No characters assigned to step matrix: ", names(CladisticMatrix$Topper$StepMatrices)[i], ".", sep = ""))
      
      # Build step matrix lines:
      step_matrix_lines <- paste(c(step_matrix_top, costs, ";", paste("smatrix + ", names(CladisticMatrix$Topper$StepMatrices)[i], " ", paste(hits - 1, collapse = " "), ";", sep = "")), collapse = "\n")
      
      # Add step matrix lines to all step matrix lines:
      all_step_matrix_lines <- c(all_step_matrix_lines, step_matrix_lines)
      
    }
    
    # Make step matrix block:
    StepMatrixBlock <- paste(paste(c(all_step_matrix_lines, paste("ccode ( ", paste(sort(global_hits - 1), collapse = " "), ";", sep = "")), collapse = "\n"), "\n", sep = "")
    
  # If no step matri(ces):
  } else {
    
    # Create empty step matrix block:
    StepMatrixBlock <- ""
    
  }

  # Build ccode block:
  CCodeBlock <- paste("ccode ", paste(paste(ifelse(ifelse(Ordering == "cont", "ord", Ordering) == "ord", "+", "-"), ifelse(Weights == 0, "]", "["), "/", ifelse(Ordering == "cont", "1", ifelse(Weights == 0, 1, Weights)), " ", paste(1:sum(NCharacters) - 1), sep = ""), collapse = " "), " ;\n", sep = "")
  
  # Build full string with all blocks together:
  FullString <- paste("taxname=;\nmxram 4096;\ntaxname +", max(nchar(rownames(CladisticMatrix$Matrix_1$Matrix))), ";\nnstates num 32;\nxread\n", HeaderBlock, sum(NCharacters), " ", NTaxa, "\n", MatrixBlock, ";\n", CCodeBlock, StepMatrixBlock, "proc/;\n", sep = "")
  
  # If adding analysis block:
  if(add.analysis.block) {
    
    # If there are few enough taxa for an exact solution:
    if(NTaxa <= 24) {
      
      # Use implicit enumeration:
      AnalysisBlock <- c("collapse [;", "ienum;")
      
    # If there are too many taxa (requiring heuristics):
    } else {
      
      # First iteration with new tech (where scratch.tre is created):
      AnalysisBlock <- c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "tsave scratch.tre;", "save;", "tsave /;")
      
      # Iterations 2-20 (where trees are appended to scratch.tre):
      AnalysisBlock <- c(AnalysisBlock, rep(c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "tsave scratch.tre +;", "save;", "tsave /;"), 19))
      
      # Read in new technology trees (scratch.tre) and finish with a heuristic (tbr) search for all mpts:
      AnalysisBlock <- c(AnalysisBlock, c("hold 100000;", "shortread scratch.tre;", "bbreak=tbr;"))
      
    }
    
    # Get stripped file name for use in export lines:
    out.file <- strsplit(strsplit(filename, "/")[[1]][length(strsplit(filename, "/")[[1]])], "\\.")[[1]][1]
    
    # Make name for strict consensus and MPTs tree:
    strict.name <- paste("export -", out.file, "tntmpts_plus_strict.nex;", sep = "")
    
    # Make MRP file name:
    mrp.name <- c("mrp;", paste("export ", out.file, "mrp.nex;", sep = ""))
    
    # Add analysis block to output file:
    FullString <- gsub("proc/;\n", paste(paste(AnalysisBlock, collapse = "\n"), "nelsen*;", strict.name, paste(mrp.name, collapse = "\n"), "\nproc/;\n"), FullString)
    
  }
  
  # Write to file:
  write(FullString, filename)
  
}
