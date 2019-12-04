#' Writes out a morphological #NEXUS data file
#'
#' @description
#'
#' Writes out a morphological data file in #NEXUS format.
#'
#' @param CladisticMatrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#' @param filename The file name to write to. Should end in \code{.nex}.
#'
#' @details
#'
#' Writes out a #NEXUS (Maddison et al. 1997) data file representing the distribution of characters in a set of taxa. Data must be in the format created by importing data with \link{ReadMorphNexus}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{WriteMorphTNT}
#'
#' @references Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. Systematic Biology, 46, 590-621.
#'
#' @examples
#'
#' # Write out Michaux 1989 to current working directory:
#' WriteMorphNexus(CladisticMatrix = Michaux1989, filename = "Michaux1989.nex")
#'
#' # Remove file when finished:
#' file.remove("Michaux1989.nex")
#'
#' @export WriteMorphNexus
WriteMorphNexus <- function(CladisticMatrix, filename) {
  
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
    
    # Return just the newly formatted matrix (now a vector):
    return(DataMatrix$Matrix)
    
  }
  
  # Subfunction for collapsing strings of character numbers (i.e., for weight and ordering in assumptions block):
  Zipper <- function(x) {
    
    # Set up empty zipped version of strin:
    ZippedString <- vector(mode = "character")
    
    # As long as there are zip-able parts to x remaining:
    while(any(diff(x) == 1)) {
      
      # Find start of first zippable part:
      StartOfZip <- which(diff(x) == 1)[1]
      
      # Find length of zippable part:
      ZipLength <- rle(diff(x))$lengths[which(rle(diff(x))$values == 1)[1]]
      
      # Use length and start to find end of zip:
      EndOfZip <- StartOfZip + ZipLength
      
      # Add zippable part to string (and any preceding unzippable parts):
      ZippedString <- c(ZippedString, setdiff(x[1:EndOfZip], x[StartOfZip:EndOfZip]), paste(x[StartOfZip], x[EndOfZip], sep = "-"))
      
      # Cut x down ready to find next zippable part:
      x <- x[-(1:EndOfZip)]
      
    }
    
    # If there are characters left (or there are no zippable parts):
    if(length(x) > 0) ZippedString <- c(ZippedString, x)
    
    # Collapse to a single string with spaces:
    ZippedString <- paste(ZippedString, collapse = " ")
    
    # Return zipped string:
    return(ZippedString)
    
  }
  
  # Isolate just data blocks (i.e., CladisticMatrix without Topper):
  DataBlocks <- CladisticMatrix[2:length(CladisticMatrix)]
  
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
  
  # Set up header block (returns empty string if nothing there):
  HeaderBlock <- ifelse(nchar(CladisticMatrix$Topper$Header) > 0, paste("[", CladisticMatrix$Topper$Header, "]\n\n", sep = ""), "")

  # Set up taxa block (only required if multiple matrix blocks as sets number of taxa, will be empty string otherwise):
  TaxaBlock <- ifelse(length(DataBlocks) > 1, paste("BEGIN TAXA;\n\tDIMENSIONS NTAX=", NTaxa, ";\n\tTAXLABELS\n\t\t", paste(rownames(CladisticMatrix$Matrix_1$Matrix), collapse = " "), "\n;\nEND;\n\n", sep = ""), "")
  
  # Set up data block (only required if a single matrix block):
  DataBlock <- ifelse(length(DataBlocks) == 1, paste("BEGIN DATA;\n\tDIMENSIONS  NTAX=", NTaxa, " NCHAR=", NCharacters, " ;\n\tFORMAT DATATYPE=", Datatypes, " SYMBOLS=\" ", Symbols, "\" MISSING=", Missing, " GAP=", Gap, " ;\n", sep = ""), "")
  
  # Set up character block (including MATRIX that will begin data):
  CharacterBlock <- ifelse(rep(length(DataBlocks), length(DataBlocks)) > 1, paste("BEGIN CHARACTERS;\n\t", ifelse(nchar(BlockNames) > 0, paste("TITLE  ", BlockNames, ";\n", sep = ""), ""), "\tDIMENSIONS  NCHAR=", NCharacters, ";\n\tFORMAT  DATATYPE=", ifelse(Datatypes == "CONTINUOUS", "CONTINUOUS ", paste(Datatypes, " SYMBOLS=\" ", Symbols, "\" ", sep = "")), "MISSING=", Missing, " GAP=", Gap, " ;\nMATRIX\n\n", sep = ""), "MATRIX\n\n")
  
  # Take character block and meld with matri(ces) into matrix block(s):
  MatrixBlock <- paste(paste(CharacterBlock, unlist(lapply(DataBlocksAsTextStrings, paste, collapse = "\n")), "\n;\nEND;\n\n", sep = ""), collapse = "")
  
  # Make sure step matrices are a list if null:
  if(!is.list(CladisticMatrix$Topper$StepMatrices)) CladisticMatrix$Topper$StepMatrices <- list(NULL)
  
  # Create step matrix block:
  StepMatrixBlock <- paste(ifelse(!unlist(lapply(CladisticMatrix$Topper$StepMatrices, is.null)), paste(paste("\tUSERTYPE '", names(CladisticMatrix$Topper$StepMatrices), "' (STEPMATRIX) = ", unlist(lapply(CladisticMatrix$Topper$StepMatrices, ncol)), "\n", sep = ""), paste("\t", unlist(lapply(lapply(CladisticMatrix$Topper$StepMatrices, colnames), paste, collapse = " ")), "\n\t", sep = ""), unlist(lapply(lapply(lapply(CladisticMatrix$Topper$StepMatrices, function(x) { diag(x) <- "."; return(x) }), apply, 1, paste, collapse = " "), paste, collapse = "\n\t")), "\n\t;\n", sep = ""), ""), collapse = "")

  # Get ordering of all characters in sequence:
  Ordering <- unlist(lapply(DataBlocks, '[[', "Ordering"))
  
  # Get weights of all characters in sequence:
  Weights <- unlist(lapply(DataBlocks, '[[', "Weights"))

  # Create options block (if no block names):
  if(all(is.na(BlockNames))) OptionsBlock <- paste(ifelse(all(Ordering == "unord"), "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n", ifelse(all(Ordering == "ord"), "\tOPTIONS  DEFTYPE=ord PolyTcount=MINSTEPS ;\n", "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n")), ifelse(length(unique(Ordering)) == 1 && length(setdiff(unique(Ordering), c("ord", "unord"))) == 0, "", paste("\tTYPESET * UNTITLED  = ", paste(paste(sort(unique(Ordering)), unlist(lapply(lapply(lapply(as.list(sort(unique(Ordering))), '==', Ordering), which), Zipper)), sep = ": "), collapse = ", "), ";\n", sep = "")), collapse = "")
  
  # Create options block (if there are block names):
  if(!all(is.na(unlist(BlockNames)))) OptionsBlock <- paste(paste("\tTYPESET * UNTITLED  (CHARACTERS = ", BlockNames, ")  =  ", unlist(lapply(lapply(DataBlocks, '[[', "Ordering"), function(x) paste(paste(paste(sort(unique(x)), unlist(lapply(lapply(lapply(as.list(sort(unique(x))), '==', x), which), Zipper)), sep = ": "), collapse = ", "), sep = ""))), ";\n", sep = ""), collapse = "")
  
  # Replace cont with Squared if continuous characters present:
  if(length(grep(" cont: ", OptionsBlock)) > 0) OptionsBlock <- gsub(" cont: ", " Squared: ", OptionsBlock)
  
  # Convert continuosu character weights to one before making weights block:
  Weights[Ordering == "cont"] <- 1
  
  # Create weights block (if no block names):
  if(all(is.na(BlockNames))) WeightsBlock <- ifelse(all(Weights == 1), "", paste("\tWTSET * UNTITLED  = ", paste(paste(sort(unique(Weights)), unlist(lapply(lapply(lapply(as.list(sort(unique(Weights))), '==', Weights), which), Zipper)), sep = ": "), collapse = ", "), ";\n", sep = ""))
  
  # Create weights block (if there are block names):
  if(!all(is.na(unlist(BlockNames)))) WeightsBlock <- paste(paste("\tWTSET * UNTITLED  (CHARACTERS = ", BlockNames, ")  =  ", unlist(lapply(lapply(DataBlocks, '[[', "Weights"), function(x) paste(paste(paste(sort(unique(x)), unlist(lapply(lapply(lapply(as.list(sort(unique(x))), '==', x), which), Zipper)), sep = ": "), collapse = ", "), sep = ""))), ";\n", sep = ""), collapse = "")
  
  # Build assumptions block:
  AssumptionBlock <- paste("BEGIN ASSUMPTIONS;\n", StepMatrixBlock, OptionsBlock, WeightsBlock, "END;\n", sep = "")
  
  # Build full string with all blocks together:
  FullString <- paste("#NEXUS\n\n", HeaderBlock, TaxaBlock, DataBlock, MatrixBlock, AssumptionBlock, sep = "")

  # Write to file:
  write(FullString, filename)

}
