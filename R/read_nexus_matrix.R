#' Reads in a morphological #NEXUS data file
#'
#' @description
#'
#' Reads in a morphological data file in #NEXUS format.
#'
#' @param File A file name specified by either a variable of mode character, or a double-quoted string.
#' @param Equaliseweights Optional that overrides the weights specified in the file to make all characters truly equally weighted.
#'
#' @details
#'
#' Reads in a #NEXUS (Maddison et al. 1997) data file representing the distribution of characters (continuous, discrete, DNA etc.) in a set of taxa. Unlike \link{read.nexus.data} this function can handle polymorphisms (e.g., \code{(012)}).
#'
#' Note that the function is generally intolerant to excursions from a standard format and it is recommended your data be formatted like the \code{morphmatrix.nex} example below. However, the function also produces informative error messages if (expected) excursions are discovered.
#'
#' Previously all empty values (missing or inapplicable) were treated as NAs. But now anything coded as a "gap" now appears as an empty text string ("") in the matrix. Additionally, previously polymorphisms and uncertianties were both considered as polymorphisms with multiple states separated by an ampersand ("&"), but now polymorphisms use the ampersand ("&") and uncertainties use a slash ("/"), allowing for different treatment later and correct outputting when writing to #NEXUS format. (NB: TNT does not allow this distinction and so both polymorphisms and uncertainties will be output as polymorphisms.)
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
#' \link{build_cladistic_matrix}, \link{compactify_matrix}, \link{prune_cladistic_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @references
#'
#' Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. \emph{Systematic Biology}, \bold{46}, 590-621.
#'
#' @examples
#' 
#' # Create example matrix (NB: creates a file in the current
#' # working directory called "morphmatrix.nex"):
#' cat("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS  NTAX=5 NCHAR=5;\n\t
#' FORMAT SYMBOLS= \" 0 1 2\" MISSING=? GAP=- ;\nMATRIX\n\n
#' Taxon_1  010?0\nTaxon_2  021?0\nTaxon_3  02111\nTaxon_4  011-1
#' \nTaxon_5  001-1\n;\nEND;\n\nBEGIN ASSUMPTIONS;\n\t
#' OPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n\t
#' TYPESET * UNTITLED  = unord: 1 3-5, ord: 2;\n\t
#' WTSET * UNTITLED  = 1: 2, 2: 1 3-5;\nEND;", file = "morphmatrix.nex")
#' 
#' # Read in example matrix:
#' morph.matrix <- read_nexus_matrix("morphmatrix.nex")
#' 
#' # View example matrix in R:
#' morph.matrix
#'
#' # Remove the generated data set:
#' file.remove("morphmatrix.nex")
#' 
#' @export read_nexus_matrix
read_nexus_matrix <- function(File, Equaliseweights = FALSE) {
  
  # ADD ABILITY TO READ CHARSET LINES
  # COULD BE MULTIPLE TYPESET OR WTSET LINES, NEED TO CHECK FOR THIS!
  # ADD ABILITY TO READ AND STORE CHARACTER STATES LABELS ETC.
  # ADD ABILITY TO DEAL WITH RANGES FOR CONTINUOUS DATA
  
  # Line formatting function to be used in lapply below to deal with polymorphic characters:
  format_lines <- function(x, direction = "in") {
    
    # Split current line by character:
    CurrentString <- strsplit(x, split = "")[[1]]
    
    # Check for square brackets:
    if (any(CurrentString == "[")) stop("TNT-style square brackets, [], found in matrix. Replace with () or {} depending on whether a true polymorphism or uncertainty respectively.")
    
    # If polymorphisms (i.e., parentheses) are found:
    if (any(c(CurrentString == "{", CurrentString == "("))) {
      
      # While there are polymorphisms:
      while(any(c(CurrentString == "{", CurrentString == "("))) {
        
        # Find beginning of first polymorphism (will loop through first polymorphism until none are left):
        FirstPolymorphismBegins <- c(which(CurrentString == "{"), which(CurrentString == "("))[1]
        
        # If polymorphsim begins with a parenthesis then it must end with a parenthesis:
        if (CurrentString[FirstPolymorphismBegins] == "(") PolymorphismEndsWith <- ")"
        
        # If polymorphsim begins with a curly brace then it must end with a curly brace:
        if (CurrentString[FirstPolymorphismBegins] == "{") PolymorphismEndsWith <- "}"
        
        # Find end of first polymorphism:
        FirstPolymorphismEnds <- which(CurrentString == PolymorphismEndsWith)[1]
        
        # If a true polymorphism (two or more states observed) and reading data in use ampersand to separate states:
        if (CurrentString[FirstPolymorphismBegins] == "(" && direction == "in") PolymorphicCharacter <- paste(sort(CurrentString[(FirstPolymorphismBegins + 1):(FirstPolymorphismEnds - 1)]), collapse = "&")
        
        # If an uncertainty (two or more possible states) and reading data in use slash to separate states:
        if (CurrentString[FirstPolymorphismBegins] == "{" && direction == "in") PolymorphicCharacter <- paste(sort(CurrentString[(FirstPolymorphismBegins + 1):(FirstPolymorphismEnds - 1)]), collapse = "/")
        
        # If a true polymorphism (two or more states observed) and exporting data out use ampersand to separate states:
        if (CurrentString[FirstPolymorphismBegins] == "(" && direction == "out") PolymorphicCharacter <- paste("(", paste(sort(CurrentString[(FirstPolymorphismBegins + 1):(FirstPolymorphismEnds - 1)]), collapse = ""), ")", sep = "")
        
        # If an uncertainty (two or more possible states) and exporting data out use slash to separate states:
        if (CurrentString[FirstPolymorphismBegins] == "{" && direction == "out") PolymorphicCharacter <- paste("{", paste(sort(CurrentString[(FirstPolymorphismBegins + 1):(FirstPolymorphismEnds - 1)]), collapse = ""), "}", sep = "")
        
        # Remove now redundant characters from current string:
        CurrentString <- CurrentString[-((FirstPolymorphismBegins + 1):FirstPolymorphismEnds)]
        
        # Store polymorphic character in string:
        CurrentString[FirstPolymorphismBegins] <- PolymorphicCharacter
        
      }
      
    }
    
    # Return string with polymorphisms formatted to single characters:
    return(CurrentString)
    
  }
  
  # Subfunction to extract ordering information:
  extract_assumptions <- function(assumptionline) {
    
    # Reformat as a list:
    x <- lapply(lapply(as.list(gsub(";", "", strsplit(assumptionline, split = ", ")[[1]])), strsplit, split = ": "), unlist)
    
    # Extract just the numbers:
    Numbers <- lapply(lapply(lapply(x, '[', 2), strsplit, split = " "), unlist)
    
    # If there are hyphens (ranges) in the data:
    if (length(grep("-", Numbers)) > 0) {
      
      # Subfunction for unpacking ranges denoted by hyphens:
      unpack_ranges <- function(Numbers) {
        
        # Whilst hyphens remain in the data:
        while(length(grep("-", Numbers)) > 0) {
          
          # Find first hyphen:
          FirstHyphen <- grep("-", Numbers)[1]
          
          # Unpack and add at end whilst removing first hyphen:
          Numbers <- c(Numbers[-FirstHyphen], as.character(c(as.numeric(strsplit(Numbers[FirstHyphen], split = "-")[[1]])[1]:as.numeric(strsplit(Numbers[FirstHyphen], split = "-")[[1]])[2])))
          
        }
        
        # Convert to actual numbers and sort back into increasing order:
        Numbers <- sort(as.numeric(Numbers))
        
        # Return unpacked numbers:
        return(Numbers)
        
      }
      
      # Apply unpacking function to data:
      Numbers <- lapply(Numbers, unpack_ranges)
      
    # If there are no hyphens:
    } else {
      
      # Simply convert to numbers and ensure they are sorted in increasing order:
      Numbers <- lapply(lapply(Numbers, as.numeric), sort)
      
    }
    
    # Add character type as names to numbers list:
    names(Numbers) <- unlist(lapply(x, '[', 1))
    
    # Build empty vector to store output:
    output <- vector(mode = "character", length = max(unlist(Numbers)))
    
    # Store character type in output:
    for(i in 1:length(Numbers)) output[Numbers[[i]]] <- names(Numbers)[i]
    
    # Check all characters have some kind of ordering described:
    if (length(which(output == "")) > 0) stop(paste("The following characters have no specified ordering: ", paste(which(output == ""), collapse = ", "), ". Check ASSUMPTIONS block and edit appropriately.", sep = ""))
    
    # Return output:
    return(output)
    
  }
  
  # Subfunction to find range of states for each character:
  find_ranges <- function(x) {
    
    # Convert each column of matrix to a list of numeric values:
    x <- lapply(lapply(lapply(lapply(lapply(apply(x, 2, as.list), unlist), strsplit, split = "&|/"), unlist), as.numeric), sort)
    
    # Check for any empty columns (all NAs) and insert a dummy 0 if found:
    if (any(lapply(x, length) == 0)) x[which(lapply(x, length) == 0)] <- 0
    
    # Convert to a min-max matrix:
    x <- matrix(unlist(lapply(x, range)), ncol = 2, byrow = TRUE, dimnames = list(c(), c("Min", "Max")))
    
    # Return(x):
    return(x)
    
  }
  
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

  # Little piece of code to help deal with foreign characters:
  Sys.setlocale('LC_ALL','C')

  # Read in NEXUS file as raw text:
  X <- readLines(File, warn = FALSE)

  # Check that this is a #NEXUS file:
  if (length(grep("#NEXUS", X, ignore.case = TRUE)) == 0) stop("This is not a #NEXUS file.")

  # Replace any lower case matrix with upper case:
  if (length(grep("MATRIX", X, ignore.case = FALSE)) == 0) X <- gsub("matrix", "MATRIX", X, ignore.case = FALSE)

  # Check that the #NEXUS file includes a matrix:
  if (length(grep("MATRIX", X, ignore.case = TRUE)) == 0) stop("There is no matrix present.")

  # Are there blocks of the NEXUS file that can be removed? (i.e., trees, macclade, mesquite etc.:
  if (length(grep("begin trees|begin mesquite|begin mrbayes|begin macclade|begin treebase|begin notes|begin paup", X, ignore.case = TRUE)) > 0) {

    # Get position of all Begin tags:
    all.begins <- grep("begin ", X, ignore.case = TRUE)

    # Get position of all End tags:
    all.ends <- grep("end;|endblock;|end ;|endblock ;", X, ignore.case = TRUE)

    # Check Begin and End are in balance:
    if (length(all.begins) != length(all.ends)) stop("Begin and End tags are not equal in number.")

    # Get rows for blocks to delete:
    block.rows <- match(grep("begin trees|begin mesquite|begin mrbayes|begin macclade|begin treebase|begin notes|begin paup", X, ignore.case = TRUE), all.begins)

    # Find start lines:
    block.begins <- all.begins[block.rows]

    # Find end lines:
    block.ends <- all.ends[block.rows]

    # Incrementally remove superflouous blocks:
    for(i in length(block.begins):1) X <- X[c(c(1:(block.begins[i] - 1)), c((block.ends[i] + 1):length(X)))]

  }

  # Case if there are superflous taxon, character, or state labels:
  if (length(grep("taxlabels|charstatelabels|charlabels|statelabels", X, ignore.case = TRUE)) > 0) {

    # Find label beginning rows:
    label.ends <- label.begins <- grep("taxlabels|charstatelabels|charlabels|statelabels", X, ignore.case = TRUE)

    # For each label find the end tag that corresponds to it:
    for(i in 1:length(label.begins)) label.ends[i] <- min(grep(";", X, ignore.case = TRUE)[grep(";", X, ignore.case = TRUE) > label.begins[i]])

    # Incrementally remove superflouous labels:
    for(i in length(label.begins):1) X <- X[c(c(1:(label.begins[i] - 1)), c((label.ends[i] + 1):length(X)))]

  }

  # If there is text inside single quotes:
  if (length(grep("'", X)) > 0) {

    # Find items to replace (text not in single quotes):
    replacement.items <- unique(unlist(strsplit(gsub(paste("'[A-Z|\t| |0-9|a-z|\\?|\\-|\\(|)|\\.]*'", sep=""), "''", X[grep("'", X)]), "''")))

    # Make sure nothing "" is not in this list:
    replacement.items <- setdiff(unique(unlist(strsplit(replacement.items, " "))), "")

    # Make block that isolates lines with single quotes:
    lines.to.edit <- X[grep("'", X)]

    # Remove text outside single quotes:
    for(i in replacement.items[order(nchar(replacement.items), decreasing=TRUE)]) lines.to.edit <- gsub(i, "", lines.to.edit, fixed=T)

    # Remove double spaces:
    while(length(grep("  ", lines.to.edit)) > 0) lines.to.edit <- gsub("  ", " ", lines.to.edit, fixed=T)

    # Now isolate names within single quotes:
    names.in.single.quotes <- unique(gsub("'", "", strsplit(gdata::trim(paste(lines.to.edit, collapse="")), "' '")[[1]]))

    # Make sure nothing "" is not in this list:
    names.in.single.quotes <- which(!names.in.single.quotes == "")

    # Replace names in single quotes with underscored versions without single quotes or parentheses:
    for(i in names.in.single.quotes) X <- gsub(i, gsub("\\(|)", "", gsub(" ", "_", i)), X, fixed=T)

    # Remove all single quotes from text file:
    X <- gsub("'", "", X, fixed = TRUE)

  }

  # Get rid of spaces around equals signs to make later regular expresssions work:
  while(length(grep(" = |= | =", X))) X <- gsub(" = |= | =", "=", X)

  # Remove weird characters (may need to add to the list or not if Sys.setlocale fixes the issues):
  X <- gsub("\x94|\x93|\xd5|\xd4|\xd3|\xd2|'", "", X)

  # Replace tabs with spaces and trim leading and trailing spaces from each line:
  X <- apply(matrix(gsub("\t", " ", X)), 2, gdata::trim)

  # Delete any empty lines (if present):
  if (length(which(X == "")) > 0) X <- X[-which(X == "")]
  
  # If header text is present:
  if (length(grep("\\[", X)) > 0) {

    # Work out beginning and ending lines for text:
    textlines <- apply(cbind(setdiff(grep("\\[", X), grep("\\[[0-9:A-Z:a-z]{1}\\]", X)), setdiff(grep("\\]", X), grep("\\[[0-9:A-Z:a-z]{1}\\]", X))), 1, paste, collapse = ":")

    # Convert beginning and endings to numerics for finding in vector:
    lines.to.delete <- textlines <- eval(parse(text = paste("c(", paste(textlines, collapse = ","), ")", sep = "")))
    
    # Find only lines beginning with a square bracket to avoid issue with polymorphic use of square brackets (can be caught later):
    LinesToKeep <- which(unlist(lapply(strsplit(X[textlines], split = ""), '[', 1)) == "[")
    
    # Reduce textlines appropriately:
    textlines <- textlines[LinesToKeep]
    
    # Reduce lines to delete appropriately:
    lines.to.delete <- lines.to.delete[LinesToKeep]
  
    # Grab text and store:
    textlines <- paste(gsub("\\[|\\]", "", X[textlines]), collapse = "\n")
  
    # Delete text lines to avoid regular expression errors later:
    X <- X[-lines.to.delete]

  # If header text is absent:
  } else {
    
    # Store empty text string:
    textlines <- ""

  }

  # Little test for valid NEXUS formatting for number of taxa:
  if (length(grep("ntax", X, ignore.case = TRUE)) == 0) stop("Number of taxa not defined.")

  # Little test for valid NEXUS formatting for number of characters:
  if (length(grep("nchar", X, ignore.case = TRUE)) == 0) stop("Number of characters not defined.")
  
  # Get number of taxa:
  NTax <- as.numeric(strsplit(strsplit(gsub(";|\t", "", X[grep("NTAX=|Ntax=|ntax=", X)]), "NTAX=|Ntax=|ntax=")[[1]][2], " ")[[1]][1])
  
  # Create empty list to store matrix block(s):
  MatrixBlockList <- list()

  # If there are character blocks (rather than a data block):
  if (length(grep("begin characters", X, ignore.case = TRUE)) > 0) {
    
    # Find line where character block(s) begin:
    CharacterBlockBegins <- grep("begin characters", X, ignore.case = TRUE)
    
    # Find line where matching character block(s) end:
    CharacterBlockEnds <- grep("end;", X, ignore.case = TRUE)[unlist(lapply(lapply(lapply(as.list(CharacterBlockBegins), '<', grep("end;", X, ignore.case = TRUE)), which), min))]
    
    # For each character block:
    for(i in 1:length(CharacterBlockBegins)) {
      
      # isolate ith character block:
      CurrentCharacterBlock <- X[CharacterBlockBegins[i]:CharacterBlockEnds[i]]
      
      # Locate dataype line:
      datatypeLine <- grep("datatype=", CurrentCharacterBlock, ignore.case = TRUE)
      
      # Check there is a dataype:
      if (length(datatypeLine) == 0) stop("Character datatype not specified.")
      
      # Check there are not multiple datatypes for a single character block:
      if (length(datatypeLine) > 1) stop("Multiple character datatypes specified.")
      
      # Get ith datatype:
      Currentdatatype <- strsplit(strsplit(CurrentCharacterBlock[datatypeLine], "datatype=|DATATYPE=")[[1]][2], " ")[[c(1,1)]]
      
      # Add ith character block to list:
      MatrixBlockList[[i]] <- CurrentCharacterBlock
      
      # Add datatype to names for ith matrix block:
      names(MatrixBlockList)[i] <- Currentdatatype

    }
    
  }
  
  # If there is a data block:
  if (length(grep("begin data", X, ignore.case = TRUE)) > 0) {
    
    # Get beginning line of data block:
    BeginDataLine <- grep("begin data", X, ignore.case = TRUE)
    
    # Get end line of data block:
    EndDataLine <- grep("end;", X, ignore.case = TRUE)[unlist(lapply(lapply(lapply(as.list(BeginDataLine), '<', grep("end;", X, ignore.case = TRUE)), which), min))]
    
    # Add data block to list:
    MatrixBlockList[[1]] <- X[BeginDataLine:EndDataLine]
    
    # Locate datatype line:
    datatypeLine <- grep("datatype=", MatrixBlockList[[1]], ignore.case = TRUE)
    
    # If there is a stated datatype:
    if (length(datatypeLine) > 0) {
      
      # Store stated datatype as list name:
      names(MatrixBlockList)[1] <- strsplit(strsplit(MatrixBlockList[[1]][datatypeLine], "datatype=|DATATYPE=")[[1]][2], " ")[[c(1,1)]]
      
    # If there is no stated datatype:
    } else {
      
      # Set datatype as "STANDARD":
      names(MatrixBlockList)[1] <- "STANDARD"
      
    }

  }
  
  # If there are any interleaved block(s):
  if (length(unlist(lapply(MatrixBlockList, grep, pattern = "INTERLEAVE"))) > 0) {
    
    # Get interleaved blocks:
    InterleavedBlocks <- which(unlist(lapply(MatrixBlockList, grep, pattern = "INTERLEAVE") > 0))
    
    # For each interleaved block:
    for(i in InterleavedBlocks) {
      
      # Isolate current matrix block:
      CurrentBlock <- MatrixBlockList[[i]]
      
      # Remove interleave text:
      CurrentBlock <- gsub("INTERLEAVE ", "", CurrentBlock)
      
      # Find line matrix begins on:
      MatrixBegins <- which(toupper(CurrentBlock) == "MATRIX") + 1
      
      # Find line matrix ends:
      MatrixEnds <- rev(which(CurrentBlock == ";")) - 1
      
      # Get number of blocks in interleaved matrix:
      NumberOfBlocks <- length(MatrixBegins:MatrixEnds) / NTax
      
      # Check is divisible by number of taxa:
      if ((NumberOfBlocks %% 1) != 0) stop("Interleaved matrix has greater or fewer lines than a multiple of number of taxa.")
      
      # Empty lst to store interleaved blocks:
      InterleavedBlockList <- list()
      
      # For each interleaved block:
      for(j in 1:NumberOfBlocks) {
        
        # Get end of block:
        BlockEnds <- j * NTax
        
        # Get beginning of block:
        BlockBegins <- BlockEnds - NTax + 1
        
        # Isolate block:
        JBlock <- CurrentBlock[MatrixBegins:MatrixEnds][BlockBegins:BlockEnds]
        
        # Remove any double spaces:
        while(length(grep("  ", JBlock)) > 0) JBlock <- gsub("  ", " ", JBlock)
        
        # Check for rogue spaces and stop and warn if found:
        if (any(unlist(lapply(strsplit(JBlock, split = " "), length)) > 2)) stop("Rogue spaces found in interleaved block. Check taxon names and polymorphisms and remove.")
        
        # Store jth block in interleaved blocks:
        InterleavedBlockList[[j]] <- JBlock
        
      }
      
      # Reformat as matrices:
      InterleavedBlockList <- lapply(lapply(lapply(InterleavedBlockList, strsplit, split = " "), unlist), matrix, ncol = 2, byrow = TRUE)
      
      # Combine as single large matrix:
      InterleavedBlockList <- matrix(unlist(lapply(InterleavedBlockList, '[', , 2, drop = FALSE)), nrow = NTax, byrow = FALSE, dimnames = list(InterleavedBlockList[[1]][, 1], c()))
      
      # Collapse matrix to single lines:
      InterleavedBlockList <- apply(cbind(cbind(rownames(InterleavedBlockList), rep(" ", NTax)), InterleavedBlockList), 1, paste, collapse = "")
      
      # Overwrite matrix block with now de-interleaved version:
      MatrixBlockList[[i]] <- c(MatrixBlockList[[i]][1:(MatrixBegins - 1)], InterleavedBlockList, MatrixBlockList[[i]][(MatrixEnds + 1):length(MatrixBlockList[[i]])])
      
    }
    
  }
  
  # If any matrix blocks are of type "MIXED":
  if (any(lapply(lapply(lapply(strsplit(names(MatrixBlockList), split = ""), '[', 1:5), paste, collapse = ""), "toupper") == "MIXED")) {
    
    # Get any mixed datatype block(s):
    MixedBlocks <- which(lapply(lapply(lapply(strsplit(names(MatrixBlockList), split = ""), '[', 1:5), paste, collapse = ""), "toupper") == "MIXED")
    
    # For each mixed block:
    for(i in MixedBlocks) {
      
      # Isolate current block:
      CurrentBlock <- MatrixBlockList[[i]]
      
      # Reduce to MATRIX block:
      CurrentBlock <- CurrentBlock[(which(toupper(CurrentBlock) == "MATRIX") + 1):(rev(which(CurrentBlock == ";")) - 1)]
      
      # Remove any double spaces:
      while(length(grep("  ", CurrentBlock))) CurrentBlock <- gsub("  ", " ", CurrentBlock)
      
      # Check for any rogue spaces:
      if (any(unlist(lapply(strsplit(CurrentBlock, split = " "), length)) > 2)) stop("Rogue spaces found in matrix block. Check taxon names and polymorphisms.")
      
      # Isolate taxon names:
      TaxonNames <- matrix(unlist(strsplit(CurrentBlock, split = " ")), ncol = 2, nrow = NTax, byrow = TRUE)[, 1]
      
      # Reformat data block as list:
      CurrentBlock <- as.list(matrix(unlist(strsplit(CurrentBlock, split = " ")), ncol = 2, nrow = NTax, byrow = TRUE)[, 2])
      
      # Convert data block into individual character vectors:
      CurrentBlock <- lapply(CurrentBlock, format_lines, direction = "out")
      
      # Add taxon names to list:
      names(CurrentBlock) <- TaxonNames
      
      # Isolate components of mixed block:
      MixedComponents <- strsplit(names(MatrixBlockList)[i], "\\(|\\)")[[1]][2]
      
      # Isolate further (remove spaces so components can be isolated):
      MixedComponents <- strsplit(gsub(" ", "", MixedComponents), split = ",")[[1]]
      
      # Format as list:
      MixedComponents <- lapply(strsplit(MixedComponents, split = ":"), toupper)
      
      # Get total number of characters in block:
      Totalcharacters <- max(as.numeric(unlist(lapply(lapply(MixedComponents, '[', 2), strsplit, split = "-")))) - min(as.numeric(unlist(lapply(lapply(MixedComponents, '[', 2), strsplit, split = "-")))) + 1
      
      # Check all rows have the same (correct) number of characters:
      if (any(!unlist(lapply(CurrentBlock, length)) == Totalcharacters)) stop("Some lines of matrix have too many or too few characters.")
      
      # For each component:
      for(j in 1:length(MixedComponents)) {
        
        # Isolate current partition as matrix:
        CurrentMatrix <- matrix(unlist(lapply(CurrentBlock, '[', as.numeric(strsplit(MixedComponents[[j]][2], split = "-")[[1]])[1]:as.numeric(strsplit(MixedComponents[[j]][2], split = "-")[[1]])[2])), nrow = NTax, byrow = TRUE, dimnames = list(names(CurrentBlock), c()))
        
        # Add component of matrix to end of list:
        MatrixBlockList[[(length(MatrixBlockList) + 1)]] <- c(toupper(MatrixBlockList[[i]][1:which(toupper(MatrixBlockList[[i]]) == "MATRIX")]), apply(cbind(cbind(rownames(CurrentMatrix), rep(" ", NTax)), CurrentMatrix), 1, paste, collapse = ""), toupper(MatrixBlockList[[i]][rev(which(MatrixBlockList[[i]] == ";"))[1]:length(MatrixBlockList[[i]])]))
        
        # Remove names:
        names(MatrixBlockList[[length(MatrixBlockList)]]) <- NULL
        
        # Add datatype name as name to new list item:
        names(MatrixBlockList)[length(MatrixBlockList)] <- MixedComponents[[j]][1]
        
        # Get number of characters in partition:
        charactersInPartition <- max(as.numeric(strsplit(MixedComponents[[j]][2], split = "-")[[1]])) - min(as.numeric(strsplit(MixedComponents[[j]][2], split = "-")[[1]])) + 1
        
        # Update number of characters for parition:
        MatrixBlockList[[length(MatrixBlockList)]] <- gsub(pattern = paste("NCHAR=", Totalcharacters, sep = ""), replacement = paste("NCHAR=", charactersInPartition, sep = ""), MatrixBlockList[[length(MatrixBlockList)]], fixed = TRUE)
        
        # Update datatype for partition:
        MatrixBlockList[[length(MatrixBlockList)]] <- gsub(pattern = paste("DATATYPE=", toupper(names(MatrixBlockList)[i]), sep = ""), replacement = paste("DATATYPE=", MixedComponents[[j]][1], sep = ""), MatrixBlockList[[length(MatrixBlockList)]], fixed = TRUE)
        
      }
      
    }
    
    # Remove now redundant mixed blocks-:
    MatrixBlockList <- MatrixBlockList[-MixedBlocks]

  }
  
  # Will need to know datatypes before proceeding so check for nonstandard ones now:
  Nonstandarddatatypes <- setdiff(toupper(names(MatrixBlockList)), c("CONTINUOUS", "DNA", "NUCLEOTIDE", "PROTEIN", "RESTRICTION", "RNA", "STANDARD"))
  
  # Stop if non-standard datatype(s) found:
  if (length(Nonstandarddatatypes) > 0) stop("Non-standard datatype found (i.e., neither CONTINUOUS, DNA, NUCLEOTIDE, PROTEIN, RESTRICTION, RNA, or STANDARD).")
  
  # Create empty block names vector:
  block_names <- rep(NA, length(MatrixBlockList))
  
  # For each matrix block:
  for(i in 1:length(MatrixBlockList)) {
    
    # Isolate ith block:
    CurrentBlock <- MatrixBlockList[[i]]
    
    # Get title line (if found) - should work even if there is a taxon whose name begins "TITLE":
    TitleLine <- intersect(which(unlist(lapply(lapply(strsplit(CurrentBlock, split = ""), '[', 1:5), paste, collapse = "")) == "TITLE"), 1:which(unlist(lapply(lapply(strsplit(CurrentBlock, split = ""), '[', 1:5), paste, collapse = "")) == "MATRI")[1])
    
    # Store block name:
    if (length(TitleLine) > 0) block_names[i] <- gsub(";", "", rev(strsplit(CurrentBlock[TitleLine], split = " ")[[1]])[1])
    
  }

  # Get number of characters in each block:
  NChars <- lapply(lapply(lapply(lapply(lapply(lapply(lapply(lapply(lapply(MatrixBlockList, function(x) x[grep("NCHAR=|Nchar=|nchar=", x)]), gsub, pattern = ";|\t", replacement = ""), strsplit, split = "NCHAR=|Nchar=|nchar="), '[[', 1), '[', 2), strsplit, split = " "), '[[', 1), '[', 1), as.numeric)
  
  # Check matri(ces) have dimension:
  if (!(all(NChars > 0) && NTax > 0)) stop("One or more matrix blocks have no dimensions (i.e., zero characters or taxa).")
  
  # Function to get symbols from matrix block:
  get_symbols <- function(X) {
    
    # If symbols are specified in the file:
    if (length(grep("symbols", X, ignore.case = TRUE)) > 0) {
      
      # Get initial set of symbols:
      symbols <- strsplit(strsplit(strsplit(X[grep("symbols", X, ignore.case = TRUE)], "SYMBOLS=|symbols=|symbols=")[[1]][2], "\"")[[1]][2], " |")[[1]]
      
      # Collapse to just the symbols themselves:
      symbols <- symbols[nchar(symbols) == 1]
      
      # Special case of a tilda indicating a range:
      if (length(symbols) == 3 && symbols[2] == "~") symbols <- as.character(as.numeric(symbols[1]):as.numeric(symbols[3]))
      
    # If symbols are not specified in the file:
    } else {
      
      # Find datatype row (if specified):
      datatypeRow <- grep("datatype", X, ignore.case = TRUE)
      
      # If a datatype was specified:
      if (length(datatypeRow) > 0) {
        
        # Store datatype:
        datatype <- strsplit(strsplit(toupper(X[datatypeRow]), split = "DATATYPE=")[[1]][2], split = " ")[[1]][1]
        
      # If datatype is not specified:
      } else {
        
        # Set datatype as STANDARD:
        datatype <- "STANDARD"
        
      }
      
      # If CONTINUOUS use default symbols::
      if (datatype == "CONTINUOUS") symbols <- NULL
      
      # If DNA use default symbols::
      if (datatype == "DNA") symbols <- c("A", "C", "G", "T", "R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N")
      
      # If NUCLEOTIDE use default symbols::
      if (datatype == "NUCLEOTIDE") symbols <- c("A", "C", "G", "T", "R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N")
      
      # If PROTEIN use default symbols::
      if (datatype == "PROTEIN") symbols <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
      
      # If RESTRICTION use default symbols::
      if (datatype == "RESTRICTION") symbols <- c("0", "1")
      
      # If STANDARD use default symbols::
      if (datatype == "RNA") symbols <- c("A", "C", "G", "U", "R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N")
      
      # If STANDARD use default symbols::
      if (datatype == "STANDARD") symbols <- c(c(0:9), LETTERS[1:22])
      
    }
    
    # Return symbols:
    return(symbols)
    
  }
  
  # Get symbols from each block:
  symbols <- lapply(MatrixBlockList, get_symbols)
  
  # Get missing character function:
  get_missing <- function(X) {
    
    # If the missing character is specified:
    if (length(grep("missing", X, ignore.case = TRUE)) > 0) {
      
      # Get missing character:
      missing <- strsplit(strsplit(X[grep("missing", X, ignore.case = TRUE)][1], "MISSING=|missing=|missing=")[[1]][2], "")[[1]][1]
      
    # If the missing character is not specified:
    } else {
      
      # Set as default symbol of the question mark:
      missing <- "?"
      
    }
    
    # Return gap character:
    return(missing)

  }
  
  # Get missing symbol(s):
  missing <- lapply(MatrixBlockList, get_missing)
  
  # Get gap character function:
  get_gap <- function(X) {
    
    # If the gap character is specified:
    if (length(grep("gap", X[1:grep("MATRIX", X)[1]], ignore.case = TRUE)) > 0) {
      
      # Get gap character:
      gap <- strsplit(strsplit(X[grep("gap", X, ignore.case = TRUE)][1], "GAP=|gap=|gap=")[[1]][2], "")[[1]][1]
      
    # If the gap character is not specified:
    } else {
      
      # Set as default symbol of the dash:
      gap <- "-"
      
    }
    
    # Return gap character:
    return(gap)
    
  }
  
  # Get gap symbol(s):
  gap <- lapply(MatrixBlockList, get_gap)
  
  # Find first line of each matrix:
  MatrixStartLines <- lapply(lapply(lapply(MatrixBlockList, '==', "MATRIX"), which), '+', 1)
  
  # Find alst line of each matrix:
  MatrixEndLines <- lapply(lapply(lapply(MatrixBlockList, '==', ";"), which), '-', 1)
  
  # Cut down matrices to just matrix block:
  for(i in 1:length(MatrixBlockList)) MatrixBlockList[[i]] <- MatrixBlockList[[i]][MatrixStartLines[[i]]:MatrixEndLines[[i]]]
  
  # Stop if incorrect number of taxa found:
  if (any(lapply(MatrixBlockList, length) != NTax)) stop("Some matrix block(s) have too many or too few taxa. Check for linebreaks or incorrect NTAX value.")
  
  # If any blocks are of type "CONTINUOUS":
  if (any(toupper(names(MatrixBlockList)) == "CONTINUOUS")) {
    
    # Get block(s) of type "CONTINUOUS":
    ContinuousBlocks <- which(toupper(names(MatrixBlockList)) == "CONTINUOUS")
    
    # For each block of type "CONTINUOUS":
    for(i in ContinuousBlocks) {
      
      # Isolate ith continuous matrix block:
      CurrentBlock <- MatrixBlockList[[i]]
      
      # While there are double spaces (that will complicate strsplit later) replace these with single spaces:
      while(length(grep("  ", CurrentBlock)) > 0) CurrentBlock <- gsub("  ", " ", CurrentBlock)
      
      # Check all rows are equal in length and stop and warn if not:
      if (length(unique(unlist(lapply(strsplit(CurrentBlock, split = " "), length)))) > 1) stop("Continuous block has missing codings for some taxa (check for missing spaces or spaces in taxon names).")
      
      # Format data as matrix with rownames as taxa:
      CurrentBlock <- matrix(unlist(strsplit(CurrentBlock, split = " ")), nrow = NTax, byrow = TRUE, dimnames = list(unlist(lapply(strsplit(CurrentBlock, split = " "), '[', 1)), c()))
      
      # Remove names (first column) from matrix:
      CurrentBlock <- CurrentBlock[, -1, drop = FALSE]
      
      # Replace missing values with NA (if there are any):
      if (length(grep(missing[[ContinuousBlocks]], CurrentBlock)) > 0) CurrentBlock <- gsub(missing[[ContinuousBlocks]], NA, CurrentBlock, fixed = TRUE)
      
      # Store newly formatted matrix back into matrix block list:
      MatrixBlockList[[i]] <- CurrentBlock
      
    }
    
  }
  
  # If there are non-continuous blocks:
  if (length(setdiff(toupper(names(MatrixBlockList)), "CONTINUOUS")) > 0) {
    
    # Get non-continuous block number(s):
    NonContinuousBlocks <- which(!toupper(names(MatrixBlockList)) == "CONTINUOUS")
    
    # For each non-continuous block:
    for(i in NonContinuousBlocks) {
      
      # Isolate ith non-continuous matrix block:
      CurrentBlock <- MatrixBlockList[[i]]
      
      # While there are double spaces (that will complicate strsplit later) replace these with single spaces:
      while(length(grep("  ", CurrentBlock)) > 0) CurrentBlock <- gsub("  ", " ", CurrentBlock)
      
      # Check for rogue spaces (should just be one between taxon name and actual data (i.e., none in taxon names or polymorphisms):
      if (any(unlist(lapply(strsplit(CurrentBlock, split = " "), length)) > 2)) stop("Rogue spaces found (check taxon names and polymorphisms and remove before trying again).")
      
      # Reformat as one-column matrix with taxa as row names:
      CurrentBlock <- matrix(unlist(strsplit(CurrentBlock, split = " ")), ncol = 2, byrow = TRUE, dimnames = list(lapply(strsplit(CurrentBlock, split = " "), '[', 1), c()))[, 2, drop = FALSE]
      
      # Overwrite curentmatrix block with newly formatted matrix:
      MatrixBlockList[[i]] <- CurrentBlock
      
      # Format current block as a list (loses taxon names for now hence storage line above):
      CurrentBlock <- as.list(CurrentBlock)
      
      # Apply line formatting function to all rows:
      CurrentBlock <- lapply(CurrentBlock, format_lines, direction = "in")
      
      # If any rows have too many or too few characters:
      if (any(lapply(CurrentBlock, length) != NChars[[i]])) {
        
        # Isolate problem rows:
        ProblemRows <- which(unlist(lapply(CurrentBlock, length)) != NChars[[i]])
        
        # Stop and warn user:
        stop(paste("The following tax(a) have too many or too few characters: ", paste(rownames(MatrixBlockList[[i]])[ProblemRows], collapse = ", "), ". Check rows or NCHAR.", sep = ""))
        
      }
      
      # Store formatted matrix in matrix block list:
      MatrixBlockList[[i]] <- matrix(unlist(CurrentBlock), nrow = NTax, ncol = NChars[[i]], byrow = TRUE, dimnames = list(rownames(MatrixBlockList[[i]]), c()))
      
      # Find any unexpected (rogue) characters in matrix:
      Roguecharacters <- setdiff(unique(unlist(strsplit(unique(as.vector(MatrixBlockList[[i]])), split = "&|/"))), c(gap[[i]], missing[[i]], symbols[[i]]))
      
      # If any rogue characters were found stop and warn user:
      if (length(Roguecharacters) > 0) stop(paste("The following rogue characters were found: ", paste(Roguecharacters, collapse = ", "), ". Consider adding these to SYMBOLS or check if they were intended to be replaced with polymorphisms in the source matrix.", sep = ""))
      
      # Replace missing values with NA (if there are any):
      if (length(grep(missing[[i]], MatrixBlockList[[i]])) > 0) MatrixBlockList[[i]] <- gsub(missing[[i]], NA, MatrixBlockList[[i]], fixed = TRUE)

    }
    
  }
  
  # Get rownames from each matrix block (to check for any issues with these now they are isolated):
  RowNames <- lapply(MatrixBlockList, rownames)
  
  # Check for cuplicate taxon names and stop and warn if found:
  if (any(unlist(lapply(RowNames, duplicated)))) stop(paste("The following taxon name(s) are duplicated:", paste(unlist(RowNames)[which(unlist(lapply(RowNames, duplicated)))], collapse = ", "), ". Rename so that all taxon names are unique.", sep = ""))
  
  # Check names match across matrix blocks and stop and warn if not:
  if (length(unlist(lapply(RowNames, setdiff, unique(unlist(RowNames))))) > 0) stop(paste("Taxa do no match across matrx blocks. Check the following names:", paste(unlist(lapply(RowNames, setdiff, unique(unlist(RowNames)))), collapse = ", "), ".", sep = ""))
  
  # Store taxon names:
  TaxonNames <- unique(unlist(RowNames))
  
  # Find and special (non-underscore or alphanumeric) characters in taxon names:
  SpecialcharactersInTaxonNames <- unique(strsplit(paste(gsub("[:A-Z:a-z:0-9:_:]", "", TaxonNames), collapse = ""), split = "")[[1]])
  
  # If any special characters are found stop and warn user:
  if (length(SpecialcharactersInTaxonNames) > 0) stop(paste("The following special characters were found in the taxon names: ", paste(SpecialcharactersInTaxonNames, collapse = ", "), ". Remove or replace and try again.", sep = ""))
  
  # Get rid of spaces around dashes to make later regular expresssions work:
  while(length(grep(" - |- | -", X))) X <- gsub(" - |- | -", "-", X)
  
  # Create null variable for step matrices:
  step_matrices <- NULL
  
  # Special case if there are user-defined characters, i.e. step matrices:
  if (length(grep("USERTYPE", X, ignore.case = TRUE)) > 0) {
    
    # Get rows corresponding to start of stepmatrices:
    StepMatrixRows <- grep("USERTYPE", X, ignore.case = TRUE)
    
    # Create empty list to store step matrices:
    step_matrices <- list()
    
    # Little check in case of abnormally large number of step matrices:
    if (length(StepMatrixRows) > length(LETTERS)) stop("Not enough letters to store the number of different step matrices.")
    
    # For each step matrix:
    for(i in StepMatrixRows) {
      
      # Grab text block corresponding to step matrix:
      StepMatrixBlock <- X[i:(i + as.numeric(strsplit(X[i], "\\(STEPMATRIX\\)=")[[1]][2]) + 1)]
      
      # Remove [] labels for rows:
      StepMatrixBlock <- gsub("\\[[:A-Z:]\\] |\\[[:0-9:]\\] ", "", StepMatrixBlock)
      
      # Get the step matrix as a matrix (might still need to think about how to define e.g. infinity):
      StepMatrix <- gsub("\\.", "0", matrix(unlist(strsplit(StepMatrixBlock[3:length(StepMatrixBlock)], " ")), ncol = as.numeric(strsplit(X[i], "\\(STEPMATRIX\\)=")[[1]][2]), byrow = TRUE))
      
      # Add row and column names (assumes they start at zero and climb from there - i.e., should match symbols in order from 0 to N - 1):
      rownames(StepMatrix) <- colnames(StepMatrix) <- as.character(0:(ncol(StepMatrix) - 1))
      
      # Add step matrix to list:
      step_matrices[[length(step_matrices) + 1]] <- StepMatrix
      
      # Find step matrix name:
      StepMatrixName <- strsplit(strsplit(X[i], "USERTYPE ")[[1]][2], " ")[[1]][1]
      
      # Replace step matrix name with step_N:
      X <- gsub(StepMatrixName, paste("step_", LETTERS[length(step_matrices)], sep = ""), X)
      
      # Use step matrix name (step_N) in list for later calling:
      names(step_matrices)[length(step_matrices)] <- paste("step_", LETTERS[length(step_matrices)], sep="")
      
    }
    
  }
  
  # Set default weights as 1:
  weights <- lapply(lapply(MatrixBlockList, ncol), rep, x = 1)
  
  # Set default ordering as unordered:
  ordering <- lapply(lapply(MatrixBlockList, ncol), rep, x = "unord")
  
  # For each matrix block:
  for(i in 1:length(MatrixBlockList)) {
    
    # For each symbol replace with a number from 0 to N - 1 symbols (unless continuous which would be NULL for symbols):
    if (!is.null(symbols[[i]][1])) for(j in 1:length(symbols[[i]])) MatrixBlockList[[i]] <- gsub(symbols[[i]][j], as.character(j - 1), MatrixBlockList[[i]], fixed = TRUE)
    
    # Convert gap character(s) into empty text strings (""):
    MatrixBlockList[[i]] <- gsub(pattern = gap[[i]], replacement = "", MatrixBlockList[[i]], fixed = TRUE)
    
  }
  
  # Get minimum and maximum values for each character in each matrix:
  MinMaxMatrixList <- lapply(MatrixBlockList, find_ranges)
  
  # Now min-max is known need to check for continuous characters to amek sure default weights are all effectively 1:
  if (any(names(MatrixBlockList) == "CONTINUOUS")) {
    
    # Get numbers of continuous blocks:
    ContinuousBlocks <- which(names(MatrixBlockList) == "CONTINUOUS")
    
    # For each continuous blocks set weights as reciprocal of difference between min and max (i.e., effectively setting all weights as one):
    for(i in ContinuousBlocks) weights[[i]] <- as.numeric(gsub(Inf, 1, 1 / (MinMaxMatrixList[[i]][, "Max"] - MinMaxMatrixList[[i]][, "Min"])))
    
  }

  # Set default ordering as unord:
  Defaultordering <- "unord"
  
  # If deafult ordering is specified, store it:
  if (length(grep("DEFTYPE", toupper(X))) > 0) Defaultordering <- strsplit(strsplit(X[grep("deftype", X, ignore.case = TRUE)], "DEFTYPE=|Deftype=|deftype=")[[1]][2], " ")[[1]][1]
  
  # If default ordering is ordered then update ordering as "ord":
  if (Defaultordering == "ord") ordering <- lapply(lapply(MatrixBlockList, ncol), rep, x = "ord")
  
  # If one or more blocks are of type CONTINUOUS:
  if (any(toupper(names(MatrixBlockList)) == "CONTINUOUS")) {
    
    # Get continuous blocks:
    ContinuousBlocks <- which(toupper(names(MatrixBlockList)) == "CONTINUOUS")
    
    # Update continuous ordering to "cont" to denote characters are tp be treated as continuous:
    for(i in ContinuousBlocks) ordering[[i]] <- gsub("unord", "cont", ordering[[i]])
    
  }
  
  # Check assumptions block is supplied (to catch issue of character information appearing in a different block type (e.g., MRBAYES):
  if (length(grep("begin assumptions", X, ignore.case = TRUE)) == 0) stop("No assumptions block specified. NB: Claddis can not read ordering information stored in another type of block (e.g., MRBAYES).")
  
  # Find any TYPESET lines:
  TypesetLines <- X[lapply(lapply(strsplit(X, split = ""), '[', 1:7), paste, collapse = "") == "TYPESET"]
  
  # If there are typeset line(s):
  if (length(TypesetLines) > 0) {
    
    # Remove any double spaces found (fot easier strsplit later):
    while(length(grep("  ", TypesetLines))) TypesetLines <- gsub("  ", " ", TypesetLines)
    
    # If there are any block names (i.e., names that may be used to differentiate assumptions by character block):
    if (any(!is.na(block_names))) {
      
      # Get numbers of blocks with labels:
      LabelledBlockNumbers <- which(!is.na(block_names))
      
      # For each labelled block:
      for(i in LabelledBlockNumbers) {
        
        # Build current label text (to look for in
        CurrentLabelText <- paste("(CHARACTERS=", block_names[i], ")=", sep = "")
        
        # Find line that contains label (if used):
        LabelMatch <- grep(CurrentLabelText, TypesetLines, fixed = TRUE)
        
        # If such a line was found:
        if (length(LabelMatch) > 0) {
          
          # Isolate ordering information:
          orderingInformation <- strsplit(TypesetLines[LabelMatch], split = CurrentLabelText, fixed = TRUE)[[1]][2]
          
          # Extract ordering information:
          orderingExtracted <- extract_assumptions(orderingInformation)
          
          # Store ordering information for block in ordering of block (i.e., numbered from 1 in block not 1 in whole NEXUS file):
          ordering[[i]] <- orderingExtracted
        
        }
        
      }
      
    # If there are not any block names (or if some blocks lack names):
    } else {
      
      # Isolate ordering information:
      orderingInformation <- strsplit(TypesetLines, split = "UNTITLED=", fixed = TRUE)[[1]][2]
      
      # Extract ordering information:
      orderingExtracted <- extract_assumptions(orderingInformation)
      
      # Get lengths of ordering for each block:
      orderingLengths <- unlist(lapply(ordering, length))
      
      # For each block of the matrix list:
      for(i in 1:length(orderingLengths)) {
        
        # Store part of extracted ordering corresponding to ith block:
        ordering[[i]][1:orderingLengths[i]] <-  orderingExtracted[1:orderingLengths[i]]
        
        # Remove already transferred ordering from extracted ready for next block:
        orderingExtracted <- orderingExtracted[-(1:orderingLengths[i])]
        
      }
      
    }
  
  }
  
  # Check for characters without an ordering status and stop and warn user if found:
  if (any(is.na(unlist(ordering)))) stop(paste("The following characters have undefined ordering: ", paste(which(is.na(unlist(ordering))), collapse = ", "), ". Add their ordering to the assumptions block and try again.", sep =""))
  
  # Convert any "Squared" ordering to "cont" for continuous:
  ordering <- lapply(ordering, gsub, pattern = "Squared", replacement = "cont", ignore.case = TRUE)
  
  # Look for any non-standard ordering (i.e., not cont, ord, unord, or Step_X):
  Nonstandardordering <- setdiff(unique(unlist(ordering)), c("cont", "ord", "unord", names(step_matrices)))
  
  # If any non-standard ordering is found stop and warn user:
  if (length(Nonstandardordering) > 0) stop(paste("The following non-standard character ordering(s) were found: ", paste(Nonstandardordering, collapse = ", "), ". These should be one of type \"cont\", \"ord\", \"unord\", or step matrix.", sep = ""))
  
  # Find any WTSET lines:
  weightsetLines <- X[lapply(lapply(strsplit(X, split = ""), '[', 1:5), paste, collapse = "") == "WTSET"]
  
  # If there are weightset line(s):
  if (length(weightsetLines) > 0) {
    
    # Remove any double spaces found (fot easier strsplit later):
    while(length(grep("  ", weightsetLines))) weightsetLines <- gsub("  ", " ", weightsetLines)
    
    # If there are any block names (i.e., names that may be used to differentiate assumptions by character block):
    if (any(!is.na(block_names))) {
      
      # Get numbers of blocks with labels:
      LabelledBlockNumbers <- which(!is.na(block_names))
      
      # For each labelled block:
      for(i in LabelledBlockNumbers) {
        
        # Build current label text (to look for in
        CurrentLabelText <- paste("(CHARACTERS=", block_names[i], ")=", sep = "")
        
        # Find line that contains label (if used):
        LabelMatch <- grep(CurrentLabelText, weightsetLines, fixed = TRUE)
        
        # If such a line was found:
        if (length(LabelMatch) > 0) {
          
          # Isolate weights information:
          weightsInformation <- strsplit(weightsetLines[LabelMatch], split = CurrentLabelText, fixed = TRUE)[[1]][2]
          
          # Extract weights information:
          weightsExtracted <- extract_assumptions(weightsInformation)
          
          # Store weights information for block in weights of block (i.e., numbered from 1 in block not 1 in whole NEXUS file):
          weights[[i]] <- weightsExtracted
          
        }
        
      }
      
    # If there are not any block names (or if some blocks lack names):
    } else {
      
      # Isolate weights information:
      weightsInformation <- strsplit(weightsetLines, split = "UNTITLED=", fixed = TRUE)[[1]][2]
      
      # Extract weights information:
      weightsExtracted <- extract_assumptions(weightsInformation)
      
      # Get lengths of weights for each block:
      weightsLengths <- unlist(lapply(weights, length))
      
      # For each block of the matrix list:
      for(i in 1:length(weightsLengths)) {
        
        # Store part of extracted weights corresponding to ith block:
        weights[[i]][1:weightsLengths[i]] <-  weightsExtracted[1:weightsLengths[i]]
        
        # Remove already transferred weights from extracted ready for next block:
        weightsExtracted <- weightsExtracted[-(1:weightsLengths[i])]
        
      }
      
    }
    
  }
  
  # Convert weights to numeric:
  weights <- lapply(weights, as.numeric)
  
  # If equalising weights:
  if (Equaliseweights) {
    
    # Get starting weights by taking differences for each character (will take reciprocal later for true weight):
    Startingweights <- lapply(lapply(MinMaxMatrixList, apply, 1, diff), function(x) x)
    
    # If there are continuous characters:
    if (any(names(Startingweights) == "CONTINUOUS")) {
      
      # Get numbers of continuous blocks:
      ContinuousBlocks <- which(names(Startingweights) == "CONTINUOUS")
      
      # Set weights for continuous blocks to null:
      for(i in ContinuousBlocks) Startingweights[[i]] <- vector(mode = "numeric")
      
    }
    
    # If there are any unordered characters:
    if (any(unlist(ordering) == "unord")) {
      
      # Get numbers of blocks with unordered characters:
      BlocksWithUnorderedcharacters <- which(unlist(lapply(lapply(ordering, '==', "unord"), sum)) > 0)
      
      # Convert all unordered characters to weight one:
      for(i in BlocksWithUnorderedcharacters) Startingweights[[i]][which(ordering[[i]] == "unord")] <- 1
      
    }
    
    # If there are any step matrices specified:
    if (any(unlist(lapply(lapply(ordering, grep, pattern = "step_"), length)) > 0)) {
      
      # Get numbers of blocks with step matrix characters:
      BlocksWithStepMatrixcharacters <- which(unlist(lapply(lapply(ordering, grep, pattern = "step_"), length)) > 0)
      
      # for each block with step matrix characters:
      for(i in BlocksWithStepMatrixcharacters) {
        
        # Get step matrix characters:
        Stepcharacters <- grep("step_", ordering[[i]])
        
        # Set step matrix character as maximum possible value in step matrix (again, should be reciprocal, but that will happen later:
        for(j in Stepcharacters) Startingweights[[i]][j] <- max(as.numeric(step_matrices[[ordering[[i]][j]]]))
        
      }
      
    }
    
    # Take reciprocal of weights so they are actual weights:
    Startingweights <- lapply(Startingweights, function(x) 1 / x)
    
    # Get product weight (to multiply all weights by):
    ProductWeight <- prod(unique(unlist(lapply(Startingweights, function(x) round(1 / x)))))
    
    # Multiply weighting product through all starting weights:
    Startingweights <- lapply(Startingweights, function(x) x * ProductWeight)
    
    # Get factors of every weight currently applied:
    AllFactorsCombined <- sort(unlist(lapply(as.list(unique(unlist(Startingweights))), get_all_factors)))
    
    # Get largest common factor of all weights:
    LargestCommonFactor <- max(rle(AllFactorsCombined)$values[rle(AllFactorsCombined)$lengths == length(unique(unlist(Startingweights)))])
    
    # As long as the largest common factor is greater than 1:
    while(LargestCommonFactor > 1) {
      
      # Update starting weights by dividing through by current largest facto:
      Startingweights <- lapply(Startingweights, function(x) x / LargestCommonFactor)
      
      # Get factors of every weight currently applied:
      AllFactorsCombined <- sort(unlist(lapply(as.list(unique(unlist(Startingweights))), get_all_factors)))
      
      # Get largest common factor of all weights:
      LargestCommonFactor <- max(rle(AllFactorsCombined)$values[rle(AllFactorsCombined)$lengths == length(unique(unlist(Startingweights)))])
      
    }
    
    # If there are any constant characters:
    if (any(unlist(lapply(lapply(lapply(lapply(MinMaxMatrixList, apply, 1, diff), '==', 0), which), length)) > 0)) {
      
      # Get numbers of blocks with constant characters:
      BlocksWithConstantcharacters <- which(unlist(lapply(lapply(lapply(lapply(MinMaxMatrixList, apply, 1, diff), '==', 0), which), length)) > 0)
      
      # For each block with constant characters set weight to zero:
      for(i in BlocksWithConstantcharacters) Startingweights[[i]][which(apply(MinMaxMatrixList[[i]], 1, diff) == 0)] <- 0
      
    }
    
    # If there are continuous characters:
    if (any(names(Startingweights) == "CONTINUOUS")) {
      
      # Get numbers of continuous blocks:
      ContinuousBlocks <- which(names(Startingweights) == "CONTINUOUS")
      
      # Set weights for continuous blocks as original weights:
      for(i in ContinuousBlocks) Startingweights[[i]] <- weights[[i]]
      
    }
    
    # Update weights:
    weights <- Startingweights
    
  }
  
  # Create top list:
  TopList <- list(header = textlines, step_matrices = step_matrices)
  
  # Start to compile output starting with top list:
  Output <- list(topper = TopList)
  
  # For each matrix block:
  for(i in 1:length(MatrixBlockList)) {
    
    # Create sublist for character information:
    characters <- list(symbols = symbols[[i]], missing = missing[[i]], gap = gap[[i]])
    
    # Build list for current block:
    Block <- list(block_name = block_names[[i]], datatype = names(MatrixBlockList)[i], matrix = MatrixBlockList[[i]], ordering = ordering[[i]], weights = weights[[i]], minimum_values = MinMaxMatrixList[[i]][, "Min"], maximum_values = MinMaxMatrixList[[i]][, "Max"], characters = characters)
    
    # Store current block in output:
    Output[[(i + 1)]] <- Block
    
    # Add name to current matrix (using a number to differentiate them):
    names(Output)[(i + 1)] <- paste("matrix_", i, sep = "")
    
  }

  # Return output:
  return(invisible(Output))

}
