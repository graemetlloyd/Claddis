#' Reads in a morphological #NEXUS data file
#'
#' @description
#'
#' Reads in a morphological data file in #NEXUS format.
#'
#' @param file_name The file name or path of the #NEXUS file.
#' @param equalize_weights Optional that overrides the weights specified in the file to make all characters truly equally weighted.
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
#' # Create example matrix
#' example_matrix <- paste("#NEXUS", "", "BEGIN DATA;",
#'                         "\tDIMENSIONS  NTAX=5 NCHAR=5;",
#'                         "\tFORMAT SYMBOLS= \" 0 1 2\" MISSING=? GAP=- ;",
#'                         "MATRIX", "", "Taxon_1  010?0", "Taxon_2  021?0",
#'                         "Taxon_3  02111", "Taxon_4  011-1",
#'                         "Taxon_5  001-1", ";", "END;", "",
#'                         "BEGIN ASSUMPTIONS;",
#'                         "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;",
#'                         "\tTYPESET * UNTITLED  = unord: 1 3-5, ord: 2;",
#'                         "\tWTSET * UNTITLED  = 1: 2, 2: 1 3-5;",
#'                         "END;", sep = "\n")
#'
#' # Write example matrix to current working directory called
#' # "morphmatrix.nex":
#' cat(example_matrix, file = "morphmatrix.nex")
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
read_nexus_matrix <- function(file_name, equalize_weights = FALSE) {
  
  # ADD ABILITY TO READ CHARSET LINES
  # COULD BE MULTIPLE TYPESET OR WTSET LINES, NEED TO CHECK FOR THIS!
  # ADD ABILITY TO READ AND STORE CHARACTER STATES LABELS ETC.
  # ADD ABILITY TO DEAL WITH RANGES FOR CONTINUOUS DATA
  # ALLOW READING OF TREES IF FOUND AND STORING IN TOPPER BUT MAYBE MAKE THIS OPTIONAL WITH DEFAULT TO NOT
  
  # Line formatting function to be used in lapply below to deal with polymorphic characters:
  format_lines <- function(x, direction = "in") {
    
    # Split current line by character:
    current_string <- strsplit(x, split = "")[[1]]
    
    # Check for square brackets:
    if (any(current_string == "[")) stop("TNT-style square brackets, [], found in matrix. Replace with () or {} depending on whether a true polymorphism or uncertainty respectively.")
    
    # If polymorphisms (i.e., parentheses) are found:
    if (any(c(current_string == "{", current_string == "("))) {
      
      # While there are polymorphisms:
      while(any(c(current_string == "{", current_string == "("))) {
        
        # Find beginning of first polymorphism (will loop through first polymorphism until none are left):
        first_polymorphism_begins <- c(which(x = current_string == "{"), which(x = current_string == "("))[1]
        
        # If polymorphsim begins with a parenthesis then it must end with a parenthesis:
        if (current_string[first_polymorphism_begins] == "(") first_polymorphism_ends_with <- ")"
        
        # If polymorphsim begins with a curly brace then it must end with a curly brace:
        if (current_string[first_polymorphism_begins] == "{") first_polymorphism_ends_with <- "}"
        
        # Find end of first polymorphism:
        first_polymorphism_ends <- which(x = current_string == first_polymorphism_ends_with)[1]
        
        # If a true polymorphism (two or more states observed) and reading data in use ampersand to separate states:
        if (current_string[first_polymorphism_begins] == "(" && direction == "in") polymorphic_character <- paste(sort(x = current_string[(first_polymorphism_begins + 1):(first_polymorphism_ends - 1)]), collapse = "&")
        
        # If an uncertainty (two or more possible states) and reading data in use slash to separate states:
        if (current_string[first_polymorphism_begins] == "{" && direction == "in") polymorphic_character <- paste(sort(x = current_string[(first_polymorphism_begins + 1):(first_polymorphism_ends - 1)]), collapse = "/")
        
        # If a true polymorphism (two or more states observed) and exporting data out use ampersand to separate states:
        if (current_string[first_polymorphism_begins] == "(" && direction == "out") polymorphic_character <- paste("(", paste(sort(x = current_string[(first_polymorphism_begins + 1):(first_polymorphism_ends - 1)]), collapse = ""), ")", sep = "")
        
        # If an uncertainty (two or more possible states) and exporting data out use slash to separate states:
        if (current_string[first_polymorphism_begins] == "{" && direction == "out") polymorphic_character <- paste("{", paste(sort(x = current_string[(first_polymorphism_begins + 1):(first_polymorphism_ends - 1)]), collapse = ""), "}", sep = "")
        
        # Remove now redundant characters from current string:
        current_string <- current_string[-((first_polymorphism_begins + 1):first_polymorphism_ends)]
        
        # Store polymorphic character in string:
        current_string[first_polymorphism_begins] <- polymorphic_character
        
      }
      
    }
    
    # Return string with polymorphisms formatted to single characters:
    return(current_string)
    
  }
  
  # Subfunction to extract ordering information:
  extract_assumptions <- function(assumption_line) {
    
    # Reformat as a list:
    x <- lapply(X = lapply(X = as.list(x = gsub(pattern = ";", replacement = "", x = strsplit(assumption_line, split = ", ")[[1]])), strsplit, split = ": "), unlist)
    
    # Extract just the numbers:
    character_numbers <- lapply(X = lapply(X = lapply(X = x, '[', 2), strsplit, split = " "), unlist)
    
    # If there are hyphens (ranges) in the data:
    if (length(x = grep("-", character_numbers)) > 0) {
      
      # Subfunction for unpacking ranges denoted by hyphens:
      unpack_ranges <- function(character_numbers) {
        
        # Whilst hyphens remain in the data:
        while(length(x = grep("-", character_numbers)) > 0) {
          
          # Find first hyphen:
          first_hyphen <- grep("-", character_numbers)[1]
          
          # Unpack and add at end whilst removing first hyphen:
          character_numbers <- c(character_numbers[-first_hyphen], as.character(c(as.numeric(strsplit(character_numbers[first_hyphen], split = "-")[[1]])[1]:as.numeric(strsplit(character_numbers[first_hyphen], split = "-")[[1]])[2])))
          
        }
        
        # Return character numbers:
        sort(x = as.numeric(character_numbers))
      }
      
      # Apply unpacking function to data:
      character_numbers <- lapply(X = character_numbers, unpack_ranges)
      
    # If there are no hyphens:
    } else {
      
      # Simply convert to numbers and ensure they are sorted in increasing order:
      character_numbers <- lapply(X = lapply(X = character_numbers, as.numeric), sort)
      
    }
    
    # Add character type as names to numbers list:
    names(character_numbers) <- unlist(x = lapply(X = x, '[', 1))
    
    # Build empty vector to store output:
    output <- vector(mode = "character", length = max(unlist(x = character_numbers)))
    
    # Store character type in output:
    for(i in 1:length(x = character_numbers)) output[character_numbers[[i]]] <- names(character_numbers)[i]
    
    # Check all characters have some kind of ordering described:
    if (length(x = which(x = output == "")) > 0) stop(paste("The following characters have no specified ordering: ", paste(which(x = output == ""), collapse = ", "), ". Check ASSUMPTIONS block and edit appropriately.", sep = ""))
    
    # Return output:
    output
  }
  
  # Subfunction to find range of states for each character:
  find_ranges <- function(x) {
    
    # Convert each column of matrix to a list of numeric values:
    x <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = apply(x, 2, as.list), unlist), strsplit, split = "&|/"), unlist), as.numeric), sort)
    
    # Check for any empty columns (all NAs) and insert a dummy 0 if found:
    if (any(lapply(X = x, length) == 0)) x[which(x = lapply(X = x, length) == 0)] <- 0
    
    # Convert to a min-max matrix:
    x <- matrix(unlist(x = lapply(X = x, range)), ncol = 2, byrow = TRUE, dimnames = list(c(), c("min", "max")))
    
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
  Sys.setlocale('LC_ALL', 'C')

  # Read in NEXUS file as raw text:
  raw_nexus <- readLines(file_name, warn = FALSE)

  # Check that this is a #NEXUS file:
  if (length(x = grep("#NEXUS", raw_nexus, ignore.case = TRUE)) == 0) stop("This is not a #NEXUS file.")

  # Replace any lower case matrix with upper case:
  if (length(x = grep("MATRIX", raw_nexus, ignore.case = FALSE)) == 0) raw_nexus <- gsub(pattern = "matrix", replacement = "MATRIX", x = raw_nexus, ignore.case = FALSE)

  # Check that the #NEXUS file includes a matrix:
  if (length(x = grep("MATRIX", raw_nexus, ignore.case = TRUE)) == 0) stop("There is no matrix present.")

  # Are there blocks of the NEXUS file that can be removed? (i.e., trees, macclade, mesquite etc.:
  if (length(x = grep("begin trees|begin mesquite|begin mrbayes|begin macclade|begin treebase|begin notes|begin paup", raw_nexus, ignore.case = TRUE)) > 0) {

    # Get position of all Begin tags:
    all_begins <- grep("begin ", raw_nexus, ignore.case = TRUE)

    # Get position of all End tags:
    all_ends <- grep("end;|endblock;|end ;|endblock ;", raw_nexus, ignore.case = TRUE)

    # Check Begin and End are in balance:
    if (length(x = all_begins) != length(x = all_ends)) stop("Begin and End tags are not equal in number.")

    # Get rows for blocks to delete:
    block_rows <- match(grep("begin trees|begin mesquite|begin mrbayes|begin macclade|begin treebase|begin notes|begin paup", raw_nexus, ignore.case = TRUE), all_begins)

    # Find start lines:
    block_begins <- all_begins[block_rows]

    # Find end lines:
    block_ends <- all_ends[block_rows]

    # Incrementally remove superflouous blocks:
    for(i in length(x = block_begins):1) raw_nexus <- raw_nexus[c(c(1:(block_begins[i] - 1)), c((block_ends[i] + 1):length(x = raw_nexus)))]

  }

  # Case if there are superflous taxon, character, or state labels:
  if (length(x = grep("taxlabels|charstatelabels|charlabels|statelabels", raw_nexus, ignore.case = TRUE)) > 0) {

    # Find label beginning rows:
    label_ends <- label_begins <- grep("taxlabels|charstatelabels|charlabels|statelabels", raw_nexus, ignore.case = TRUE)

    # For each label find the end tag that corresponds to it:
    for(i in 1:length(x = label_begins)) label_ends[i] <- min(grep(";", raw_nexus, ignore.case = TRUE)[grep(";", raw_nexus, ignore.case = TRUE) > label_begins[i]])

    # Incrementally remove superflouous labels:
    for(i in length(x = label_begins):1) raw_nexus <- raw_nexus[c(c(1:(label_begins[i] - 1)), c((label_ends[i] + 1):length(x = raw_nexus)))]

  }

  # If there is text inside single quotes:
  if (length(x = grep("'", raw_nexus)) > 0) {

    # Find items to replace (text not in single quotes):
    replacement_items <- unique(x = unlist(x = strsplit(gsub(pattern = paste("'[A-Z|\t| |0-9|a-z|\\?|\\-|\\(|)|\\.]*'", sep = ""), replacement = "''", x = raw_nexus[grep("'", raw_nexus)]), "''")))

    # Make sure nothing "" is not in this list:
    replacement_items <- setdiff(x = unique(x = unlist(x = strsplit(replacement_items, " "))), y = "")

    # Make block that isolates lines with single quotes:
    lines_to_edit <- raw_nexus[grep("'", raw_nexus)]

    # Remove text outside single quotes:
    for(i in replacement_items[order(nchar(x = replacement_items), decreasing = TRUE)]) lines_to_edit <- gsub(pattern = i, replacement = "", x = lines_to_edit, fixed = TRUE)

    # Remove double spaces:
    while(length(x = grep("  ", lines_to_edit)) > 0) lines_to_edit <- gsub(pattern = "  ", replacement = " ", x = lines_to_edit, fixed = TRUE)

    # Now isolate names within single quotes:
    names_in_single_quotes <- unique(x = gsub(pattern = "'", replacement = "", x = strsplit(trim_marginal_whitespace(paste(lines_to_edit, collapse = "")), "' '")[[1]]))

    # Make sure nothing "" is not in this list:
    names_in_single_quotes <- which(x = !names_in_single_quotes == "")

    # Replace names in single quotes with underscored versions without single quotes or parentheses:
    for(i in names_in_single_quotes) raw_nexus <- gsub(pattern = i, replacement = gsub(pattern = "\\(|)", replacement = "", x = gsub(pattern = " ", replacement = "_", x = i)), x = raw_nexus, fixed = TRUE)

    # Remove all single quotes from text file:
    raw_nexus <- gsub(pattern = "'", replacement = "", x = raw_nexus, fixed = TRUE)

  }

  # Get rid of spaces around equals signs to make later regular expresssions work:
  while(length(x = grep(" = |= | =", raw_nexus))) raw_nexus <- gsub(pattern = " = |= | =", replacement = "=", x = raw_nexus)

  # Remove weird characters (may need to add to the list or not if Sys.setlocale fixes the issues):
  raw_nexus <- gsub(pattern = "\x94|\x93|\xd5|\xd4|\xd3|\xd2|'", replacement = "", x = raw_nexus)
  
  # Replace all tabs with spaces:
  raw_nexus <- gsub(pattern = "\t", replacement = " ", x = raw_nexus)

  # Replace tabs with spaces and trim leading and trailing spaces from each line:
  raw_nexus <- unlist(x = lapply(X = as.list(raw_nexus), FUN = trim_marginal_whitespace))

  # Delete any empty lines (if present):
  if (length(x = which(x = raw_nexus == "")) > 0) raw_nexus <- raw_nexus[-which(x = raw_nexus == "")]
  
  # If header text is present:
  if (length(x = grep("\\[", raw_nexus)) > 0) {

    # Work out beginning and ending lines for text:
    text_lines <- apply(cbind(setdiff(x = grep("\\[", raw_nexus), grep("\\[[0-9:A-Z:a-z]{1}\\]", raw_nexus)), y = setdiff(x = grep("\\]", raw_nexus), y = grep("\\[[0-9:A-Z:a-z]{1}\\]", raw_nexus))), 1, paste, collapse = ":")

    # Convert beginning and endings to numerics for finding in vector:
    lines_to_delete <- text_lines <- eval(parse(text = paste("c(", paste(text_lines, collapse = ","), ")", sep = "")))
    
    # Find only lines beginning with a square bracket to avoid issue with polymorphic use of square brackets (can be caught later):
    lines_to_keep <- which(x = unlist(x = lapply(X = strsplit(raw_nexus[text_lines], split = ""), '[', 1)) == "[")
    
    # Reduce text_lines appropriately:
    text_lines <- text_lines[lines_to_keep]
    
    # Reduce lines to delete appropriately:
    lines_to_delete <- lines_to_delete[lines_to_keep]
  
    # Grab text and store:
    text_lines <- paste(gsub(pattern = "\\[|\\]", replacement = "", x = raw_nexus[text_lines]), collapse = "\n")
  
    # Delete text lines to avoid regular expression errors later:
    raw_nexus <- raw_nexus[-lines_to_delete]

  # If header text is absent:
  } else {
    
    # Store empty text string:
    text_lines <- ""

  }

  # Little test for valid NEXUS formatting for number of taxa:
  if (length(x = grep("ntax", raw_nexus, ignore.case = TRUE)) == 0) stop("Number of taxa not defined.")

  # Little test for valid NEXUS formatting for number of characters:
  if (length(x = grep("nchar", raw_nexus, ignore.case = TRUE)) == 0) stop("Number of characters not defined.")
  
  # Get number of taxa:
  n_taxa <- as.numeric(strsplit(strsplit(gsub(pattern = ";|\t", replacement = "", x = raw_nexus[grep("NTAX=|Ntax=|ntax=", raw_nexus)]), "NTAX=|Ntax=|ntax=")[[1]][2], " ")[[1]][1])
  
  # Create empty list to store matrix block(s):
  matrix_block_list <- list()

  # If there are character blocks (rather than a data block):
  if (length(x = grep("begin characters", raw_nexus, ignore.case = TRUE)) > 0) {
    
    # Find line where character block(s) begin:
    character_block_begins <- grep("begin characters", raw_nexus, ignore.case = TRUE)
    
    # Find line where matching character block(s) end:
    character_block_ends <- grep("end;", raw_nexus, ignore.case = TRUE)[unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = character_block_begins), '<', grep("end;", raw_nexus, ignore.case = TRUE)), which), min))]
    
    # For each character block:
    for(i in 1:length(x = character_block_begins)) {
      
      # isolate ith character block:
      current_character_block <- raw_nexus[character_block_begins[i]:character_block_ends[i]]
      
      # Locate dataype line:
      datatype_line <- grep("datatype=", current_character_block, ignore.case = TRUE)
      
      # Check there is a dataype:
      if (length(x = datatype_line) == 0) stop("Character datatype not specified.")
      
      # Check there are not multiple datatypes for a single character block:
      if (length(x = datatype_line) > 1) stop("Multiple character datatypes specified.")
      
      # Get ith datatype:
      current_datatype <- strsplit(strsplit(current_character_block[datatype_line], "datatype=|DATATYPE=")[[1]][2], " ")[[c(1,1)]]
      
      # Add ith character block to list:
      matrix_block_list[[i]] <- current_character_block
      
      # Add datatype to names for ith matrix block:
      names(matrix_block_list)[i] <- current_datatype

    }
    
  }
  
  # If there is a data block:
  if (length(x = grep("begin data", raw_nexus, ignore.case = TRUE)) > 0) {
    
    # Get beginning line of data block:
    begin_data_line <- grep("begin data", raw_nexus, ignore.case = TRUE)
    
    # Get end line of data block:
    end_data_line <- grep("end;", raw_nexus, ignore.case = TRUE)[unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = begin_data_line), '<', grep("end;", raw_nexus, ignore.case = TRUE)), which), min))]
    
    # Add data block to list:
    matrix_block_list[[1]] <- raw_nexus[begin_data_line:end_data_line]
    
    # Locate datatype line:
    datatype_line <- grep("datatype=", matrix_block_list[[1]], ignore.case = TRUE)
    
    # If there is a stated datatype:
    if (length(x = datatype_line) > 0) {
      
      # Store stated datatype as list name:
      names(matrix_block_list)[1] <- strsplit(strsplit(matrix_block_list[[1]][datatype_line], "datatype=|DATATYPE=")[[1]][2], " ")[[c(1,1)]]
      
    # If there is no stated datatype:
    } else {
      
      # Set datatype as "STANDARD":
      names(matrix_block_list)[1] <- "STANDARD"
      
    }

  }
  
  # If there are any interleaved block(s):
  if (length(x = unlist(x = lapply(X = matrix_block_list, grep, pattern = "INTERLEAVE"))) > 0) {
    
    # Get interleaved blocks:
    interleaved_blocks <- which(x = unlist(x = lapply(X = matrix_block_list, grep, pattern = "INTERLEAVE") > 0))
    
    # For each interleaved block:
    for(i in interleaved_blocks) {
      
      # Isolate current matrix block:
      current_block <- matrix_block_list[[i]]
      
      # Remove interleave text:
      current_block <- gsub(pattern = "INTERLEAVE ", replacement = "", x = current_block)
      
      # Find line matrix begins on:
      matrix_begins <- which(x = toupper(current_block) == "MATRIX") + 1
      
      # Find line matrix ends:
      matrix_ends <- rev(which(x = current_block == ";")) - 1
      
      # Get number of blocks in interleaved matrix:
      n_blocks <- length(x = matrix_begins:matrix_ends) / n_taxa
      
      # Check is divisible by number of taxa:
      if ((n_blocks %% 1) != 0) stop("Interleaved matrix has greater or fewer lines than a multiple of number of taxa.")
      
      # Empty lst to store interleaved blocks:
      interleaved_block_list <- list()
      
      # For each interleaved block:
      for(j in 1:n_blocks) {
        
        # Get end of block:
        block_ends <- j * n_taxa
        
        # Get beginning of block:
        block_begins <- block_ends - n_taxa + 1
        
        # Isolate block:
        block_j <- current_block[matrix_begins:matrix_ends][block_begins:block_ends]
        
        # Remove any double spaces:
        while(length(x = grep("  ", block_j)) > 0) block_j <- gsub(pattern = "  ", replacement = " ", x = block_j)
        
        # Check for rogue spaces and stop and warn if found:
        if (any(unlist(x = lapply(X = strsplit(block_j, split = " "), length)) > 2)) stop("Rogue spaces found in interleaved block. Check taxon names and polymorphisms and remove.")
        
        # Store jth block in interleaved blocks:
        interleaved_block_list[[j]] <- block_j
        
      }
      
      # Reformat as matrices:
      interleaved_block_list <- lapply(X = lapply(X = lapply(X = interleaved_block_list, strsplit, split = " "), unlist), matrix, ncol = 2, byrow = TRUE)
      
      # Combine as single large matrix:
      interleaved_block_list <- matrix(unlist(x = lapply(X = interleaved_block_list, '[', , 2, drop = FALSE)), nrow = n_taxa, byrow = FALSE, dimnames = list(interleaved_block_list[[1]][, 1], c()))
      
      # Collapse matrix to single lines:
      interleaved_block_list <- apply(cbind(cbind(rownames(x = interleaved_block_list), rep(" ", n_taxa)), interleaved_block_list), 1, paste, collapse = "")
      
      # Overwrite matrix block with now de-interleaved version:
      matrix_block_list[[i]] <- c(matrix_block_list[[i]][1:(matrix_begins - 1)], interleaved_block_list, matrix_block_list[[i]][(matrix_ends + 1):length(x = matrix_block_list[[i]])])
      
    }
    
  }
  
  # If any matrix blocks are of type "MIXED":
  if (any(lapply(X = lapply(X = lapply(X = strsplit(names(matrix_block_list), split = ""), '[', 1:5), paste, collapse = ""), "toupper") == "MIXED")) {
    
    # Get any mixed datatype block(s):
    mixed_blocks <- which(x = lapply(X = lapply(X = lapply(X = strsplit(names(matrix_block_list), split = ""), '[', 1:5), paste, collapse = ""), "toupper") == "MIXED")
    
    # For each mixed block:
    for(i in mixed_blocks) {
      
      # Isolate current block:
      current_block <- matrix_block_list[[i]]
      
      # Reduce to MATRIX block:
      current_block <- current_block[(which(x = toupper(current_block) == "MATRIX") + 1):(rev(which(x = current_block == ";")) - 1)]
      
      # Remove any double spaces:
      while(length(x = grep("  ", current_block))) current_block <- gsub(pattern = "  ", replacement = " ", x = current_block)
      
      # Check for any rogue spaces:
      if (any(unlist(x = lapply(X = strsplit(current_block, split = " "), length)) > 2)) stop("Rogue spaces found in matrix block. Check taxon names and polymorphisms.")
      
      # Isolate taxon names:
      taxon_names <- matrix(unlist(x = strsplit(current_block, split = " ")), ncol = 2, nrow = n_taxa, byrow = TRUE)[, 1]
      
      # Reformat data block as list:
      current_block <- as.list(x = matrix(unlist(x = strsplit(current_block, split = " ")), ncol = 2, nrow = n_taxa, byrow = TRUE)[, 2])
      
      # Convert data block into individual character vectors:
      current_block <- lapply(X = current_block, format_lines, direction = "out")
      
      # Add taxon names to list:
      names(current_block) <- taxon_names
      
      # Isolate components of mixed block:
      mixed_components <- strsplit(names(matrix_block_list)[i], "\\(|\\)")[[1]][2]
      
      # Isolate further (remove spaces so components can be isolated):
      mixed_components <- strsplit(gsub(pattern = " ", replacement = "", x = mixed_components), split = ",")[[1]]
      
      # Format as list:
      mixed_components <- lapply(X = strsplit(mixed_components, split = ":"), toupper)
      
      # Get total number of characters in block:
      total_characters <- max(as.numeric(unlist(x = lapply(X = lapply(X = mixed_components, '[', 2), strsplit, split = "-")))) - min(as.numeric(unlist(x = lapply(X = lapply(X = mixed_components, '[', 2), strsplit, split = "-")))) + 1
      
      # Check all rows have the same (correct) number of characters:
      if (any(!unlist(x = lapply(X = current_block, length)) == total_characters)) stop("Some lines of matrix have too many or too few characters.")
      
      # For each component:
      for(j in 1:length(x = mixed_components)) {
        
        # Isolate current partition as matrix:
        current_matrix <- matrix(unlist(x = lapply(X = current_block, '[', as.numeric(strsplit(mixed_components[[j]][2], split = "-")[[1]])[1]:as.numeric(strsplit(mixed_components[[j]][2], split = "-")[[1]])[2])), nrow = n_taxa, byrow = TRUE, dimnames = list(names(current_block), c()))
        
        # Add component of matrix to end of list:
        matrix_block_list[[(length(x = matrix_block_list) + 1)]] <- c(toupper(matrix_block_list[[i]][1:which(x = toupper(matrix_block_list[[i]]) == "MATRIX")]), apply(cbind(cbind(rownames(x = current_matrix), rep(" ", n_taxa)), current_matrix), 1, paste, collapse = ""), toupper(matrix_block_list[[i]][rev(which(x = matrix_block_list[[i]] == ";"))[1]:length(x = matrix_block_list[[i]])]))
        
        # Remove names:
        names(matrix_block_list[[length(x = matrix_block_list)]]) <- NULL
        
        # Add datatype name as name to new list item:
        names(matrix_block_list)[length(x = matrix_block_list)] <- mixed_components[[j]][1]
        
        # Get number of characters in partition:
        characters_in_partition <- max(as.numeric(strsplit(mixed_components[[j]][2], split = "-")[[1]])) - min(as.numeric(strsplit(mixed_components[[j]][2], split = "-")[[1]])) + 1
        
        # Update number of characters for parition:
        matrix_block_list[[length(x = matrix_block_list)]] <- gsub(pattern = paste("NCHAR=", total_characters, sep = ""), replacement = paste("NCHAR=", characters_in_partition, sep = ""), x = matrix_block_list[[length(x = matrix_block_list)]], fixed = TRUE)
        
        # Update datatype for partition:
        matrix_block_list[[length(x = matrix_block_list)]] <- gsub(pattern = paste("DATATYPE=", toupper(names(matrix_block_list)[i]), sep = ""), replacement = paste("DATATYPE=", mixed_components[[j]][1], sep = ""), x = matrix_block_list[[length(x = matrix_block_list)]], fixed = TRUE)
        
      }
      
    }
    
    # Remove now redundant mixed blocks-:
    matrix_block_list <- matrix_block_list[-mixed_blocks]

  }
  
  # Will need to know datatypes before proceeding so check for nonstandard ones now:
  nonstandard_datatypes <- setdiff(x = toupper(names(matrix_block_list)), y = c("CONTINUOUS", "DNA", "NUCLEOTIDE", "PROTEIN", "RESTRICTION", "RNA", "STANDARD"))
  
  # Stop if non-standard datatype(s) found:
  if (length(x = nonstandard_datatypes) > 0) stop("Non-standard datatype found (i.e., neither CONTINUOUS, DNA, NUCLEOTIDE, PROTEIN, RESTRICTION, RNA, or STANDARD).")
  
  # Create empty block names vector:
  block_names <- rep(NA, length(x = matrix_block_list))
  
  # For each matrix block:
  for(i in 1:length(x = matrix_block_list)) {
    
    # Isolate ith block:
    current_block <- matrix_block_list[[i]]
    
    # Get title line (if found) - should work even if there is a taxon whose name begins "TITLE":
    title_line <- intersect(which(x = unlist(x = lapply(X = lapply(X = strsplit(current_block, split = ""), '[', 1:5), paste, collapse = "")) == "TITLE"), 1:which(x = unlist(x = lapply(X = lapply(X = strsplit(current_block, split = ""), '[', 1:5), paste, collapse = "")) == "MATRI")[1])
    
    # Store block name:
    if (length(x = title_line) > 0) block_names[i] <- gsub(pattern = ";", replacement = "", x = rev(strsplit(current_block[title_line], split = " ")[[1]])[1])
    
  }

  # Get number of characters in each block:
  n_characters <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = matrix_block_list, function(x) x[grep("NCHAR=|Nchar=|nchar=", x)]), gsub, pattern = ";|\t", replacement = ""), strsplit, split = "NCHAR=|Nchar=|nchar="), '[[', 1), '[', 2), strsplit, split = " "), '[[', 1), '[', 1), as.numeric)
  
  # Check matri(ces) have dimension:
  if (!(all(n_characters > 0) && n_taxa > 0)) stop("One or more matrix blocks have no dimensions (i.e., zero characters or taxa).")
  
  # Function to get symbols from matrix block:
  get_symbols <- function(x) {
    
    # If symbols are specified in the file:
    if (length(x = grep("symbols", x, ignore.case = TRUE)) > 0) {
      
      # Get initial set of symbols:
      symbols <- strsplit(strsplit(strsplit(x[grep("symbols", x, ignore.case = TRUE)], "SYMBOLS=|symbols=|symbols=")[[1]][2], "\"")[[1]][2], " |")[[1]]
      
      # Collapse to just the symbols themselves:
      symbols <- symbols[nchar(x = symbols) == 1]
      
      # Special case of a tilda indicating a range:
      if (length(x = symbols) == 3 && symbols[2] == "~") symbols <- as.character(as.numeric(symbols[1]):as.numeric(symbols[3]))
      
    # If symbols are not specified in the file:
    } else {
      
      # Find datatype row (if specified):
      datatype_row <- grep("datatype", x, ignore.case = TRUE)
      
      # If a datatype was specified:
      if (length(x = datatype_row) > 0) {
        
        # Store datatype:
        datatype <- strsplit(strsplit(toupper(x[datatype_row]), split = "DATATYPE=")[[1]][2], split = " ")[[1]][1]
        
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
  symbols <- lapply(X = matrix_block_list, get_symbols)
  
  # Get missing character function:
  get_missing <- function(x) {
    
    # If the missing character is specified:
    if (length(x = grep("missing", x, ignore.case = TRUE)) > 0) {
      
      # Get missing character:
      missing <- strsplit(strsplit(x[grep("missing", x, ignore.case = TRUE)][1], "MISSING=|missing=|missing=")[[1]][2], "")[[1]][1]
      
    # If the missing character is not specified:
    } else {
      
      # Set as default symbol of the question mark:
      missing <- "?"
      
    }
    
    # Return gap character:
    return(missing)

  }
  
  # Get missing symbol(s):
  missing <- lapply(X = matrix_block_list, get_missing)
  
  # Get gap character function:
  get_gap <- function(x) {
    
    # If the gap character is specified:
    if (length(x = grep("gap", x[1:grep("MATRIX", x)[1]], ignore.case = TRUE)) > 0) {
      
      # Get gap character:
      gap <- strsplit(strsplit(x[grep("gap", x, ignore.case = TRUE)][1], "GAP=|gap=|gap=")[[1]][2], "")[[1]][1]
      
    # If the gap character is not specified:
    } else {
      
      # Set as default symbol of the dash:
      gap <- "-"
      
    }
    
    # Return gap character:
    return(gap)
    
  }
  
  # Get gap symbol(s):
  gap <- lapply(X = matrix_block_list, get_gap)
  
  # Find first line of each matrix:
  matrix_start_lines <- lapply(X = lapply(X = lapply(X = matrix_block_list, '==', "MATRIX"), which), '+', 1)
  
  # Find alst line of each matrix:
  matrix_end_lines <- lapply(X = lapply(X = lapply(X = matrix_block_list, '==', ";"), which), '-', 1)
  
  # Cut down matrices to just matrix block:
  for(i in 1:length(x = matrix_block_list)) matrix_block_list[[i]] <- matrix_block_list[[i]][matrix_start_lines[[i]]:matrix_end_lines[[i]]]
  
  # Stop if incorrect number of taxa found:
  if (any(lapply(X = matrix_block_list, length) != n_taxa)) stop("Some matrix block(s) have too many or too few taxa. Check for linebreaks or incorrect NTAX value.")
  
  # If any blocks are of type "CONTINUOUS":
  if (any(toupper(names(matrix_block_list)) == "CONTINUOUS")) {
    
    # Get block(s) of type "CONTINUOUS":
    continuous_blocks <- which(x = toupper(names(matrix_block_list)) == "CONTINUOUS")
    
    # For each block of type "CONTINUOUS":
    for(i in continuous_blocks) {
      
      # Isolate ith continuous matrix block:
      current_block <- matrix_block_list[[i]]
      
      # While there are double spaces (that will complicate strsplit later) replace these with single spaces:
      while(length(x = grep("  ", current_block)) > 0) current_block <- gsub(pattern = "  ", replacement = " ", x = current_block)
      
      # Check all rows are equal in length and stop and warn if not:
      if (length(x = unique(x = unlist(x = lapply(X = strsplit(current_block, split = " "), length)))) > 1) stop("Continuous block has missing codings for some taxa (check for missing spaces or spaces in taxon names).")
      
      # Format data as matrix with rownames as taxa:
      current_block <- matrix(unlist(x = strsplit(current_block, split = " ")), nrow = n_taxa, byrow = TRUE, dimnames = list(unlist(x = lapply(X = strsplit(current_block, split = " "), '[', 1)), c()))
      
      # Remove names (first column) from matrix:
      current_block <- current_block[, -1, drop = FALSE]
      
      # Replace missing values with NA (if there are any):
      if (length(x = grep(missing[[continuous_blocks]], current_block)) > 0) current_block <- gsub(pattern = missing[[continuous_blocks]], replacement = NA, x = current_block, fixed = TRUE)
      
      # Store newly formatted matrix back into matrix block list:
      matrix_block_list[[i]] <- current_block
      
    }
    
  }
  
  # If there are non-continuous blocks:
  if (length(x = setdiff(x = toupper(names(matrix_block_list)), y = "CONTINUOUS")) > 0) {
    
    # Get non-continuous block number(s):
    noncontinuous_blocks <- which(x = !toupper(names(matrix_block_list)) == "CONTINUOUS")
    
    # For each non-continuous block:
    for(i in noncontinuous_blocks) {
      
      # Isolate ith non-continuous matrix block:
      current_block <- matrix_block_list[[i]]
      
      # While there are double spaces (that will complicate strsplit later) replace these with single spaces:
      while(length(x = grep("  ", current_block)) > 0) current_block <- gsub(pattern = "  ", replacement = " ", x = current_block)
      
      # Check for rogue spaces (should just be one between taxon name and actual data (i.e., none in taxon names or polymorphisms):
      if (any(unlist(x = lapply(X = strsplit(current_block, split = " "), length)) > 2)) stop("Rogue spaces found (check taxon names and polymorphisms and remove before trying again).")
      
      # Reformat as one-column matrix with taxa as row names:
      current_block <- matrix(unlist(x = strsplit(current_block, split = " ")), ncol = 2, byrow = TRUE, dimnames = list(lapply(X = strsplit(current_block, split = " "), '[', 1), c()))[, 2, drop = FALSE]
      
      # Overwrite curentmatrix block with newly formatted matrix:
      matrix_block_list[[i]] <- current_block
      
      # Format current block as a list (loses taxon names for now hence storage line above):
      current_block <- as.list(x = current_block)
      
      # Apply line formatting function to all rows:
      current_block <- lapply(X = current_block, format_lines, direction = "in")
      
      # If any rows have too many or too few characters:
      if (any(lapply(X = current_block, length) != n_characters[[i]])) {
        
        # Isolate problem rows:
        problem_rows <- which(x = unlist(x = lapply(X = current_block, length)) != n_characters[[i]])
        
        # Stop and warn user:
        stop(paste("The following tax(a) have too many or too few characters: ", paste(rownames(x = matrix_block_list[[i]])[problem_rows], collapse = ", "), ". Check rows or NCHAR.", sep = ""))
        
      }
      
      # Store formatted matrix in matrix block list:
      matrix_block_list[[i]] <- matrix(unlist(x = current_block), nrow = n_taxa, ncol = n_characters[[i]], byrow = TRUE, dimnames = list(rownames(x = matrix_block_list[[i]]), c()))
      
      # Find any unexpected (rogue) characters in matrix:
      rogue_characters <- setdiff(x = unique(x = unlist(x = strsplit(unique(x = as.vector(matrix_block_list[[i]])), split = "&|/"))), y = c(gap[[i]], missing[[i]], symbols[[i]]))
      
      # If any rogue characters were found stop and warn user:
      if (length(x = rogue_characters) > 0) stop(paste("The following rogue characters were found: ", paste(rogue_characters, collapse = ", "), ". Consider adding these to SYMBOLS or check if they were intended to be replaced with polymorphisms in the source matrix.", sep = ""))
      
      # Replace missing values with NA (if there are any):
      if (length(x = grep(missing[[i]], matrix_block_list[[i]])) > 0) matrix_block_list[[i]] <- gsub(pattern = missing[[i]], replacement = NA, x = matrix_block_list[[i]], fixed = TRUE)

    }
    
  }
  
  # Get rownames from each matrix block (to check for any issues with these now they are isolated):
  row_names <- lapply(X = matrix_block_list, rownames)
  
  # Check for cuplicate taxon names and stop and warn if found:
  if (any(unlist(x = lapply(X = row_names, duplicated)))) stop(paste("The following taxon name(s) are duplicated:", paste(unlist(x = row_names)[which(x = unlist(x = lapply(X = row_names, duplicated)))], collapse = ", "), ". Rename so that all taxon names are unique.", sep = ""))
  
  # Check names match across matrix blocks and stop and warn if not:
  if (length(x = unlist(x = lapply(X = row_names, setdiff, unique(x = unlist(x = row_names))))) > 0) stop(paste("Taxa do no match across matrx blocks. Check the following names:", paste(unlist(x = lapply(X = row_names, setdiff, unique(x = unlist(x = row_names)))), collapse = ", "), ".", sep = ""))
  
  # Store taxon names:
  taxon_names <- unique(x = unlist(x = row_names))
  
  # Find and special (non-underscore or alphanumeric) characters in taxon names:
  special_characters_in_names <- unique(x = strsplit(paste(gsub(pattern = "[:A-Z:a-z:0-9:_:]", replacement = "", x = taxon_names), collapse = ""), split = "")[[1]])
  
  # If any special characters are found stop and warn user:
  if (length(x = special_characters_in_names) > 0) stop(paste("The following special characters were found in the taxon names: ", paste(special_characters_in_names, collapse = ", "), ". Remove or replace and try again.", sep = ""))
  
  # Get rid of spaces around dashes to make later regular expresssions work:
  while(length(x = grep(" - |- | -", raw_nexus))) raw_nexus <- gsub(pattern = " - |- | -", replacement = "-", x = raw_nexus)
  
  # Create null variable for step matrices:
  step_matrices <- NULL
  
  # Special case if there are user-defined characters, i.e. step matrices:
  if (length(x = grep("USERTYPE", raw_nexus, ignore.case = TRUE)) > 0) {
    
    # Get rows corresponding to start of stepmatrices:
    step_matrix_rows <- grep("USERTYPE", raw_nexus, ignore.case = TRUE)
    
    # Create empty list to store step matrices:
    step_matrices <- list()
    
    # Little check in case of abnormally large number of step matrices:
    if (length(x = step_matrix_rows) > length(x = LETTERS)) stop("Not enough letters to store the number of different step matrices.")
    
    # For each step matrix:
    for(i in step_matrix_rows) {
      
      # Grab text block corresponding to step matrix:
      step_matrix_block <- raw_nexus[i:(i + as.numeric(strsplit(raw_nexus[i], "\\(STEPMATRIX\\)=")[[1]][2]) + 1)]
      
      # Remove [] labels for rows:
      step_matrix_block <- gsub(pattern = "\\[[:A-Z:]\\] |\\[[:0-9:]\\] ", replacement = "", x = step_matrix_block)
      
      # Get the step matrix as a matrix (might still need to think about how to define e.g. infinity):
      step_matrix <- gsub(pattern = "\\.", replacement = "0", x = matrix(unlist(x = strsplit(step_matrix_block[3:length(x = step_matrix_block)], " ")), ncol = as.numeric(strsplit(raw_nexus[i], "\\(STEPMATRIX\\)=")[[1]][2]), byrow = TRUE))
      
      # Add row and column names (assumes they start at zero and climb from there - i.e., should match symbols in order from 0 to N - 1):
      rownames(x = step_matrix) <- colnames(x = step_matrix) <- as.character(0:(ncol(step_matrix) - 1))
      
      # Add step matrix to list:
      step_matrices[[length(x = step_matrices) + 1]] <- step_matrix
      
      # Find step matrix name:
      step_matrix_name <- strsplit(strsplit(raw_nexus[i], "USERTYPE ")[[1]][2], " ")[[1]][1]
      
      # Replace step matrix name with step_N:
      raw_nexus <- gsub(pattern = step_matrix_name, replacement = paste("step_", LETTERS[length(x = step_matrices)], sep = ""), x = raw_nexus)
      
      # Use step matrix name (step_N) in list for later calling:
      names(step_matrices)[length(x = step_matrices)] <- paste("step_", LETTERS[length(x = step_matrices)], sep = "")
      
    }
    
  }
  
  # Set default weights as 1:
  character_weights <- lapply(X = lapply(X = matrix_block_list, ncol), rep, x = 1)
  
  # Set default ordering as unordered:
  ordering <- lapply(X = lapply(X = matrix_block_list, ncol), rep, x = "unord")
  
  # For each matrix block:
  for(i in 1:length(x = matrix_block_list)) {
    
    # For each symbol replace with a number from 0 to N - 1 symbols (unless continuous which would be NULL for symbols):
    if (!is.null(symbols[[i]][1])) for(j in 1:length(x = symbols[[i]])) matrix_block_list[[i]] <- gsub(pattern = symbols[[i]][j], replacement = as.character(j - 1), x = matrix_block_list[[i]], fixed = TRUE)
    
    # Convert gap character(s) into empty text strings (""):
    matrix_block_list[[i]] <- gsub(pattern = gap[[i]], replacement = "", matrix_block_list[[i]], fixed = TRUE)
    
  }
  
  # Get minimum and maximum values for each character in each matrix:
  min_max_matrix_list <- lapply(X = matrix_block_list, find_ranges)
  
  # Now min-max is known need to check for continuous characters to make sure default character_weights are all effectively 1:
  if (any(names(matrix_block_list) == "CONTINUOUS")) {
    
    # Get numbers of continuous blocks:
    continuous_blocks <- which(x = names(matrix_block_list) == "CONTINUOUS")
    
    # For each continuous blocks set weights as reciprocal of difference between min and max (i.e., effectively setting all weights as one):
    for(i in continuous_blocks) character_weights[[i]] <- as.numeric(gsub(pattern = Inf, replacement = 1, x = 1 / (min_max_matrix_list[[i]][, "max"] - min_max_matrix_list[[i]][, "min"])))
    
  }

  # Set default ordering as unord:
  default_ordering <- "unord"
  
  # If deafult ordering is specified, store it:
  if (length(x = grep("DEFTYPE", toupper(raw_nexus))) > 0) default_ordering <- strsplit(strsplit(raw_nexus[grep("deftype", raw_nexus, ignore.case = TRUE)], "DEFTYPE=|Deftype=|deftype=")[[1]][2], " ")[[1]][1]
  
  # If default ordering is ordered then update ordering as "ord":
  if (default_ordering == "ord") ordering <- lapply(X = lapply(X = matrix_block_list, ncol), rep, x = "ord")
  
  # If one or more blocks are of type CONTINUOUS:
  if (any(toupper(names(matrix_block_list)) == "CONTINUOUS")) {
    
    # Get continuous blocks:
    continuous_blocks <- which(x = toupper(names(matrix_block_list)) == "CONTINUOUS")
    
    # Update continuous ordering to "cont" to denote characters are tp be treated as continuous:
    for(i in continuous_blocks) ordering[[i]] <- gsub(pattern = "unord", replacement = "cont", x = ordering[[i]])
    
  }
  
  # Check assumptions block is supplied (to catch issue of character information appearing in a different block type (e.g., MRBAYES):
  if (length(x = grep("begin assumptions", raw_nexus, ignore.case = TRUE)) == 0) stop("No assumptions block specified. NB: Claddis can not read ordering information stored in another type of block (e.g., MRBAYES).")
  
  # Find any TYPESET lines:
  typeset_lines <- raw_nexus[lapply(X = lapply(X = strsplit(raw_nexus, split = ""), '[', 1:7), paste, collapse = "") == "TYPESET"]
  
  # If there are typeset line(s):
  if (length(x = typeset_lines) > 0) {
    
    # Remove any double spaces found (fot easier strsplit later):
    while(length(x = grep("  ", typeset_lines))) typeset_lines <- gsub(pattern = "  ", replacement = " ", x = typeset_lines)
    
    # If there are any block names (i.e., names that may be used to differentiate assumptions by character block):
    if (any(!is.na(block_names))) {
      
      # Get numbers of blocks with labels:
      labelled_block_numbers <- which(x = !is.na(block_names))
      
      # For each labelled block:
      for(i in labelled_block_numbers) {
        
        # Build current label text (to look for in
        current_label_text <- paste("(CHARACTERS=", block_names[i], ")=", sep = "")
        
        # Find line that contains label (if used):
        label_match <- grep(current_label_text, typeset_lines, fixed = TRUE)
        
        # If such a line was found:
        if (length(x = label_match) > 0) {
          
          # Isolate ordering information:
          ordering_information <- strsplit(typeset_lines[label_match], split = current_label_text, fixed = TRUE)[[1]][2]
          
          # Extract ordering information:
          ordering_extracted <- extract_assumptions(assumption_line = ordering_information)
          
          # Store ordering information for block in ordering of block (i.e., numbered from 1 in block not 1 in whole NEXUS file):
          ordering[[i]] <- ordering_extracted
        
        }
        
      }
      
    # If there are not any block names (or if some blocks lack names):
    } else {
      
      # Isolate ordering information:
      ordering_information <- strsplit(typeset_lines, split = "UNTITLED=", fixed = TRUE)[[1]][2]
      
      # Extract ordering information:
      ordering_extracted <- extract_assumptions(assumption_line = ordering_information)
      
      # Get lengths of ordering for each block:
      ordering_lengths <- unlist(x = lapply(X = ordering, length))
      
      # For each block of the matrix list:
      for(i in 1:length(x = ordering_lengths)) {
        
        # Store part of extracted ordering corresponding to ith block:
        ordering[[i]][1:ordering_lengths[i]] <-  ordering_extracted[1:ordering_lengths[i]]
        
        # Remove already transferred ordering from extracted ready for next block:
        ordering_extracted <- ordering_extracted[-(1:ordering_lengths[i])]
        
      }
      
    }
  
  }
  
  # Check for characters without an ordering status and stop and warn user if found:
  if (any(is.na(unlist(x = ordering)))) stop(paste("The following characters have undefined ordering: ", paste(which(x = is.na(unlist(x = ordering))), collapse = ", "), ". Add their ordering to the assumptions block and try again.", sep =""))
  
  # Convert any "Squared" ordering to "cont" for continuous:
  ordering <- lapply(X = ordering, gsub, pattern = "Squared", replacement = "cont", ignore.case = TRUE)
  
  # Look for any non-standard ordering (i.e., not cont, ord, unord, or Step_X):
  nonstandard_ordering <- setdiff(x = unique(x = unlist(x = ordering)), y = c("cont", "ord", "unord", names(step_matrices)))
  
  # If any non-standard ordering is found stop and warn user:
  if (length(x = nonstandard_ordering) > 0) stop(paste("The following non-standard character ordering(s) were found: ", paste(nonstandard_ordering, collapse = ", "), ". These should be one of type \"cont\", \"ord\", \"unord\", or step matrix.", sep = ""))
  
  # Find any WTSET lines:
  weightset_lines <- raw_nexus[lapply(X = lapply(X = strsplit(raw_nexus, split = ""), '[', 1:5), paste, collapse = "") == "WTSET"]
  
  # If there are weightset line(s):
  if (length(x = weightset_lines) > 0) {
    
    # Remove any double spaces found (fot easier strsplit later):
    while(length(x = grep("  ", weightset_lines))) weightset_lines <- gsub(pattern = "  ", replacement = " ", x = weightset_lines)
    
    # If there are any block names (i.e., names that may be used to differentiate assumptions by character block):
    if (any(!is.na(block_names))) {
      
      # Get numbers of blocks with labels:
      labelled_block_numbers <- which(x = !is.na(block_names))
      
      # For each labelled block:
      for(i in labelled_block_numbers) {
        
        # Build current label text (to look for in
        current_label_text <- paste("(CHARACTERS=", block_names[i], ")=", sep = "")
        
        # Find line that contains label (if used):
        label_match <- grep(current_label_text, weightset_lines, fixed = TRUE)
        
        # If such a line was found:
        if (length(x = label_match) > 0) {
          
          # Isolate weights information:
          weighting_information <- strsplit(weightset_lines[label_match], split = current_label_text, fixed = TRUE)[[1]][2]
          
          # Extract weights information:
          weighting_extracted <- extract_assumptions(assumption_line = weighting_information)
          
          # Store weights information for block in weights of block (i.e., numbered from 1 in block not 1 in whole NEXUS file):
          character_weights[[i]] <- weighting_extracted
          
        }
        
      }
      
    # If there are not any block names (or if some blocks lack names):
    } else {
      
      # Isolate weights information:
      weighting_information <- strsplit(weightset_lines, split = "UNTITLED=", fixed = TRUE)[[1]][2]
      
      # Extract weights information:
      weighting_extracted <- extract_assumptions(assumption_line = weighting_information)
      
      # Get lengths of weights for each block:
      weighting_lengths <- unlist(x = lapply(X = character_weights, length))
      
      # For each block of the matrix list:
      for(i in 1:length(x = weighting_lengths)) {
        
        # Store part of extracted weights corresponding to ith block:
        character_weights[[i]][1:weighting_lengths[i]] <-  weighting_extracted[1:weighting_lengths[i]]
        
        # Remove already transferred weights from extracted ready for next block:
        weighting_extracted <- weighting_extracted[-(1:weighting_lengths[i])]
        
      }
      
    }
    
  }
  
  # Convert weights to numeric:
  character_weights <- lapply(X = character_weights, as.numeric)
  
  # If equalising weights:
  if (equalize_weights) {
    
    # Get starting weights by taking differences for each character (will take reciprocal later for true weight):
    starting_weights <- lapply(X = lapply(X = min_max_matrix_list, apply, 1, diff), function(x) x)
    
    # If there are continuous characters:
    if (any(names(starting_weights) == "CONTINUOUS")) {
      
      # Get numbers of continuous blocks:
      continuous_blocks <- which(x = names(starting_weights) == "CONTINUOUS")
      
      # Set weights for continuous blocks to null:
      for(i in continuous_blocks) starting_weights[[i]] <- vector(mode = "numeric")
      
    }
    
    # If there are any unordered characters:
    if (any(unlist(x = ordering) == "unord")) {
      
      # Get numbers of blocks with unordered characters:
      blocks_with_unordered_characters <- which(x = unlist(x = lapply(X = lapply(X = ordering, '==', "unord"), sum)) > 0)
      
      # Convert all unordered characters to weight one:
      for(i in blocks_with_unordered_characters) starting_weights[[i]][which(x = ordering[[i]] == "unord")] <- 1
      
    }
    
    # If there are any step matrices specified:
    if (any(unlist(x = lapply(X = lapply(X = ordering, grep, pattern = "step_"), length)) > 0)) {
      
      # Get numbers of blocks with step matrix characters:
      blocks_with_stepmatrix_characters <- which(x = unlist(x = lapply(X = lapply(X = ordering, grep, pattern = "step_"), length)) > 0)
      
      # for each block with step matrix characters:
      for(i in blocks_with_stepmatrix_characters) {
        
        # Get step matrix characters:
        step_characters <- grep("step_", ordering[[i]])
        
        # Set step matrix character as maximum possible value in step matrix (again, should be reciprocal, but that will happen later:
        for(j in step_characters) starting_weights[[i]][j] <- max(as.numeric(step_matrices[[ordering[[i]][j]]]))
        
      }
      
    }
    
    # Take reciprocal of weights so they are actual weights:
    starting_weights <- lapply(X = starting_weights, function(x) 1 / x)
    
    # Get product weight (to multiply all weights by):
    product_weight <- prod(unique(x = unlist(x = lapply(X = starting_weights, function(x) round(1 / x)))))
    
    # Multiply weighting product through all starting weights:
    starting_weights <- lapply(X = starting_weights, function(x) x * product_weight)
    
    # Get factors of every weight currently applied:
    all_factors_combined <- sort(x = unlist(x = lapply(X = as.list(x = unique(x = unlist(x = starting_weights))), get_all_factors)))
    
    # Get largest common factor of all weights:
    largest_common_factor <- max(rle(all_factors_combined)$values[rle(all_factors_combined)$lengths == length(x = unique(x = unlist(x = starting_weights)))])
    
    # As long as the largest common factor is greater than 1:
    while(largest_common_factor > 1) {
      
      # Update starting weights by dividing through by current largest facto:
      starting_weights <- lapply(X = starting_weights, function(x) x / largest_common_factor)
      
      # Get factors of every weight currently applied:
      all_factors_combined <- sort(x = unlist(x = lapply(X = as.list(x = unique(x = unlist(x = starting_weights))), get_all_factors)))
      
      # Get largest common factor of all weights:
      largest_common_factor <- max(rle(all_factors_combined)$values[rle(all_factors_combined)$lengths == length(x = unique(x = unlist(x = starting_weights)))])
      
    }
    
    # If there are any constant characters:
    if (any(unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = min_max_matrix_list, apply, 1, diff), '==', 0), which), length)) > 0)) {
      
      # Get numbers of blocks with constant characters:
      blocks_with_constant_characters <- which(x = unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = min_max_matrix_list, apply, 1, diff), '==', 0), which), length)) > 0)
      
      # For each block with constant characters set weight to zero:
      for(i in blocks_with_constant_characters) starting_weights[[i]][which(x = apply(min_max_matrix_list[[i]], 1, diff) == 0)] <- 0
      
    }
    
    # If there are continuous characters:
    if (any(names(starting_weights) == "CONTINUOUS")) {
      
      # Get numbers of continuous blocks:
      continuous_blocks <- which(x = names(starting_weights) == "CONTINUOUS")
      
      # Set weights for continuous blocks as original weights:
      for(i in continuous_blocks) starting_weights[[i]] <- character_weights[[i]]
      
    }
    
    # Update weights:
    character_weights <- starting_weights
    
  }
  
  # Create top list:
  top_list <- list(header = text_lines, step_matrices = step_matrices)
  
  # Start to compile output starting with top list:
  cladistic_matrix <- list(topper = top_list)
  
  # For each matrix block:
  for(i in 1:length(x = matrix_block_list)) {
    
    # Create sublist for character information:
    characters <- list(symbols = symbols[[i]], missing = missing[[i]], gap = gap[[i]])
    
    # Build list for current block:
    current_block <- list(block_name = block_names[[i]], datatype = names(matrix_block_list)[i], matrix = matrix_block_list[[i]], ordering = ordering[[i]], character_weights = character_weights[[i]], minimum_values = min_max_matrix_list[[i]][, "min"], maximum_values = min_max_matrix_list[[i]][, "max"], characters = characters)
    
    # Store current block in output:
    cladistic_matrix[[(i + 1)]] <- current_block
    
    # Add name to current matrix (using a number to differentiate them):
    names(cladistic_matrix)[(i + 1)] <- paste("matrix_", i, sep = "")
    
  }
  
  # Assign class to cladistic_matrix:
  class(cladistic_matrix) <- "cladisticMatrix"

  # Return cladistic_matrix invisibly:
  cladistic_matrix
}
