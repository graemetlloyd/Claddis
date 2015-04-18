ReadMorphNexus <- function(file, equalise.weights=FALSE) {
    
# ADD WARNINGS FOR discrete.matrix AND continuous.matrix IF ROWS DO NOT HAVE THE RIGHT NUMBER OF CHARACTERS

  # Little piece of code to help deal with foreign characters:
  Sys.setlocale('LC_ALL','C')

  # Read in NEXUS file as raw text:
  X <- readLines(file, warn=FALSE)

  # Check that this is a #NEXUS file:
  if(length(grep("#NEXUS", X, ignore.case=T)) == 0) stop("This is not a #NEXUS file.")

  # Check that this is a #NEXUS file:
  if(length(grep("MATRIX", X, ignore.case=T)) == 0) stop("There is no matrix present.")

  # Are there blocks of the NEXUS file that can be removed? (i.e., trees, macclade, mesquite etc.:
  if(length(grep("begin trees|begin mesquite|begin mrbayes|begin macclade|begin treebase|begin notes|begin paup", X, ignore.case=T)) > 0) {

    # Get position of all Begin tags:
    all.begins <- grep("begin ", X, ignore.case=T)

    # Get position of all End tags:
    all.ends <- grep("end;|endblock;|end ;|endblock ;", X, ignore.case=T)

    # Check Begin and End are in balance:
    if(length(all.begins) != length(all.ends)) stop("Begin and End tags are not equal in number.")

    # Get rows for blocks to delete:
    block.rows <- match(grep("begin trees|begin mesquite|begin mrbayes|begin macclade|begin treebase|begin notes|begin paup", X, ignore.case=T), all.begins)

    # Find start lines:
    block.begins <- all.begins[block.rows]

    # Find end lines:
    block.ends <- all.ends[block.rows]

    # Incerementally remove superflouous blocks:
    for(i in length(block.begins):1) X <- X[c(c(1:(block.begins[i] - 1)), c((block.ends[i] + 1):length(X)))]

  }

  # Case if there are superflous taxon, character, or state labels:
  if(length(grep("taxlabels|charstatelabels|charlabels|statelabels", X, ignore.case=T)) > 0) {

    # Find label beginning rows:
    label.ends <- label.begins <- grep("taxlabels|charstatelabels|charlabels|statelabels", X, ignore.case=T)

    # For each label find the end tag that corresponds to it:
    for(i in 1:length(label.begins)) label.ends[i] <- min(grep(";", X, ignore.case=T)[grep(";", X, ignore.case=T) > label.begins[i]])

    # Incerementally remove superflouous labels:
    for(i in length(label.begins):1) X <- X[c(c(1:(label.begins[i] - 1)), c((label.ends[i] + 1):length(X)))]

  }

  # If there is text inside single quotes:
  if(length(grep("'", X)) > 0) {

    # Find items to replace (text not in single quotes):
    replacement.items <- unique(unlist(strsplit(gsub(paste("'[A-Z|\t| |0-9|a-z|\\?|\\-|\\(|)|\\.]*'", sep=""), "''", X[grep("'", X)]), "''")))

    # Make sure nothing "" is not in this list:
    replacement.items <- setdiff(unique(unlist(strsplit(replacement.items, " "))), "")

    # Need to add escape characters for regular expression to work:
    #replacement.items <- gsub("\\?", "\\\\?", replacement.items)

    # Make block that isolates lines with single quotes:
    lines.to.edit <- X[grep("'", X)]

    # Remove text outside single quotes:
    for(i in replacement.items[order(nchar(replacement.items), decreasing=TRUE)]) lines.to.edit <- gsub(i, "", lines.to.edit, fixed=T)

    # Remove double spaces:
    while(length(grep("  ", lines.to.edit)) > 0) lines.to.edit <- gsub("  ", " ", lines.to.edit, fixed=T)

    # Now isolate names within single quotes:
    names.in.single.quotes <- unique(gsub("'", "", strsplit(trim(paste(lines.to.edit, collapse="")), "' '")[[1]]))

    # Make sure nothing "" is not in this list:
    names.in.single.quotes <- which(!names.in.single.quotes == "")

    # Replace names in single quotes with underscored versions without single quotes or parentheses:
    for(i in names.in.single.quotes) X <- gsub(i, gsub("\\(|)", "", gsub(" ", "_", i)), X, fixed=T)

    # Remove all single quotes from text file:
    X <- gsub("'", "", X, fixed=T)

  }

  # Get rid of spaces around equals signs to make later regular expresssions work:
  while(length(grep(" = |= | =", X))) X <- gsub(" = |= | =", "=", X)

  # Remove weird characters (may need to add to the list or not if Sys.setlocale fixes the issues):
  X <- gsub("\x94|\x93|\xd5|\xd4|\xd3|\xd2|'", "", X)

  # Replace tabs wth spaces and trim leading and trailing spaces from each line:
  X <- apply(matrix(gsub("\t", " ", X)), 2, trim)

  # Delete any empty lines (if present):
  if(length(grep(TRUE, X == "")) > 0) X <- X[-grep(TRUE, X == "")]

  # If header text is present:
  if(length(grep("\\[", X)) > 0) {

    # Work out beginning and ending lines for text:
    textlines <- apply(cbind(setdiff(grep("\\[", X), grep("\\[[0-9:A-Z:a-z]{1}\\]", X)), setdiff(grep("\\]", X), grep("\\[[0-9:A-Z:a-z]{1}\\]", X))), 1, paste, collapse=":")

    # Convert beginning and endings to numerics for finding in vector:
    lines.to.delete <- textlines <- eval(parse(text=paste("c(", paste(textlines, collapse=","), ")", sep="")))  
  
    # Grab text and store:
    textlines <- paste(gsub("\\[|\\]", "", X[textlines]), collapse="\n")
  
    # Delete text lines to avoid regular expression errors later:
    X <- X[-lines.to.delete]

  # If header text is absent:
  } else {
    
    # Store empty text string:
    textlines <- ""

  }

  # Little test for valid NEXUS formatting for number of taxa:
  if(length(grep("ntax", X, ignore.case=T)) == 0) stop("Number of taxa not defined.")

  # Little test for valid NEXUS formatting for number of characters:
  if(length(grep("nchar", X, ignore.case=T)) == 0) stop("Number of characters not defined.")

  # Get number of taxa:
  ntax <- as.numeric(strsplit(strsplit(gsub(";|\t", "", X[grep("NTAX=|Ntax=|ntax=", X)]), "NTAX=|Ntax=|ntax=")[[1]][2], " ")[[1]][1])

  # Grab and format nchar line(s):
  nchar <- matrix(unlist(strsplit(gsub(";|\t", "", X[grep("NCHAR=|Nchar=|nchar=", X)]), "NCHAR=|Nchar=|nchar=")), ncol=2, byrow=TRUE)[,2]

  # Isolate just the number(s) of characters:
  for(i in 1:length(nchar)) nchar[i] <- strsplit(nchar[i], " ")[[1]][1]

  # Convert to numeric:
  nchar <- as.numeric(nchar)

  # Check matrix has size:
  if((nchar * ntax) == 0) stop("Matrix has no dimensions.")

  # If symbols are specified in the file:
  if(length(grep("symbols", X, ignore.case=T)) > 0) {

    # Get initial set of symbols:
    symbols <- strsplit(strsplit(strsplit(X[grep("symbols", X, ignore.case=T)], "SYMBOLS=|Symbols=|symbols=")[[1]][2], "\"")[[1]][2], " |")[[1]]
    
    # Collapse to just the symbols themselves:
    symbols <- symbols[nchar(symbols) == 1]

  # If symbols are not specified in the file:
  } else {
  
    # Use a default set of symbols:
    symbols <- c(c(0:9), LETTERS[1:22])

  }

# ADD BIT HERE FOR STANDARD AND DNA AND PROTEIN ETC.

  # If the missing character is specified:
  if(length(grep("missing", X, ignore.case=T)) > 0) {

    # Get missing character:
    missing <- strsplit(strsplit(X[grep("missing", X, ignore.case=T)][1], "MISSING=|Missing=|missing=")[[1]][2], "")[[1]][1]

  # If the missing character is not specified:
  } else {

    # Set as default symbol of the question mark:
    missing <- "?"

  }

  # If the gap character is specified:
  if(length(grep("gap", X, ignore.case=T)) > 0) {

    # Get gap character:
    gap <- strsplit(strsplit(X[grep("gap", X, ignore.case=T)][1], "GAP=|Gap=|gap=")[[1]][2], "")[[1]][1]

  # If the gap character is not specified:
  } else {

    # Set as default symbol of the dash:
    gap <- "-"

  }

  # Get start of matrix block(s):
  matrix.startlines <- setdiff(setdiff(grep("matrix", X, ignore.case=T), grep("stepmatrix", X, ignore.case=T)), grep(";", X))

  # Get end of matrix blocks:
  matrix.endlines <- setdiff(grep(";", X), grep("end", X, ignore.case=T))[setdiff(grep(";", X), grep("end", X, ignore.case=T)) > matrix.startlines]

  # For each block of matrix:
  for(i in length(matrix.startlines):1) {

    # If lines break (and not interleaved:
    if((matrix.endlines[i] - matrix.startlines[i] - 1) > ntax && length(grep("interleave", X, ignore.case=T)) == 0) {

      # Set current matrix block:
      current.matrixblock <- X[(matrix.startlines[i] + 1):(matrix.endlines[i] - 1)]

      # Remove all character symbols to find broken lines:
      naked.matrixblock <- gsub(paste(c(symbols, missing, gap, "\\{", "\\(", ")", "}", "\\?", "[0-9]", "\\."), collapse="|"), "", current.matrixblock)

      # Case if this seems to have identified the problem (i.e., once line breaks are removed number of lines is equal to number of taxa):
      if(sum(nchar(naked.matrixblock) > 0) == ntax) {

        # Get start rows for removing line breaks:
        start.rows <- which(nchar(naked.matrixblock) > 0)

        # Create empty block to store new matrixblock:
        new.matrixblock <- vector(mode="character")

        # For each taxon:
        for(j in 1:length(start.rows)) {

          # If taxon is not the last taxon get the block corresponding to it:
          if(j < length(start.rows)) taxon.block <- current.matrixblock[start.rows[j]:(start.rows[min(which(start.rows > start.rows[j]))] - 1)]

          # If taxon is the last taxon get the block corresponding to it:
          if(j == length(start.rows)) taxon.block <- current.matrixblock[start.rows[j]:length(current.matrixblock)]

          # Ensure there are spaces at the end of the taxon name:
          if(length(grep(" ", taxon.block[1])) == 0) taxon.block[1] <- paste(taxon.block[1], "  ", sep="")

          # Make into single line and add to new matrix block:
          new.matrixblock <- c(new.matrixblock, paste(taxon.block, collapse=""))

        }

        # Replace matrix block with line breaks with one without:
        X <- c(X[1:matrix.startlines[i]], new.matrixblock, X[matrix.endlines[i]:length(X)])

      # Case if apparent line breaks do not add up to nmber of taxa:
      } else {

        # Stop and print warning.:
        stop("Matrix block(s) have line breaks. Ensure each taxon is on a single line.")

      }

    }

  }

  # if there are interelaved lines:
  if(length(grep("interleave", X, ignore.case=T)) > 0) {

    # Get start of matrix block(s):
    matrix.startlines <- setdiff(setdiff(grep("matrix", X, ignore.case=T), grep("stepmatrix", X, ignore.case=T)), grep(";", X))

    # Get end of matrix block(s):
    matrix.endlines <- setdiff(grep(";", X), grep("end", X, ignore.case=T))[setdiff(grep(";", X), grep("end", X, ignore.case=T)) > matrix.startlines]

    # For each block:
    for(i in length(matrix.startlines):1) {

      # Case if more lines than expected, number of lines is divisible by number of taxa and is a repetition of the number of taxa (i.e., is interleaved):
      if((matrix.endlines[i] - matrix.startlines[i] - 1) > ntax && (matrix.endlines[i] - matrix.startlines[i] - 1) %% ntax == 0 && (matrix.endlines[i] - matrix.startlines[i] - 1) / ntax > 1) {

        # Get number of repeated taxa (i.e. how many lines is matrix broken over):
        number.of.repeats <- (matrix.endlines[i] - matrix.startlines[i] - 1) / ntax

        # Get initial start lines:
        startlines <- 1:number.of.repeats * ntax - ntax + 1

        # Get intiial end lines:
        endlines <- 1:number.of.repeats * ntax

        # Get actual start line:
        startline <- (matrix.startlines[i] + 1)

        # Get actual end line:
        endline <- (matrix.endlines[i] - 1)

        # Update startlines:
        startlines <- startlines + startline - 1

        # Update endlines:
        endlines <- endlines + startline - 1

        # For each start line:
        for(j in 1:length(startlines)) {

          # If first start line create new block:
          if(j == 1) new.block <- X[startlines[j]:endlines[j]]

          # For subsequent start lines add to new block with matrix begin and end tags:
          if(j > 1) new.block <- c(new.block, ";","END;", "MATRIX", X[startlines[j]:endlines[j]])

        }

        # Update file with ne blocks:
        X <- c(X[1:(startline - 1)], new.block, X[(endline + 1):length(X)])

      }

    }

  }

  # First sort of getting rows that correspond to the actual data matri(ces):
  matrixblocks <- matrix(rep(1:ntax, length(setdiff(setdiff(grep("matrix", X, ignore.case=T), grep("stepmatrix", X, ignore.case=T)), grep(";", X)))), ncol=length(setdiff(setdiff(grep("matrix", X, ignore.case=T), grep("stepmatrix", X, ignore.case=T)), grep(";", X))))
  
  # Add start values:
  for(i in 1:length(setdiff(setdiff(grep("matrix", X, ignore.case=T), grep("stepmatrix", X, ignore.case=T)), grep(";", X)))) matrixblocks[, i] <- matrixblocks[, i] + (setdiff(setdiff(grep("matrix", X, ignore.case=T), grep("stepmatrix", X, ignore.case=T)), grep(";", X))[i])

  # Convert to vector for use:
  matrixblocks <- sort(matrixblocks)

  # Get matrix block:
  matrixblock <- X[matrixblocks]

  # Little check for length violation:
  if(length(sort(match(";", matrixblock))) > 0) stop("Number of lines of matrix is less than the number of taxa. Check line breaks and NTAX.")

  # Collapse all double spaces to single spaces in order to extract row names:
  while(length(grep("  ", matrixblock)) > 0) matrixblock <- gsub("  ", " ", matrixblock)

  # Little function for use in lapply on the line below:
  getfirstitem <- function(x) return(x[1])

  # Get row (taxon) names:
  row.names <- unlist(lapply(strsplit(matrixblock, " "), getfirstitem))

  # Strip out row (taxon) names:
  matrixblock <- gsub(gsub("\\?", "\\\\?", paste(row.names, " ", sep="", collapse="|")), "", matrixblock)

  # Relaces curly braces with parentheses for polymorphisms:
  matrixblock <- gsub("}", ")", (gsub("\\{", "(", matrixblock)))

  # If there are (potential) polymorphisms in the matrix block:
  if(length(grep("\\([0-9:A-Z:a-z]{1}[ |,|/|\\&]{1}", matrixblock)) > 0) {
  
    # First collapse lines with polymorphisms to just polymorphisms:
    polymorphism.combos <- gsub(paste(")[", paste(c(symbols, missing, gap, "\\?"), collapse="|"), "]*", sep=""), ")", gsub(paste("[", paste(c(symbols, missing, gap, "\\?"), collapse="|"), "]*\\(", sep=""), "\\(", matrixblock))
  
    # Now get just the polymorphisms
    polymorphism.combos <- unique(strsplit(paste(polymorphism.combos[grep(")", polymorphism.combos)], collapse=""), ")|\\(")[[1]])[-1]

    # Get the symbols that break up the polymorphism (e.g., space, comma, slash, ampersand):
    split.symbols <- unique(strsplit(paste(gsub(paste(symbols, collapse="|"), "", polymorphism.combos), collapse=""), "")[[1]])

    # Get replacements for polymorphism combinations:
    polymorphism.combo.replacements <- gsub(paste(split.symbols, collapse="|"), "", polymorphism.combos)

    # Ensure all polymorphisms combinations containing a space, ampersand, slash or comma are reduced to just the regular polymorphism:
    for(i in 1:length(polymorphism.combos)) matrixblock <- gsub(paste("\\(", polymorphism.combos[i], ")", sep=""), paste("\\(", polymorphism.combo.replacements[i], ")", sep=""), matrixblock)

  }

# THE BIT BELOW NEEDS TO BE MADE MORE LIKE THE BIT ABOVE (ALTHOUGH STRAY NAMES MAKES THIS TRICKY)

  # If there are spaces in every row of the matrix block (SHOULD ALSO ONLY BE DISCRETE MATRIX BLOCK(S)):
  if(length(grep(" ", matrixblock)) == length(matrixblock)) {

    # Little function to generate all possible character combinations (to help find those with spaces between that can be removed):
    CharCombs <- function(chars) {

      # Store combinations (both directions, i.e. 0 1 and 1 0):
      combos <- apply(rbind(t(combn(chars, 2)), cbind(t(combn(chars, 2))[, 2], t(combn(chars, 2))[, 1])), 1, paste, collapse=" ")

      # Add repeats (e.g., 0 0):
      combos <- c(apply(rbind(chars, chars), 2, paste, collapse=" "), combos)

      # Output combinations:
      return(combos)

    }

    # Find all possible character combinations for matrix:
    character.combinations <- CharCombs(gsub("\\?", "\\\\?", c(symbols, gap, missing, "\\(", ")")))

    # For each character combination:
    for(i in character.combinations) {

      # Whilst it exists, replace it without space gap:
      while(length(grep(i, matrixblock)) > 0) matrixblock <- gsub(i, gsub(" ", "", i), matrixblock)

    }

  }

  # Check for any remaining stray row name parts (i.e., detached by spaces from the rest of the name):
  while(length(grep(" ", matrixblock)) > 0) {

    # Identify problem rows:
    problem.rows <- grep(" ", matrixblock)

    # Identify stray name parts:
    stray.name.parts <- unlist(lapply(strsplit(matrixblock[problem.rows], " "), getfirstitem))

    # Reunite stray name parts with rest of name:
    row.names[problem.rows] <- paste(row.names[problem.rows], stray.name.parts, sep="_")

    # Remove stray name parts from matrixblock:
    matrixblock <- gsub(gsub("\\(", "\\\\(", paste(paste(stray.name.parts, " ", sep=""), collapse="|")), "", matrixblock)

  }

  # Create empty matrix to store character-taxon matrix:
  MATRIX <- matrix(nrow=ntax, ncol=sum(nchar))

  # If there are multiple matrix block elements:
  if(length(row.names) > ntax) {

    # Establish row order of first block:
    row.order <- row.names[1:ntax]

    # Ensure subsequent blocks share same order:
    for(i in 1:(length(row.names) / ntax) - 1) matrixblock[((ntax * i) + 1):(ntax * (i + 1))] <- matrixblock[((ntax * i) + 1):(ntax * (i + 1))][match(row.order, row.names[((ntax * i) + 1):(ntax * (i + 1))])]

  }

  # Vector to store taxon names:
  MATRIXrn <- row.names[1:ntax]
    
  # Add taxon names to character-taxon matrix:
  rownames(MATRIX) <- MATRIXrn

  # Convert matrix block to actual matrix with nrows equal to number of taxa:
  matrixblock <- matrix(matrixblock, nrow=ntax)

  # Get discrete columns:
  discrete.columns <- setdiff(1:ncol(matrixblock), grep(TRUE, unlist(lapply(apply(matrixblock, 2, grep, pattern="\\."), length)) > 0))

  # Collapse discrete columns and repeat:
  matrixblock[, discrete.columns] <- apply(as.matrix(matrixblock[, discrete.columns], nrow=ntax), 1, paste, collapse="")

  # If repeats:
  if(length(discrete.columns) > 1) {

    # Get rid of repeats:
    matrixblock <- matrixblock[, -discrete.columns[2:length(discrete.columns)]]

    # Ensure is still a matrix:
    if(!is.matrix(matrixblock)) matrixblock <- as.matrix(matrixblock, nrow=ntax)

  }

  # Find columns in matrix that contain continuous data:
  continuous.columns <- grep(TRUE, unlist(lapply(apply(matrixblock, 2, grep, pattern="\\."), length)) > 0)

  # Find columns in matrix that contain discrete data:
  discrete.columns <- setdiff(c(1:ncol(matrixblock)), continuous.columns)

  # Get character numbers that correspond to each column (columns are from and to):
  column.numbers <- cbind(c(1, (nchar + 1)[-length(nchar)]), cumsum(nchar))

  # If there is continuous data:
  if(length(continuous.columns) > 0) {

    # Make continuous characters into a matrix block (with columns for each character):
    continuous.matrix <- matrix(unlist(strsplit(apply(matrix(matrixblock[, continuous.columns], ncol=length(continuous.columns)), 1, paste, collapse=" "), " ")), ncol=sum(nchar[continuous.columns]), byrow=TRUE)

    # Replace missing symbol with NA:
    continuous.matrix[continuous.matrix == missing] <- NA

    # Replace gap symbol with NA:
    continuous.matrix[continuous.matrix == gap] <- NA

    # Make from-to matrix for continuous character numbers:
    continuous.character.numbers <- matrix(column.numbers[continuous.columns, ], ncol=2)

    # Modify matrix to store all character numbers in to column:
    for(i in 1:nrow(continuous.character.numbers)) continuous.character.numbers[i, 2] <- paste(continuous.character.numbers[i, 1]:continuous.character.numbers[i, 2], collapse=" ")

    # Grab character numbers and convert back into numerics:
    continuous.character.numbers <- as.numeric(strsplit(paste(continuous.character.numbers[, 2], collapse=" "), " ")[[1]])

    # Store continuous characters in full matrix:
    MATRIX[, continuous.character.numbers] <- continuous.matrix

  }

  # If there is discrete data:
  if(length(discrete.columns) > 0) {

    # Grab columns of discrete data and form into matrix:
    discrete.matrix <- matrix(matrixblock[, discrete.columns], ncol=length(discrete.columns))

    # Delete any remaining whitespace:
    discrete.matrix <- gsub("\n|\t| ", "", discrete.matrix)

    # Collapse matrix to a single column:
    discrete.matrix <- matrix(apply(discrete.matrix, 1, paste, collapse=""), ncol=1)

    # Convert into list by splitting into single characters:
    discrete.matrix <- strsplit(discrete.matrix, "")

    # Get numeric characters to replace symbols with:
    symbol.replacements <- 0:(length(symbols) - 1)

    # For each row of the discrete matrix:
    for(i in 1:length(discrete.matrix)) {
  
      # For each symbol:
      for(j in 1:length(symbols)) {
  
        # Replace symbol with appropriate numerical equivalent:
        discrete.matrix[[i]][discrete.matrix[[i]] == symbols[j]] <- symbol.replacements[j]
  
      }

    }
  
    # Find rows with parentheses and by extension polymorphisms:
    polymorphism.rows <- grep(TRUE, unlist(lapply(lapply(discrete.matrix, grep, pattern="\\("), length)) > 0)

    # As long as there are polymorphisms:
    if(length(polymorphism.rows) > 0) {
  
      # For eaach row containing polymorphisms:
      for(i in polymorphism.rows) {

        # Store ith row as active row:
        active.row <- discrete.matrix[[i]]

        # Continue modifying ith row until all polymorphisms are specified correctly:
        while(length(grep("\\(", active.row)) > 0) {

          # Define columns that correspond to the first polymorphisms encountered:
          column.hits <- match("(", active.row):match(")", active.row)

          # Define polymorphisms and store:
          active.row[column.hits[1]] <- paste(active.row[column.hits[2:(length(column.hits) - 1)]], collapse="&")

          # Delete now redundnat additional characters:
          active.row <- active.row[-column.hits[2:length(column.hits)]]

        }

        # Store row back in discrete matrix:
        discrete.matrix[[i]] <- active.row

      }

    }
  
    # Little test to check the number of characters on each line is the same as nchar:
    if(length(grep(TRUE, unlist(lapply(discrete.matrix, length)) != nchar)) > 0) stop("Some lines have too many or too few characters.")

    # Convert into matrix format:
    discrete.matrix <- matrix(unlist(discrete.matrix), byrow=TRUE, nrow=ntax)

    # Replace missing symbol with NA:
    discrete.matrix[discrete.matrix == missing] <- NA

    # Replace gap symbol with NA:
    discrete.matrix[discrete.matrix == gap] <- NA

    # Make from-to matrix for discrete character numbers:
    discrete.character.numbers <- matrix(column.numbers[discrete.columns, ], ncol=2)

    # Modify matrix to store all character numbers in to column:
    for(i in 1:nrow(discrete.character.numbers)) discrete.character.numbers[i, 2] <- paste(discrete.character.numbers[i, 1]:discrete.character.numbers[i, 2], collapse=" ")

    # Grab character numbers and convert back into numerics:
    discrete.character.numbers <- as.numeric(strsplit(paste(discrete.character.numbers[, 2], collapse=" "), " ")[[1]])

    # Store discrete characters in full matrix:
    MATRIX[, discrete.character.numbers] <- discrete.matrix

  }

  # Get rid of spaces around dashes to make later regular expresssions work:
  while(length(grep(" - |- | -", X))) X <- gsub(" - |- | -", "-", X)

  # If minimum level of ordering is specified:
  if(length(grep("DEFTYPE", X)) > 0) {

    # Set default ordering of characters:
    default.ordering <- strsplit(strsplit(X[grep("deftype", X, ignore.case=T)], "DEFTYPE=|Deftype=|deftype=")[[1]][2], " ")[[1]][1]

  # If minimum level of ordering is not specified:
  } else {

    # Set default ordering to unordered:
    default.ordering <- "unord"

  }

# BELOW NEEDS SOME KIND OF IGNORE.CASE MODIFIER

  # If default ordering is ordered:
  if(default.ordering == "ord" || default.ordering == "Ord" || default.ordering == "ORD") {

    # Set all characters as ordered:
    ordering <- rep("ord", sum(nchar))

  # If default ordering is unordered:
  } else {

    # Set all characters as unordered:
    ordering <- rep("unord", sum(nchar))

  }

  # If there are continuous characters specify these:
  if(length(continuous.columns) > 0) ordering[continuous.character.numbers] <- "CONT"

  # Create null variable for step matrices:
  step.matrices <- NULL

  # Special case if there are user-defined characters, i.e. step matrices:
  if(length(grep("USERTYPE", X, ignore.case=T)) > 0) {

    # Get rows corresponding to start of stepmatrices:
    step.matrix.rows <- grep("USERTYPE", X, ignore.case=T)
  
    # Create empty list to store step matrices:
    step.matrices <- list()

    # Little check in case of abnormally large number of step matrices:
    if(length(step.matrix.rows) > length(LETTERS)) stop("Not enough letters to store the number of different step matrices.")
  
    # For each step matrix:
    for(i in step.matrix.rows) {
  
      # Grab text block corresponding to step matrix:
      step.matrix.block <- X[i:(i + as.numeric(strsplit(X[i], "\\(STEPMATRIX\\)=")[[1]][2]) + 1)]
  
      # Get the step matrix as a matrix (might still need to think about how to define e.g. infinity):
      step.matrix <- gsub("\\.", "0", matrix(unlist(strsplit(step.matrix.block[3:length(step.matrix.block)], " ")), ncol=as.numeric(strsplit(X[i], "\\(STEPMATRIX\\)=")[[1]][2]) + 1, byrow=TRUE)[, -1])

      # Add row and column names:
      rownames(step.matrix) <- colnames(step.matrix) <- symbols[1:nrow(step.matrix)]

      # Add step matrix to list:
      step.matrices[[length(step.matrices) + 1]] <- step.matrix

      # Find step matrix name:
      step.matrix.name <- strsplit(strsplit(X[i], "USERTYPE ")[[1]][2], " ")[[1]][1]

      # Replace stpe matrix name with step_N:
      X <- gsub(step.matrix.name, paste("step_", LETTERS[length(step.matrices)], sep=""), X)

      # Use step matrix name (step_N) in list for later calling:
      names(step.matrices)[length(step.matrices)] <- paste("step_", LETTERS[length(step.matrices)], sep="")

    }

  }
  
  # Deal with ordering, if present:
  if(length(grep("TYPESET", X, ignore.case=T)) > 0) {
        
    # Find appropriate information line:
    ordering.line <- X[grep("TYPESET", X, ignore.case=T)]

    # Make sure ordering is on single line or rest of code will not work:
    if(length(grep(";", ordering.line)) == 0) stop("Ordering breaks over multiple lines. Place on single line.")

    # Get clean text:
    ordering.line <- trim(gsub(",|;", "", strsplit(ordering.line, "=")[[1]][2]))

    # Number of character types specified:
    ntype <- nchar(ordering.line) - nchar(gsub(":", "", ordering.line))
  
    # Get list of character types:
    types <- gsub("\\", "", strsplit(gsub("[0-9]|-| ", "", ordering.line), ":")[[1]][1:ntype], fixed=TRUE)
  
    # Get initial list of numbers of characters of each type:
    ordering.numbers <- strsplit(ordering.line, paste(paste(types, ": ", sep=""), collapse="|"))[[1]]

    # Conditional to deal with weird slashes in numbering:
    if(length(grep("\\3", ordering.numbers, fixed=T)) > 0) {

      # Identify effected rows:
      effected.rows <- grep("\\3", ordering.numbers, fixed=T)

      # For each effected row:
      for(i in 1:length(effected.rows)) {

        # Isolate data:
        effected.row <- ordering.numbers[effected.rows[i]]

        # Split by spaces:
        effected.row <- strsplit(effected.row, " ")[[1]]

        # For each isolated element fix the issue:
        for(j in grep("\\", effected.row, fixed=T)) effected.row[j] <- gsub("-", " ", gsub("\\3", "", effected.row[j], fixed=T))

        # Update ordering numbers:
        ordering.numbers[effected.rows[i]] <- paste(effected.row, collapse=" ")

      }

    }

    # Modify further ready for parsing and evaluating:
    ordering.numbers <- paste("c(", gsub(" ", ",", trim(gsub("-", ":", ordering.numbers[nchar(ordering.numbers) > 0]))), ")", sep="")
  
    # For each character type:
    for(i in 1:length(types)) {
  
      # Update ordering vector with character type:
      ordering[eval(parse(text=ordering.numbers[i]))] <- types[i]

    }

  }

  # Set default weight set of all equal (and 1):
  weights <- rep(1, sum(nchar))
    
  # Deal with differential weighting, if present:
  if(length(grep("WTSET", X, ignore.case=T)) > 0) {

    # Grab weighting line:
    weighting.line <- gsub(",|;", "", strsplit(X[grep("WTSET", X, ignore.case=T)], "=")[[1]][2])
  
    # Isolate weight values:
    weight.values <- gsub(":", "", strsplit(weighting.line, " ")[[1]][grep(":", strsplit(weighting.line, " ")[[1]])])
  
    # Get initial list of numbers of characters of each type:
    weighting.numbers <- strsplit(weighting.line, paste(paste(weight.values, ": ", sep=""), collapse="|"))[[1]]
  
    # Modify further ready for parsinga nd evaluating:
    weighting.numbers <- paste("c(", gsub(" ", ",", trim(gsub("-", ":", weighting.numbers[nchar(weighting.numbers) > 0]))), ")", sep="")
  
    # For each weight value:
    for(i in 1:length(weight.values)) {
  
      # Store updated weights:
      weights[eval(parse(text=weighting.numbers[i]))] <- as.numeric(weight.values[i])
  
    }
  
  }

  # Vectors to store minimum and maximum values:
  min.vals <- max.vals <- vector(length=nchar)
    
  # For each character find the minimum and maximum value:
  for(i in 1:nchar) {
        
    # List non-missing values for character:
    non.nas <- sort(MATRIX[, i])
    
    # As long as there are non-NAs (i.e., there is not missing data for all taxa):
    if(length(non.nas) > 0) {
  
      # Case if polymorphisms are available:
      if(length(grep("&", non.nas)) > 0) {
            
        # Get polymorphisms:
        polymorphisms <- grep("&", non.nas)
            
        # For each polymorphism:
        for(j in polymorphisms) {
                
          # Add split values to end of non.nas vector:
          non.nas <- c(non.nas, strsplit(non.nas[j], "&")[[1]])

        }
            
        # Remove non-split values from non.nas vector:
        non.nas <- non.nas[-polymorphisms]

      }
        
      # Find maxima:
      max.vals[i] <- max(as.numeric(non.nas))

      # Find minima:
      min.vals[i] <- min(as.numeric(non.nas))

    # Case if no data for any taxon:
    } else {

      # Set both minimum and maximum values as zero:
      max.vals[i] <- min.vals[i] <- 0

    }

  }

  # If wanting to equalise weights (NB: all continuous characters will continue to be weighted one):
  if(equalise.weights) {

    # Get starting weights:
    weights <- apply(rbind(c(max.vals - min.vals), rep(1, nchar)), 2, max)
  
    # Updaye weights for unordered characters:
    weights[ordering == "unord"] <- 1
  
    # Update weights for ordered characters:
    weights[ordering == "ord"] <- 1 / weights[ordering == "ord"]
  
    # If there are step matrices:
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

  # Remove characters not favoured by phylogenetic software (periods, spaces, slashes, dashes and plusses):
  rownames(MATRIX) <- gsub("\\.", "_dot_", gsub("-", "_dash_", gsub("/", "_slash_", gsub("\\+", "_plus_", gsub(" ", "_", rownames(MATRIX))))))
    
  # Create output formatted data:
  result <- list(textlines, MATRIX, ordering, weights, max.vals, min.vals, step.matrices, symbols)

  # Add names to result:
  names(result) <- c("header", "matrix", "ordering", "weights", "max.vals", "min.vals", "step.matrices", "symbols")

  # Return result:
  return(result)

}
