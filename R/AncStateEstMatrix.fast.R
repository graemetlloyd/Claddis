#' Ancestral State Estimation
#' 
#' Given a tree and a cladistic matrix uses likelihood to estimate the ancestral states for every character.
#' 
#' Uses the \link{rerootingMethod} (Yang et al. 1995) as implemented in the \link{phytools} package to make ancestral state estimates. Here these are collapsed to the most likely state, or if two or more states are most likely, a polymorphism of the most likely states. This is the method used by Brusatte et al. (2014).
#' 
#' @param morph.matrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param tree A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{morph.matrix}.
#' @param estimate.allchars An optional that allows the user to make estimates for all ancestral values. The default will only make estimates for nodes that link coded terminals.
#' @param estimate.tips An optional that allows the user to make estimates for tip values. The default only makes estimates for internal nodes.
#'
#' @return \item{anc.lik.matrix}{A matrix of nodes (hypothetical ancestors; rows) against characters (columns) listing the reconstructed ancestral states.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. Current Biology, 24, 2386-2392.
#' 
#' Yang, Z., Kumar, S. and Nei, M., 1995. A new method of inference of ancestral nucleotide and amino acid sequences. Genetics, 141, 1641-1650.
#'
#' @examples
#' 
#' # Set random seed:
#' set.seed(17)
#' 
#' # Generate a random tree for the Michaux data set:
#' tree <- rtree(nrow(Michaux1989$matrix))
#' 
#' # Update taxon names to match those in the data matrix:
#' tree$tip.label <- rownames(Michaux1989$matrix)
#' 
#' # Set root time by making youngest taxon extant:
#' tree$root.time <- max(diag(vcv(tree)))
#' 
#' # Estimate ancestral states:
#' AncStateEstMatrix(Michaux1989, tree)
#' 
#' @export AncStateEstMatrix
AncStateEstMatrix.fast <- function(morph.matrix, tree, estimate.allchars = FALSE, estimate.tips = FALSE, uncertainty.threshold = 0) {

  # Catch problem with trees with no branch lengths:
  if(is.null(tree$edge.length)) stop("ERROR:\n Tree must have branch lengths.")

  # Catch problem with polytomies:
  if(tree$Nnode < (Ntip(tree) - 1)) stop("ERROR:\n Tree must be fully bifurcating.")

  # Catch problem with zero-length branches:
  if(any(tree$edge.length == 0)) stop("ERROR:\n Tree must not have zero-length branches.")

  #TG: This bit is actually solved when creating the character trees (chartree) that are forced to have no labels
  #TG: One advantage of leaving the node.labels on the tree is to allow them to be returned in the output matrix
  # # Remove node labels from tree (causes bug downstream):
  # tree$node.label <- NULL

  # Collapse matrix to vectors for each character (state and ordering combination):
  collapse.matrix <- apply(rbind(morph.matrix$matrix, morph.matrix$ordering), 2, paste, collapse="")

  # Find just unique characters (no point repeating ancestral state reconstruction if codings and ordering are identical):
  unique.characters <- match(unique(collapse.matrix), collapse.matrix)

  # Case if only estimating states for ancestral nodes:
  if(estimate.tips == FALSE) {
    
    # Create ancestral storage matrix:
    anc.lik.matrix <- matrix(nrow=Nnode(tree), ncol=length(morph.matrix$matrix[1, ]))
    
    # Label matrix to record ancestral state estimates:
    rownames(anc.lik.matrix) <- c((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
    
  # Case if also estimating states for tip nodes:
  } else {

    # Create ancestral storage matrix (including tips):
    anc.lik.matrix <- matrix(nrow=Ntip(tree) + Nnode(tree), ncol=length(morph.matrix$matrix[1, ]))

    # Label matrix to record ancestral state estimates:
    rownames(anc.lik.matrix) <- c(tree$tip.label, (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
    
  }

  #~~~~~~~
  #TG: Isolating the character matrix as a list per characters
  morph.matrix.character <- as.list(data.frame(morph.matrix$matrix, stringsAsFactors = FALSE))
  #TG: Isolating the ordering as a list per characters
  morph.matrix.ordering <- as.list(morph.matrix$ordering)
  #TG: Adding characters IDs
  names(morph.matrix.character) <- seq(from = 1, to = ncol(morph.matrix$matrix))
  names(morph.matrix.ordering) <- seq(from = 1, to = ncol(morph.matrix$matrix))

  #~~~~~
  #TG: Getting max and min values as a list (each element corresponds to a character and contains respectively the min and the max val)
  minmaxvals <- as.list(data.frame(t(cbind(morph.matrix$min.vals, morph.matrix$max.vals)), stringsAsFactors = FALSE))
  names(minmaxvals) <- seq(from = 1, to = ncol(morph.matrix$matrix))

  #~~~~~
  #TG: Dealing with the simple characters first (i.e. all NA or constant)
  #TG: these two following parts excludes characters for which ancestral states should not be calculated.
  #TG: It excludes the NA and/or constant characters from the morph.matrix.character/morph.matrix.ordering/minmaxvals lists allowing these lists to be treated without constant characters bits.

  #~~~~~
  #TG: Function for setting constant states
  set.constant.state <- function(state, tree) {
    constant_character <- rep(state, (Ntip(tree) + Nnode(tree)))
    if(is.null(tree$node.label)) {
      names(constant_character) <- c(tree$tip.label, (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
    } else {
      names(constant_character) <- c(tree$tip.label, tree$node.label)
    }
    return(constant_character)
  }

  #TG: Dealing with characters that have only NAs
  NA.characters <- unlist(lapply(morph.matrix.character, function(X) all(is.na(X))))
  #TG: If any characters have NAs, remove them from the lists and store their ID in the NA.characters object
  if(any(NA.characters)) {
    #TG: Creating the list of constant characters
    NA.temp <- rep(NA, length(which(NA.characters)))
    names(NA.temp) <- names(which(NA.characters))
    NA.characters <- lapply(as.list(NA.temp), set.constant.state, tree)
    #TG: remove the NA characters from the initial lists
    morph.matrix.character[names(NA.characters)] <- NULL
    morph.matrix.ordering[names(NA.characters)] <- NULL
    minmaxvals[names(NA.characters)] <- NULL
  } else {
    #TG: No NA characters
    NA.characters <- NULL
  }

  #TG: Dealing with the constant characters
  #~~~~~~
  #TG: Function for getting the constant characters (with their states)
  get.constant.characters <- function(character) {
    if(length(unique(character)) == 1) return(unique(character))
   }
  #TG: get the constant characters
  constant.characters <- unlist(lapply(morph.matrix.character, get.constant.characters))
  #TG: If any constant characters, remove them from the lists and store their ID in the constant.character object along with their state
  if(length(constant.characters) > 0) {
    #TG: Create the list of constant characters
    constant.characters <- lapply(as.list(constant.characters), set.constant.state, tree)
    #TG: remove the constant characters from the initial lists
    morph.matrix.character[names(constant.characters)] <- NULL
    morph.matrix.ordering[names(constant.characters)] <- NULL
    minmaxvals[names(constant.characters)] <- NULL    
  } else {
    #TG: no constant characters
    constant.characters <- NULL
  }


  #########
  # Former applicable characters loop
  #########
  #TG: applied version of the loop going through each applicable characters

      
  #~~~~~
  #TG: Function for treating the missing values as all possible states (if estimate.allchars == TRUE)
  set.NA.all.possible.states <- function(morph.matrix.character, minmaxvals) {
    morph.matrix.character[is.na(morph.matrix.character)] <- paste(minmaxvals[1]:minmaxvals[2], collapse="&")
    return(morph.matrix.character)
  }
  
  # If estimating states for all taxa then treat missing values as all possible states:
  if(estimate.allchars) {
    morph.matrix.character <- mapply(set.NA.all.possible.states, morph.matrix.character, minmaxvals, SIMPLIFY = FALSE)
  }

  #~~~~~
  #TG: function for finding the tips which cannot be used due to missing data:
  find.tips.to.go <- function(morph.matrix.character, morph.matrix) {
    return(rownames(morph.matrix$matrix)[which(is.na(morph.matrix.character))])
  }

  # Find tips which cannot be used due to missing data:
  tipstogo <- lapply(morph.matrix.character, find.tips.to.go, morph.matrix)
      
  # Only continue if at least three tips in pruned tree:
  #TG: detect if any tips to go is > (Ntip(tree) -2) (leading to trees with less than three tips)
  if(any(unlist(lapply(tipstogo, length)) > (Ntip(tree) -2))) {
    stop("DEBUG: not developed yet!")
    #TG: probably have to remove the characters...
  }

  #~~~~
  #TG: function for getting the character tree
  get.the.character.trees <- function(tipstogo, tree) {
    # If there are tips to remove create a pruned tree:
    if(length(tipstogo) > 0) chartree <- drop.tip(tree, tipstogo)
  
    # If there are no tips to remove just use full tree:
    if(length(tipstogo) == 0) chartree <- tree

    return(chartree)
  }

  #TG: Getting the trees for each character
  chartrees <- lapply(tipstogo, get.the.character.trees, tree)
  
  #~~~~~
  #TG: Function for getting the tip values from the character trees
  get.tip.values <- function(morph.matrix.character, chartrees, morph.matrix) {
    # Getting the tip values
    tipval <- morph.matrix.character[match(chartrees$tip.label, rownames(morph.matrix$matrix))]
    return(tipval)
  }

  # Get tip values for the pruned tree:
  tipvals <- mapply(get.tip.values, morph.matrix.character, chartrees, MoreArgs = list(morph.matrix), SIMPLIFY = FALSE)
  # Get tip names for the pruned tree:
  tipnames <- lapply(chartrees, function(X) return(X$tip.label))


  #########
  # Former non gap-coded condition (in the loop)
  #########
  #TG: applied version of the part of the applicable characters loop dealing with characters that are not gap-coded or unordered


  # # Continue if not gap-coded (or apparently gap-coded) or is unordered:
  # if(!any(diff(sort(as.numeric(unlist(strsplit(tipvals, split = "&"))))) > 1) || morph.matrix.ordering[[i]] == "unord") { #REPLACE OUT OF LOOP
  # #TG: I'm not sure about this part? The aim here is to select only characters that are not ordered and are not multi-state (i.e. no "&")?


  #~~~~~
  #TG: Function for getting the models for each characters
  set.my.models <- function(minmaxvals, morph.matrix.ordering) {
    # Set discrete character estimation model if unordered and or binary character:
    if(minmaxvals[2] - minmaxvals[1] == 1 || minmaxvals[2] - minmaxvals[1] > 1 && morph.matrix.ordering == "unord") mymodel <- "ER"

    # Set discrete character estimation model if ordered multistate character:
    if(minmaxvals[2] - minmaxvals[1] > 1 && morph.matrix.ordering == "ord") {
              
    # Create all zero matrix:
    mymodel <- matrix(0, nrow=(minmaxvals[2] - minmaxvals[1]) + 1, ncol=(minmaxvals[2] - minmaxvals[1]) + 1)
              
    # Name rows and columns as states:
    rownames(mymodel) <- colnames(mymodel) <- minmaxvals[1]:minmaxvals[2]
              
    # Enter one for all the off-diagonal diagonals (an ordered change model):
    for(j in 1:(length(mymodel[1, ]) - 1)) mymodel[j + 1, j] <- mymodel[j, j + 1] <- 1
    }

    return(mymodel)
  }
  
  # Getting all the models for each character            
  mymodels <- mapply(set.my.models, minmaxvals, morph.matrix.ordering, SIMPLIFY = FALSE)

  #~~~~~
  #TG: Creating a matrix to store probabilities of tip values:
  create.tipvals.mat <- function(character.list, tipvals, minmaxvals, tipnames) {
    # Create matrix to store probabilities of tip values:
    tipvals.mat <- matrix(0, nrow=length(tipvals[[character.list]]), ncol=minmaxvals[[character.list]][2] - minmaxvals[[character.list]][1] + 1)

    # Add rownames (tip labels):
    rownames(tipvals.mat) <- tipnames[[character.list]]

    # Add colunames (state values):
    colnames(tipvals.mat) <- minmaxvals[[character.list]][1]:minmaxvals[[character.list]][2]
    return(tipvals.mat)
  }

  #TG: set up the list of characters to deal with
  character.list <- as.list(seq(from = 1, to = length(tipvals)))
  # Create matrices to store probabilities of tip values:
  tipvals.mat <- lapply(character.list, create.tipvals.mat, tipvals, minmaxvals, tipnames)
  
  #~~~~        
  #TG: Function for filling the observed probabilities (dealing with polymorphism)
  fill.observed.probabilities <- function (tipvals.mat, tipvals) {
    #Fill all probabilities equal to one (non-polymorphisms):
    for(j in colnames(tipvals.mat)) tipvals.mat[which(tipvals == j), j] <- 1

    # If there are polymorphisms make all observed states equally probable:
    if(any(apply(tipvals.mat, 1, sum) == 0)) {
                
      # Get list of tip values with polymorphisms:
      polymorphism.values <- which(apply(tipvals.mat, 1, sum) == 0)
                
      # Go through each polymorphism:
      for(j in polymorphism.values) {
                    
        # Get list of each state:
        states <- strsplit(tipvals[j], "&")[[1]]
                    
        # Make each state equally probable in tip values matrix:
        tipvals.mat[j, states] <- 1 / length(states)
      }
    }
    return(tipvals.mat)
  }

  #TG: Fill all the observed probabilities:
  tipvals.mat <- mapply(fill.observed.probabilities, tipvals.mat, tipvals, SIMPLIFY = FALSE)

  #~~~~~
  #TG: removing node.labels function
  remove.node.label <- function(tree) {
    tree$node.label <- NULL
    return(tree)
  }

  # Remove any potential node labels on the character tree to avoid an error from rerootingMethod():
  chartrees <- lapply(chartrees, remove.node.label)

  #~~~~~
  #TG: rerootingMethod lapply version (intakes 3 lists)
  apply.rerootingMethod <- function(character.list, chartrees, tipvals.mat, mymodels) {
    state.likelihoods <- rerootingMethod(chartrees[[character.list]], tipvals.mat[[character.list]], model = mymodels[[character.list]])$marginal.anc  
    #TG: Not entirely sure about that but using "phytools::rerootingMethod" rather than "rerootingMethod" might increase the computational speed (or not...)
  }

  # Get likelihoods for each state in taxa and ancestors:
  state.likelihoods <- lapply(character.list, apply.rerootingMethod, chartrees, tipvals.mat, mymodels)
          
  #TG: For this last bit (i.e. getting the max.lik and then getting the corresponding state), I've slightly changed the algorithm.
  #TG: This is to be able to incorporate the scaled likelihood threshold + optimising speed.
  #TG: The algorithm now goes:
    #1 - Check if any scale likelihood value is higher that the threshold
    #2 - If yes, order likelihood score in decreasing order and get the first one's name (i.e. the highest)
    #3 - If no, return NA

  #~~~~~
  #TG: Function for getting the most likely character state (or NA if scale likelihood < threshold)
  get.max.lik.state <- function(likelihood.estimates, uncertainty.threshold) {
    #Check if any likelihood estimation is higher than the threshold
    if(any(likelihood.estimates > uncertainty.threshold)) {
      #Check if there is only on highest likelihood score
      if(sort(likelihood.estimates, decreasing = TRUE)[1] != sort(likelihood.estimates, decreasing = TRUE)[2]) {
        #Get the name of the state with the highest score
        max.lik.state <- names(sort(likelihood.estimates, decreasing = TRUE)[1])
      } else {
        #Get all the states with equal likelihood
        max.lik.state <- paste(names(sort(likelihood.estimates, decreasing = TRUE))[which(sort(likelihood.estimates, decreasing = TRUE) == sort(likelihood.estimates, decreasing = TRUE)[1])], collapse = "&")
      }
    } else {
      #Return NA (no estimation above the threshold)
      max.lik.state <- NA
    }
    return(max.lik.state)
  }

  #~~~~~~
  #TG: function for applying get.max.lik.state to the list of state.likelihoods
  apply.get.max.lik.state <- function(state.likelihoods, uncertainty.threshold) {
    return(apply(state.likelihoods, 1, get.max.lik.state, uncertainty.threshold))
  }

  #TG: Calculating the maximum likelihood states
  max_lik_states <- lapply(state.likelihoods, apply.get.max.lik.state, uncertainty.threshold)
  #TG: Adding characters IDs
  names(max_lik_states) <- names(morph.matrix.character)

  #########
  # Combining all the characters together
  #########

  #TG: Add the NA.characters
  if(!is.null(NA.characters)) {
    max_lik_states <- c(max_lik_states, NA.characters)
  }

  #TG: Add the constant characters
  if(!is.null(constant.characters)) {
    max_lik_states <- c(max_lik_states, constant.characters)
  }

  #TG: Re order the characters
  max_lik_states <- max_lik_states[as.character(sort(as.numeric(names(max_lik_states))))]
  #TG: Transforming the list into a matrix
  anc.lik.matrix <- matrix(unlist(max_lik_states), byrow = FALSE, ncol = length(max_lik_states))
  #TG: Adding row names
  if(is.null(tree$node.label)) {
    rownames(anc.lik.matrix) <- c(tree$tip.label, (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
  } else {
    rownames(anc.lik.matrix) <- c(tree$tip.label, tree$node.label)
  }

  # Return completed discrete ancestral character estimation matrix:
  return(anc.lik.matrix)

}

# # Testing
# require(Claddis)
# require(testthat)
# morph.matrix <- Michaux1989
# tree <- rtree(nrow(Michaux1989$matrix))
# tree$tip.label <- rownames(Michaux1989$matrix)
# tree$root.time <- max(diag(vcv(tree)))
# estimate.allchars = TRUE
# estimate.tips = TRUE
# uncertainty.threshold = 0 # An uncertainty threshold value for the scaled likelihoods below wich estimations are replaced by NAs

# #Gauthier Testing
# #Make the tree
# tree <- rtree(nrow(Michaux1989$matrix))
# tree$tip.label <- rownames(Michaux1989$matrix)
# tree$root.time <- max(diag(vcv(tree)))
# #Run them
# fast_version <- system.time(fast_mat <- AncStateEstMatrix.fast(Michaux1989, tree, estimate.allchars = TRUE, estimate.tips = TRUE, uncertainty.threshold = 0))
# norm_version <- system.time(norm_mat <- AncStateEstMatrix(Michaux1989, tree, estimate.allchars = TRUE, estimate.tips = TRUE))


#Loading the packages
library(Claddis)
library(testthat)
#Some quick function to make the trees (random)
make.tree <- function(matrix) {
  tree <- rtree(nrow(matrix$matrix))
  tree$tip.label <- rownames(matrix$matrix)
  tree$root.time <- max(diag(vcv(tree)))
  return(tree)
}

#Testing with the Michaux 1989 matrix
set.seed(123)
time_base <- system.time(results_base <- AncStateEstMatrix(Michaux1989, tree = make.tree(Michaux1989), estimate.allchars = TRUE, estimate.tips = TRUE))
set.seed(123)
time_fast <- system.time(results_fast <- AncStateEstMatrix.fast(Michaux1989, tree = make.tree(Michaux1989), estimate.allchars = TRUE, estimate.tips = TRUE))
#What's the time difference?
time_base/time_fast

#Is it actually the same output?
test_that("Are both version giving the exact same results (no rounding nor random seeds)", {
  expect_equal(results_base, results_fast)
})

#Testing with the Gauthier 1986 matrix
set.seed(123)
time_base <- system.time(results_base <- AncStateEstMatrix(Gauthier1986, tree = make.tree(Gauthier1986), estimate.allchars = TRUE, estimate.tips = TRUE))
set.seed(123)
time_fast <- system.time(results_fast <- AncStateEstMatrix.fast(Gauthier1986, tree = make.tree(Gauthier1986), estimate.allchars = TRUE, estimate.tips = TRUE))
#What's the time difference?
time_base/time_fast

#Is it actually the same output?
test_that("Are both version giving the exact same results (no rounding nor random seeds)", {
  expect_equal(results_base, results_fast)
})

#Testing with a more complex (bigger) matrix
#Some random matrix
matrix <- matrix(sample(c("0", "1", "2", "0&1", "0&2", "1&2", "0&1&2"), 5000, replace=TRUE, prob=c(0.425, 0.42, 0.12, 0.01, 0.01, 0.01, 0.005))
  , nrow=50, ncol=100, dimnames=list(c(1:50)))
#Adding 25% of missing characters
matrix[sample(1:5000, 200)] <- NA

#ordering
ordering <- sample(c("unord", "ord"), 100, replace=TRUE, prob=c(0.9,0.1))

#weights
weights <- sample(c(1:3), 100, replace=TRUE, prob=c(0.85,0.1,0.05))

#Function for making a morph.matrix like object
make.nexus<-function(matrix, header, ordering, weights) {
    nexus<-list()
    nexus$header<-header
    nexus$matrix<-matrix
    nexus$ordering<-ordering
    nexus$weights<-weights
    nexus$max.vals<-as.numeric(apply(matrix, 2, max, na.rm=TRUE))
    nexus$min.vals<-as.numeric(apply(matrix, 2, min, na.rm=TRUE))
    return(nexus)
}

morph.matrix <- make.nexus(matrix, header="example", ordering, weights)

#Testing with the random matrix
set.seed(123)
time_base <- system.time(results_base <- AncStateEstMatrix(morph.matrix, tree = make.tree(morph.matrix), estimate.allchars = TRUE, estimate.tips = TRUE))
set.seed(123)
time_fast <- system.time(results_fast <- AncStateEstMatrix.fast(morph.matrix, tree = make.tree(morph.matrix), estimate.allchars = TRUE, estimate.tips = TRUE))
#What's the time difference?
time_base/time_fast

#Is it actually the same output?
test_that("Are both version giving the exact same results (no rounding nor random seeds)", {
  expect_equal(results_base, results_fast)
})

