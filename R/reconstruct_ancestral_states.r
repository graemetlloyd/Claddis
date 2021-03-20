#' Reconstruct ancestral states using parsimony
#'
#' @description
#'
#' Given a tree, discrete states for each tip, and either a specific stepmatrix or character type, returns the most parsimonious ancestral reconstruction(s).
#'
#' @param tree A tree (phylo object).
#' @param tip_states A labelled vector of tip states. These should be discrete and match the row and column headings in \code{stepmatrix}, with labels matching the tip labels of \code{tree}.
#' @param stepmatrix Either the character type (one of \code{"ordered"} or \code{"unordered"}) or a custom-specified step matrix. If the latter this must be square with identical row and column names correspnding to the values in \code{tip_states}. The diagonal must be zero, but off-diagonal values can be any non-negative real number. I.e., they needn't be integers, and the matrix does not need to be symmetric. Note that for transitions the rows are considered the "from" values, and the columns the "to" value. Thus a cost for the transition from state "0" to state "1" would be specified by \code{stepmatrix["0", "1"]}. (This is only relevant where matrices are asymmetric in tranistion costs.)
#' @param weight The character weight (defaults to one). Setting this value ensures that the length of the tree (number of steps) is calculated accurately. Note that if a custom stepmatrix is defined and the weight is already encoded into the costs then this should be left as one.
#' @param inapplicables_as_missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default).
#'
#' @details
#'
#' .
#'
#' # Uses Swofford and Maddison 1992 general solution as covers more bases (i.e., any stepmatrix plus works for polytomies) but is also slower.
#' # Rooted trees only as direction fucking matters and roots have meaning (stepmatrices that are asymmetric)
#' # Option to fix root state (polarize) somehow
#' # From is row, to is column (directional hence requires/assumes rooted trees)
#' # Record in manual that weighted parsimony is folly, but could use to chose between MPRs (maybe?)

#' \bold{Polymorphic characters}
#'
#' # Nixon and Davis (1991) show that ignoring (true) polymorphisms leads to inaccurate (under counted) tree lengths. Hence options here.
#' # Maddison adn Maddison (1987) suggested stepatrix solutionss that allow for polymorphic ancestors.

#' # MANUAL: generalized (stepmatrix) solution as accounts for more possibilities
#' # MANUAL: really only for rooted trees (unrooted possible using same methods (see Swofford and Maddison (1992), but unclear why you would want to do this with trees without direction of time (i.e., a root).
#' # MANUAL: works with multifurcations and doesn't require a fully bifurcating tree, but does assume polytomies are therefore "hard".
#' # MANUAL: Swofford and Maddison (1992; p210): "Both [ACCTRAN and DELTRAN] start by either assuming that the ancestral state is known or by choosing a state from the set of optimal assignments to the root node" - i.e., always an arbitrainess problem at root!
#' MISSING DATA (NA) IS DEALT WITH BY SETTING ALL TIP VALUES FOR A MISSING STATE TO ZERO (AS PER SWOFFORD AND MADDISON 1992).
#' UNCERTAIN STATES (E.G., 0/1) CONSTRAIN ANCESTRAL STATES SLIGHTLY MORE THAN NA (ASSUMING THEY EXCLUDE SOME STATES. IMPLEMENTATION SAME AS SWOFFORD AND MADDISON (1992).
#' CALCULATES TREE LENGTHS!
#'
#' ASYMMETRIC TRANSITIONS REQUIRE YOU TO KNOW POLARITY A PRIORI!
#' HOW TO TREAT POLYTOMIES? MANY ISSUES RAISED IN SWOFFORD AND MADDISON (1992). ALTERNATIVE APPROACH IN LIKELIHOOD SUGGESTED BY REVELL.
#'
#' Users may also wish to refer to the more complex whole-matrix, likelihood-based function \link{estimate_ancestral_states}. Although, note that eventually parsimony will simply be an option to that function.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Maddison, D. R. and Maddison, W. P., 2000. \emph{MacClade 4}. Sinauer Associates Incoroporated, Sunderland. 492pp.
#'
#' Maddison, W. P. and Maddison, D. R., 1987. \emph{MacClade 2.1}. Distributed by the authors, Cambridge.
#'
#' Nixon, K. C. and Davis, J. I., 1991. Polymorphic taxa, missing values and cladistic analysis. \emph{Cladistics}, \bold{7}, 233-241.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. \emph{In} R. L. Mayden (ed.) Systematics, Historical Ecology, and North American Freshwater Fishes. Stanford University Press, Stanford, p187-223.
#'
#' @return
#'
#' A list with multiple components, including:
#'
#' \item{length}{The tree length (number of steps).}
#' \item{most_parsimonious_reconstructions}{A matrix where rows correspond to \emph{all} nodes (i.e., terminal and internal) in the order numbered by \code{ape}, and columns correspond to every unique most parsimonious reconstruction. I.e., if there is only one most parsimonious reconstruction there will be only one column.}
#' \item{input_tree}{The input topology used.}
#'
#' @seealso
#'
#' \link{estimate_ancestral_states}
#'
#' @examples
#'
#' # Set up the example from Swofford and Maddison 1992:
#' tree <- ape::read.tree(text = "((C,(A,B)),(D,E));")
#' tip_states <- c(A = 1, B = 2, C = 2, D = 1, E = 0)
#' stepmatrix <- matrix(
#'   data = c(0, 1, 2, 1, 0, 1, 2, 1, 0),
#'   nrow = 3,
#'   ncol = 3,
#'   dimnames = list(0:2, 0:2)
#' )
#'
#' # Look at stepmatrix to confirm it is symmetric and ordered:
#' stepmatrix
#'
#' # Reconstruct ancestral states:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )$most_parsimonious_reconstructions
#'
#' # Plot all most parsimonious reconstructions on the tree:
#' par(mfrow = c(1, ncol(x = node_estimates)))
#' for(i in 1:ncol(x = node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat using the simpler "ordered" option:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )$most_parsimonious_reconstructions
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(x = node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat but considering the character "unordered":
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "unordered"
#' )$most_parsimonious_reconstructions
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(x = node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat using a stepmatrix that penalises gains:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = matrix(
#'     data = c(0, 1, 2, 999, 0, 1, 999, 999, 0),
#'     nrow = 3,
#'     ncol = 3,
#'     dimnames = list(0:2, 0:2)
#'   )
#' )$most_parsimonious_reconstructions
#' par(mfrow = c(1, ncol(x = node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat using a stepmatrix that penalises losses:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = matrix(
#'     data = c(0, 999, 999, 1, 0, 999, 2, 1, 0),
#'     nrow = 3,
#'     ncol = 3,
#'     dimnames = list(0:2, 0:2)
#'   )
#' )$most_parsimonious_reconstructions
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Ordered with polytomies:
#' tree <- ape::read.tree(text = "((C,A,B),(D,E));") # Two non-root internal
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )$most_parsimonious_reconstructions
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(tree, show.tip.label = FALSE)
#'   tiplabels(tip_states[tree$tip.label], cex = 2)
#'   nodelabels(node_estimates[6:8, i], cex = 2)
#' }
#' tree <- ape::read.tree(text = "((C,A,B),D,E);") # One non-root internal
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )$most_parsimonious_reconstructions
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(tree, show.tip.label = FALSE)
#'   tiplabels(tip_states[tree$tip.label], cex = 2)
#'   nodelabels(node_estimates[6:7, i], cex = 2)
#' }
#' tree <- ape::read.tree(text = "(C,A,B,D,E);") # Star tree
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )$most_parsimonious_reconstructions
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(tree, show.tip.label = FALSE)
#'   tiplabels(tip_states[tree$tip.label], cex = 2)
#'   nodelabels(node_estimates[6, i], cex = 2)
#' }
#'
#' # Get tree length for an entire matrix of characters:
#' tree <- ape::read.tree(
#'   text = "(Outgroup,(Ornithischia,(Sauropodomorpha,(Ceratosauria,Procompsognathus,
#'   Liliensternus,(Carnosauria,(Ornithmimidae,Saurornitholestes,Hulsanpes,(Coelurus,
#'   Elmisauridae,(Compsognathus,(Ornitholestes,Microvenator,Caenagnathidae,
#'   (Deinonychosauria,Avialae))))))))));"
#' )
#' tip_state_matrix <- gauthier_1986$matrix_1$matrix
#' sum(x = apply(
#'   X = tip_state_matrix,
#'   MARGIN = 2, FUN = function(x) {
#'     reconstruct_ancestral_states(
#'       tree = tree,
#'       tip_states = x,
#'       stepmatrix = "unordered",
#'       weight = 1)$length
#'   }
#' ))
#'
#' @export reconstruct_ancestral_states
reconstruct_ancestral_states <- function(tree, tip_states, stepmatrix, weight = 1, inapplicables_as_missing = FALSE) {
  
  # WHOLE THING SHOULD BE REFACTORED TO DEAL WITH ENTIRE MATRIX (THE OBVIOUS END USE)
  # ALSO ALLOWS MOVING MANY CHECKS TO A SINGLE PASS OF MATRIX AND HENCE FASTER
  
  # TO DO
  # - Allow for true polymorphisms (maybe conditionals for counting changes?) This is a hard problem!
  # - How to allow for missing or uncertainty? This bit should be easy as set all states to zero, or all uncertain states to zero at tips.
  # - All (format?)/ACCTRAN/DELTRAN/(Random?)/Branch lengths as weights to help collapse?/Uncertainity (e.g., 0/1)
  # - Missing and inapplicables shpuld probably be assigned to nodes on a third pass "down" the tree (tips to roots) such that any all-descendants set gets assigned an NA/"" too. Otherwise sets false certainity by "bleeding" an ancestral value "up" the tree.
  # - Tip states as entered need preserving separately as some modification will need to happen to the variable if certain options are chosen.
  # - Make options match estimate_ancestral_states:
  #   - @param estimate_all_nodes Logical that allows the user to make estimates for all ancestral values. The default (\code{FALSE}) will only make estimates for nodes that link coded terminals (recommended).
  #   - @param estimate_tip_values Logical that allows the user to make estimates for tip values. The default (\code{FALSE}) will only makes estimates for internal nodes (recommended).
  #   - @param polymorphism_behaviour One of either "equalp" or "treatasmissing".
  #   - @param uncertainty_behaviour One of either "equalp" or "treatasmissing".

  # CHECKS TO WRITE
  # - states are discrete
  # - states in tip states are all available in stepmatrix and vice versa
  # - stepmatrix values are zero or positive (don't have to be integers?)
  # - Check there are enough states to do something with!
  # - Put checks at higher level to avoid repeating for every character in a matrix (slow?)
  
  # Reorder tip states (1 to N) and store an unmodified version:
  input_tip_states <- tip_states <- tip_states[tree$tip.label]
  
  # Store pristine input tree:
  input_tree <- tree
  
  # If there are no branch durations set these as all one:
  if(is.null(x = tree$edge.length[1])) tree$edge.length <- rep(x = 1, length.out = nrow(x = tree$edge))
  
  # If treating inapplicables as missing then replace any inapplicable tip state with NA:
  if(inapplicables_as_missing) tip_states[tip_states == ""] <- NA
  
  # If stepmatrix is not already specified as a matrix (i.e., it is simply "ordered", "unordered" etc.):
  if(!is.matrix(x = stepmatrix)) {
    
    # NEED TO DEAL WITH POLYMORPHISMS IF NOT BEING TREATED AS UNCERTAINTIES
    # NEED TO SET MORE OPTIONS LIKE DOLLO, CAMIN-SOKAL AND STRATIGRAPHIC
    
    # First check stepmatrix is of a valid type and stop and warn user if not:
    if(length(x = setdiff(x = stepmatrix, y = c("ordered", "unordered"))) > 0) stop("If not a specific matrix, then stepmatrix must be one of \"ordered\" or \"unordered\".")
    
    # Get numeric tip state values:
    tip_state_numbers <- as.numeric(x = unlist(x = strsplit(x = as.character(x = tip_states), split = "&|/")))
    
    # Get the range (min-max) of tip values:
    tip_state_range <- range(x = tip_state_numbers[!is.na(x = tip_state_numbers)])
    
    # Get the size (n rows and columns) required for the stepmatrix:
    stepmatrix_size <- diff(x = tip_state_range) + 1
    
    # Create base stepmatrix (to be altered further):
    base_stepmatrix <- matrix(data = NA, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(tip_state_range[1]:tip_state_range[2], tip_state_range[1]:tip_state_range[2]))
    
    # Add diagonal of zero as this should always be true regardless of character type:
    diag(x = base_stepmatrix) <- 0
    
    # Case if character is ordered:
    if(stepmatrix == "ordered") {
      
      # Populate upper triangle (cost is difference between states):
      base_stepmatrix[upper.tri(x = base_stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1))
      
      # Populate lower triangle (cost is difference between states):
      base_stepmatrix[lower.tri(x = base_stepmatrix)] <- unlist(x = lapply(X = as.list(x = (stepmatrix_size - 1):1), FUN = function(steps) 1:steps))
      
    }
    
    # Case if character is unordered all off-diagonal values cost one step):
    if(stepmatrix == "unordered") base_stepmatrix[is.na(x = base_stepmatrix)] <- 1
    
    # Overwrite stepmatrix with the fully formed stepmatrix:
    stepmatrix <- base_stepmatrix
    
  }
  
  # Apply character weight by multiplying through stepmatrix:
  stepmatrix <- stepmatrix
  
  # Reformat tip states as character vector:
  tip_states <- as.character(x = tip_states)
  
  # Re-add names to tip states:
  names(tip_states) <- tree$tip.label
  
  # Establish number of tips:
  tip_count <- ape::Ntip(phy = tree)
  
  # Establish total number of nodes (terminal and internal):
  node_count <- tip_count + tree$Nnode
  
  # Initialise node_values matrix:
  node_values <- matrix(data = NA, nrow = node_count, ncol = ncol(x = stepmatrix), dimnames = list(c(), colnames(stepmatrix)))
  
  # Begin by inserting Inf as default tip value:
  node_values[1:tip_count, ] <- Inf
  
  # Populate tip values (zeroes for assigned state(s):
  for(i in 1:tip_count) node_values[i, match(x = strsplit(x = tip_states[i], split = "/")[[1]], table = colnames(x = node_values))] <- 0

  # If any tips are missing then set all states for those as zero:
  if(any(x = is.na(x = tip_states))) node_values[which(x = is.na(x = tip_states)), ] <- 0
  
  # First pass (traverse tree from tips to root):
  for(needle in (tip_count + node_count - tip_count):(tip_count + 1)) {
    
    # Find decsendants of current node:
    descendants <- tree$edge[tree$edge[, 1] == needle, 2]
    
    # Calculate and store new node values:
    node_values[needle, ] <- unlist(x = lapply(X = as.list(x = colnames(x = node_values)), FUN = function(fromstate) sum(x = unlist(x = lapply(X = as.list(x = descendants), FUN = function(descendant) min(x = node_values[descendant, ] + stepmatrix[fromstate, ]))))))

  }
  
  # ADD OPTION TO FIX ROOT HERE (NAH, ALLOW OPTION TO FIX *ANY* NODE (OR COMBO OF NODES)
  # CHECK OPTION DOESN'T FORCE INFINITY!
  # CHECK Inf VALUES IN STEP MATRICES DO NOT LEAD TO ALL-Inf VALUES AT NODES!
  
  # Store tree length:
  tree_length <- min(x = node_values[tip_count + 1, ]) * weight
  
  # STOP AFTER HERE IF ONLY WANT TREE LENGTH?
  
  # WILL NEED TO MODIFY BELOW TO DEAL WITH UNCERTAINTIES AND POLYMORPHISMS AT TIPS
  
  # Build node estimates from single state results:
  node_estimates <- matrix(data = apply(X = node_values, MARGIN = 1, FUN = function(x) {x <- names(x = x[x == min(x = x)]); ifelse(test = length(x = x) == 1, x, NA)}), ncol = 1, nrow = node_count)
  
  # Update tip states with their original values:
  node_estimates[1:tip_count] <- tip_states
  
  ### ABOVE NEEDS TO BE INPUT TIP STATES NOT MODIFED VERSION! BUT LATER FOR OUTPUT? NEED TO RECORD SOMEWHERE WHAT HAS HAPPENED THUGH SO CAN PASS TO A CHARACTER MAP FUNCTION AND TREAT EVERYTHING CORRECTLY.
  
  # Only continue if there is any ambiguity at internal nodes:
  if(any(x = is.na(x = node_estimates[(tip_count + 1):node_count]))) {
    
    # Find all possible root states:
    possible_root_states <- colnames(x = node_values)[node_values[tip_count + 1, ] == min(x = node_values[tip_count + 1, ])]
    
    # Make new node estimates with every possible root state:
    node_estimates <- do.call(what = cbind, args = lapply(X = as.list(x = possible_root_states), FUN = function(x) {y <- node_estimates; y[tip_count + 1, ] <- x; y}))
      
    # For each internal node above the root to the tips:
    for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
      
      # Establish ancestor of current node:
      ancestor_node <- tree$edge[tree$edge[, 2] == needle, 1]
      
      # Reformat node_estimates as a list to enable lapply:
      node_estimates <- split(x = node_estimates, f = col(x = node_estimates))
      
      # Permute all possible values and reformat as matrix:
      node_estimates <- do.call(what = cbind, args = lapply(X = node_estimates, FUN = function(x) {
        
        # Get updated tree lengths for current node as per Swofford and Maddison 1992 second pass:
        updated_tree_lengths <- stepmatrix[x[ancestor_node], ] + node_values[needle, ]
        
        # Store all possible most parsimonious state(s) for current node:
        possible_states <- names(x = updated_tree_lengths[updated_tree_lengths == min(x = updated_tree_lengths)])
        
        # Store total number of possible states:
        n_states <- length(x = possible_states)
        
        # Create new node estimates to store possible state(s):
        new_estimates <- matrix(data = rep(x = x, times = n_states), ncol = n_states)
        
        # Add possible state(s)
        new_estimates[needle, ] <- possible_states
        
        # Return new estimates only:
        new_estimates
        
      }))
      
    }
    
  }
  
  
  
  
  
  
  
  
  
  make_all_polymorphisms <- function(single_states) unlist(x = lapply(X = as.list(x = 1:length(x = single_states)), FUN = function(x) apply(X = combn(x = sort(x = single_states), m = x), MARGIN = 2, FUN = function(y) paste(x = y, collapse = "&"))))
  
  make_stepmatrix <- function(all_states, method = c("hypercube", "manhattan")) {
    
    # Coordinate space:
    #"hypercube"
    #"hypersphere"
    #"simplex"
    
    # Distance method:
    #"manhattan"
    #"euclidean"
    #"great_circle"
    
    # Check data are not too big (>= 2^14 states) and stop and warn user if so:
    if(length(x = all_states) >= 16384) stop("Stepmatrix would be too large. Use fewer states.")
    
    # Get single states:
    single_states <- all_states[grep(pattern = "&", x = all_states, invert = TRUE)]
    
    # Create coordinate matrix and initialise with zeroes:
    state_presence_matrix <- matrix(data = 0, nrow = length(x = single_states), ncol = length(x = all_states), dimnames = list(single_states, all_states))
    
    # If using the hypercube coordinate space assign coordinates accordingly:
    if(method[1] == "hypercube") for(i in 1:ncol(x = state_presence_matrix)) state_presence_matrix[strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]], i] <- 1
    
    # If using the hypersphere or simplex coordinate space:
    if(method[1] == "hypersphere" || method[1] == "simplex") {
      
      # For each coding:
      for(i in 1:ncol(x = state_presence_matrix)) {
        
        # Isolate components of polymorphism:
        components <- strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]]
        
        # If using the hypersphere coordinate space then apply coordinates accordingly (the square root of 1/N states in polymorphism on each axis state is present):
        if(method[1] == "hypersphere") state_presence_matrix[components, i] <- sqrt(x = 1 / length(x = components))
        
        # If using the simplex coordinate space then apply coordinates accordingly (1/N states in polymorphism on each axis state is present):
        if(method[1] == "simplex") state_presence_matrix[components, i] <- 1 / length(x = components)
        
      }
      
    }
    
    # If using a manhattan or euclidean distance, calculate distance directly from coordinate-space:
    if(method[2] == "euclidean" || method[2] == "manhattan") stepmatrix <- as.matrix(x = dist(x = t(x = state_presence_matrix), method = method[2], diag = TRUE, upper = TRUE))
    
    # If using a great circle distance:
    if(method[2] == "great_circle") {
      
      # Start by calculating the euclidean distances between each point:
      stepmatrix <- as.matrix(x = dist(x = t(x = state_presence_matrix), method = "euclidean", diag = TRUE, upper = TRUE))
      
      # Treating distances as the chord length of a circle of radius one transform them to the arc length for the same:
      stepmatrix <- 2 * asin(x = stepmatrix / 2)
      
    }
    
    # Return a stepmatrix rescaled such that single state distances (e.g., 0 to 1) are one (i.e., a normal unordered character):
    stepmatrix / stepmatrix[1, 2]
    
    ### ABOVE WILL NEED TO BE MODIFIED IF USING ASYMMETRIC CHARACTERS LIKE DOLLO OR CAMINSOKAL OR STRATIGRAPHY
    
  }
  
  # KEY THING HERE IS ALLOWS POLYMORPHISMS AT INTERNAL NODES!
  



  
  # ADD SWOFFORD REF TO DESCRIPTION FILE (PROLLY NO DOI!)
  
  
  
  
  
  
  
  
  
  # Return compiled output:
  list(length = tree_length, most_parsimonious_reconstructions = node_estimates, input_tree = input_tree)

  # OUTPUT:
  # - Most of this needs to be put in different functions!
  # - Ancestral states
  # - Ambiguities?
  # - Stepmatrix used
  # - Tree used
  # - Tip states used
  # - Tree length
  # - How root value was chosen (arbitrary forced, unambiguous?)
  # - Algorithm (parsimony, ML, etc.)
  # - CI? (Character vs matrix - can they just be added? - I think not)
  # - RI?
  # - N reversals
  # - N parallelisms
  # - Distortion index for each MPR
  # - Character map?
  # - MAKE THIS A CLASS TOO!
  
}

# WRITE SOME KIND OF CHECK FOR SELF-CONSISTENCY IN STEP MATRICES. I.E., THE DIRECT PATH IS ALWAYS THE SAME LENGTH OR SHORTER THAN AN INDIRECT ONE. E.G., If 0 -> 2 is 4 steps, but 0->1->2 is only 3 steps. (Citation is MacClade 4 manual!)
# OPTION TO TREAT POLYORPHISMS AS MISSING (NOT JUST INAPPLICABLES).
# OPTION TO CONSTRAIN PARTICULAR NODES AND/OR MAKE TAXA ANCESTORS (MIGHT BE TRICKY AND A DOWN THE LINE ADDITION!) HMMM. OR HAVE CONDITIONAL THAT SETS INFINITE VAUES TO THE NON PRESENT STATES AT THIS NODE?
# STRATIGRAPHIC AS A CHARACTER TYPE ALONG WITH DOLLO AND IRRERVERSIBLE (CAMIN-SOKAL)
# WHAT IF CHARACTER IS WEIGHTED ZERO! NEED A CONDITIONAL TO DEAL WITH THIS AS STEPMATRIX CANNOT BE ALL ZEROES! MAYBE APPLY WEIGHT AT END AS MULTIPLE OF TREE LENGTH? PROBABLY SIMPLEST SOLUTION.
# MORE PUTPUT OR SUMMARY OPTIONS COULD BE N CHANGES ON EACH BRANCH ACROSS MPRS, MEAN OF SAME, MIN AND MAX OF SAME. SIMILAR FOR CHANGE TYPES I.E., STEP MATRX BUT ACTUALLY FREQUENCY OF CHANGES FOR EACH TRANSITION.
# AS RATE TENDS TOWARD INFINITY THEN RECONSTRUCTION AT A NODE TENDS TOWARDS EQUAL FREQUENCY OF EACH STATE (I.E., MAXIMAL UNCERTAINTY).
# ALLOW TREES N STEPS LONGER SOMEHOW? OR RATHER MPRs ONE STEP LONGER, OR SOME OTHER VALUE OF USERS CHOOSING.
# "This illustrates the fact that ACCTRAN and DELTRAN do not always choose a single one of the most parsimonious reconstructions." MacClade 4 manual, page 99).
# "Note that MacClade, unlike PAUP*, does not choose the lowest-valued state at the root to begin these processes. Thus MacClade's ACCTRAN and DELTRAN may not fully resolve ambiguity in the ancestral state reconstruction."
# INTERMEDIATES ARE GONNA MATTER IF MOVING TO CHARACTER MAPS. E.G., IF ONLY 0 AND 2 ARE SAMPLED BUT A CARACTER IS ORDERED THEN THERE ARE TWO CHANGES ALONG THE BRANCH NOT ONE TWO-STEP CHANGE.
# POLYMORPHISM APPROACH SHOULDN'T CONFOUND THE "NORMAL" BEHAVIOUR OF AN ORDERED OR UNORDERED CHARACTER.
# FOR SPEED KEEP STEPMATRIX OUTSIDE OF ASR FUNCTION AS WILL OFTEN REUSE SAME ONE.


# ACCTRAN
# DELTRAN
# MINF (Swofford and Maddison 1987) - I THINK NOT!
# MINSTATE *shrug emoi*
# MAXSTATE *shrug emoi*
# WEIGHTED BY BRANCH DURATION?


