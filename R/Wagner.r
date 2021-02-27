Ntips <- 10
Nstates <- 5
tree <- rtree(Ntips)
tip_states <- sample(0:(Nstates - 1), Ntips, replace = TRUE)
names(tip_states) <- tree$tip.label
#tree <- ape::read.tree(text="(A:1,(B:1,C:1):1);")
#tip_states <- c(A = 0, B = 1, C = 2)

# Kitching example:
#tree <- ape::read.tree(text="(A:1,(B:1,((C:1,D:1):1,(E:1,F:1):1):1):1);")
#tip_states <- c(A = 1, B = 0, C = 1, D = 2, E = 2, F = 4)


# one for length one for states?
# how to do missing?
# how to do inapplicable
# how to do polymorphism (tricky as must choose between "single acqusiition of additional states or one for each state)
# how to do uncertainty (as a parsimony algorithm should be minimum possible)
# generalise to non-binary trees immediately?
# wtf happens at root?
# optionally incorporate edge lengths as weighting somehow?
# do acctrana dn deltran require prior knowledge of whether states 0 or 1 would be more primitive?
estimate_ancestral_states <- function(tree, tip_states, method = "farris") {
  
  # CHECK CHARACTER IS VARIABLE -> FAST ANSWER IF NOT
  # CHECK CHARACTER IS NOT AUTAPOMORPHIC -> FAST ANSWER IF NOT
  # CHECK CHARACTER IS CODED FOR EVERY TIP
  # CHECK TREE HAS AT LEAST (TWO?) TIPS
  # OPTION TO ESTIMATE TIPS (COLLAPSE UNCERTAINTIES AND/OR POLYMORPHISMS)?
  
  # Arbtrary rooting doesn't seem to guarantee singleton states because chosen taxa could be polymorphic
  
  # Check tree is bifurcating and stop and warn user if not:
  if(!ape::is.binary(phy = tree)) stop("Function only works for fully bifurcating trees.")
  
  # Subfunction to get farris intervals (sensu Swofford and Maddison 1987):
  get_farris_intervals <- function(i, j) {
    
    # First form intersection of i and j:
    ij_intersection <- intersect(x = i, y = j)
    
    # If Farris's first rule is in effect (intersection is non-empty) return intersection:
    if(length(x = ij_intersection) > 0) return(ij_intersection)
    
    # If intersection is empty, then use equation 3 of Swofford and Maddison (1987):
    if(length(x = ij_intersection) == 0) {
      
      # Get min-maxes:
      min_i <- min(x = i)
      min_j <- min(x = j)
      max_i <- max(x = i)
      max_j <- max(x = j)
      
      # Get y and z:
      y <- median(x = c(min_i, max_i, min_j))
      z <- median(x = c(max_i, min_j, max_j))
      
      # Return Farris's second rule value:
      return(min(x = c(y, z)):max(x = c(y, z)))
      
    }
    
  }
  
  # Establish number of tips:
  tip_count <- ape::Ntip(phy = tree)
  
  # Establish total number of nodes (terminal and internal):
  node_count <- tip_count + tree$Nnode
  
  # Set up list to store values for every node:
  node_values <- as.list(x = rep(x = NA, times = node_count))
  
  # Populate tip values immediately:
  node_values[1:tip_count] <- tip_states[tree$tip.label]
  
  # Subfunction to do a Farris first pass:
  farris_first_pass <- function(tree, node_values, tip_count, node_count) {
    
    # Now for each internal node traversing tree from tips to root:
    for(needle in (tip_count + node_count - tip_count):(tip_count + 1)) {
      
      # Get immediate descendant values for current node:
      descendant_values <- node_values[tree$edge[tree$edge[, 1] == needle, 2]]
      
      # Get Farris intervals for current node:
      node_values[[needle]] <- get_farris_intervals(i = descendant_values[[1]], j = descendant_values[[2]])
      
    }
    
    # Return updated node values:
    node_values
    
  }

  # If using the Farris method:
  if(method == "farris") {
    
    # Do a first pass of the tree as per Farris rule 1:
    node_values <- farris_first_pass(tree = tree, node_values = node_values, tip_count = tip_count, node_count = node_count)
    
    # Only need to do second pass if there are any remaining ambiguiuous nodes:
    if(any(x = lapply(X = node_values[(tip_count + 1):node_count], FUN = length) > 1)) {
      
      # Now for each internal node (excluding root) traverse tree from root to tips:
      for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
        
        # Only something to do if current node state is ambiguous:
        if(length(x = node_values[[needle]]) > 1) {
          
          # Get immediate ancestor value(s) for current node:
          ancestor_values <- node_values[[tree$edge[tree$edge[, 2] == needle, 1]]]
          
          # If ancestral value is a singleton then take median of this and limits:
          if(length(x = ancestor_values) == 1) node_values[[needle]] <- median(x = c(ancestor_values, range(x = node_values[[needle]])))
          
        }
        
      }
      
    }
    
  }
  
  # If using the Farris method:
  if(method == "sm87") {
    
    # Do a first pass of the tree as per Farris rule 1:
    node_values <- farris_first_pass(tree = tree, node_values = node_values, tip_count = tip_count, node_count = node_count)

  }

  # Return node values:
  node_values
  
  # Plot output for visual inspaection when testing code:
  ape::plot.phylo(x = tree, show.tip.label = FALSE)
  ape::tiplabels(text = unlist(x = lapply(X = node_values, FUN = function(x) paste(x = x, collapse = ",")))[1:tip_count])
  ape::nodelabels(text = unlist(x = lapply(X = node_values, FUN = function(x) paste(x = x, collapse = ",")))[(tip_count + 1):(tip_count + tree$Nnode)])
  
}

estimate_ancestral_states(tree, tip_states, method = "sm87")



get_farris_intervals(i = 0, j = c(0,1))
get_farris_intervals(i = c(0:2), j = c(4:5))

i = 0:1
j = 2:3
k = 4:5

Reduce(
  f = intersect,
  x = list(
    get_farris_intervals(get_farris_intervals(j, k), i),
    get_farris_intervals(get_farris_intervals(i, k), j),
    get_farris_intervals(get_farris_intervals(i, j), k)
  )
)

i = 0:1
j = 2:3
k = 7:11
l = 5

Reduce(
  f = intersect,
  x = list(
    get_farris_intervals(get_farris_intervals(get_farris_intervals(l, k), j), i),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(i, l), j), k),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(i, k), j), l),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(l, i), k), j),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(j, l), k), i),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(j, i), k), l),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(j, k), l), i),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(i, k), l), j),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(i, j), l), k),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(l, k), i), j),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(j, l), i), k),
    get_farris_intervals(get_farris_intervals(get_farris_intervals(j, k), i), l)
  )
)


# Uses Swofford and Maddison 1992 general solution as covers more base's but is also slower.
# Rooted tres only as direction fucking matters and roots have meaning (stepmatrices that are asymmetric)
# Option to fix root state (polarize) somehow
# From is row, to is column

# ADD CI AND RI FUNCTIONS (AND OTHERS?)

# Swofford and Maddison 1992 example:
tree <- ape::read.tree(text = "((C,(A,B)),(D,E));")
tip_states <- c(A = 1, B = 2, C = 2, D = 1, E = 0)
stepmatrix <- matrix(data = c(0,1,2,1,0,1,2,1,0), nrow = 3, ncol = 3, dimnames = list(0:2, 0:2))



tree <- ape::read.tree(text = "((C,((A,B),(I,J))),(D,((E,F),(G,H))));")
tip_states <- c(A = 1, B = 2, C = 2, D = 1, E = 0, F = 0, G = 0, H = 2, I = 1, J = 1)
stepmatrix <- matrix(data = c(0,1,2,1,0,1,2,1,0), nrow = 3, ncol = 3, dimnames = list(0:2, 0:2))


# MANUAL: generlzied solution as accounts for more possibilities
# MANUAL: really only for rooted trees (unrooted possible using same methods (see Swofford and Maddison(1992), but unclear why you would want to do this with trees without direction of time (i.e., a root).
# MANUAL: works with multifirctaions and doesn't reuire a fully bifurcating tree, but does assume polytomies are therefore "hard".
# MANUAL: Swofford and Maddison (1992; p210): "Both [ACCTRAN and DELTRAN] start by either assuming that the ancestral state is known or by choosing a state from the set of optimal assignments to the root node" - i.e., always an arbitrainess problem at root!
parsimony_ancestors <- function(tree, tip_states, stepmatrix) {
  
  # TO DO
  # - Allow for polymorphisms (maybe conditionals for counting changes?) This is a hard problem!
  # - How to allow for missing or uncertainty? This bit should be easy as set all states to zero, or all uncertain states to zero at tips.
  # - ACCTRAN/DELTRAN/(Random?)/Branch lengths as weights to help collapse?/Uncertainity (e.g., 0/1)
  # - Character weights? (Just multiply through stepmatrix by weight?)
  
  # CHECKS TO WRITE
  # - states are discrete
  # - states in tip states are all available in stepmatrix and vice versa
  # - stepmatrix values are zero or positive (don't have to be integers?)
  
  # If stepmatrix is not already specified as a matrix (i.e., it is simply "ordered", "unordrered" etc.):
  if(!is.matrix(x = stepmatrix)) {
    
    # NEED TO SET MORE OPTIONS LIKE DOLLO, CAMIN-SOKAL AND....?
    
    # First check stepmatrix is of a vaid type and stop wand warn user if not:
    if(length(x = setdiff(x = stepmatrix, y = c("ordered", "unordered"))) > 0) stop("If not a specific matrix, then stepmatrix must be one of \"ordered\" or \"unordered\".")
    
    # Get tip states as numerics:
    tip_state_numbers <- as.numeric(x = tip_states)
    
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
  
  # Reorder tip states:
  tip_states <- tip_states[tree$tip.label]
  
  # Establish number of tips:
  tip_count <- ape::Ntip(phy = tree)
  
  # Establish total number of nodes (terminal and internal):
  node_count <- tip_count + tree$Nnode
  
  # Initialise node_values matrix:
  node_values <- matrix(data = NA, nrow = node_count, ncol = ncol(x = stepmatrix), dimnames = list(c(), colnames(stepmatrix)))
  
  # Begin by inserting Inf as default tip value:
  node_values[1:tip_count, ] <- Inf
  
  # Populate tip values (zeroes for assigned state(s):
  for(i in 1:tip_count) node_values[i, colnames(x = node_values) == tip_states[i]] <- 0
  
  # First pass (traverse tree from tips to root):
  for(needle in (tip_count + node_count - tip_count):(tip_count + 1)) {
    
    # Find decsendants of current node:
    descendants <- tree$edge[tree$edge[, 1] == needle, 2]
    
    # Calculate and store new node values:
    node_values[needle, ] <- unlist(x = lapply(X = as.list(x = colnames(x = node_values)), FUN = function(fromstate) sum(x = unlist(x = lapply(X = as.list(x = descendants), FUN = function(descendant) min(x = node_values[descendant, ] + stepmatrix[fromstate, ]))))))
    
  }
  
  # ADD OPTION TO FIX ROOT HERE
  # CHECK OPTION DOESN'T FORCE INFINITY!
  
  tree_length <- min(x = node_values[tip_count + 1, ])
  
  # STOP AFTER HERE IF ONLY WANT TREE LENGTH
  
  # Return
  #tree_length
  node_values
  
  
  # MAYBE NEED TO DO ALL MPRS FIRST THEN DO ACCTRAN/DELTRAN SOMEHOW?
  
  # Build matrix of node estimates:
  node_estimates <- matrix(data = apply(X = node_values, MARGIN = 1, FUN = function(x) {x <- names(x = x[x == min(x = x)]); ifelse(test = length(x = x) == 1, x, NA)}), ncol = 1, nrow = node_count)
  node_estimates[(tip_count + 1):node_count, ] <- NA
  
  # Only continue if some estimates are ambiguous:
  if(any(x = is.na(x = node_estimates))) {
    
    # Special case if root is ambiguous:
    if(is.na(x = node_estimates[tip_count + 1])) {
      
      # Find all possible root states:
      possible_root_states <- colnames(x = node_values)[node_values[tip_count + 1, ] == min(x = node_values[tip_count + 1, ])]
      
      # Make new node estimates with every possible root state:
      node_estimates <- do.call(what = cbind, args = lapply(X = as.list(x = possible_root_states), FUN = function(x) {y <- node_estimates; y[tip_count + 1, ] <- x; y}))
      
    }
    
    # Only continue if some estimates are ambiguous:
    if(any(x = is.na(x = node_estimates))) {
      
      # For each remaining ambiguous node:
      for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
        
        # Permute new tree lengths given ancestral state(s):
        permuted_new_lengths <- stepmatrix[node_estimates[tree$edge[tree$edge[, 2] == needle, 1], ], ] + matrix(data = rep(x = node_values[needle, ], times = ncol(x = node_estimates)), nrow = ncol(x = node_estimates), byrow = TRUE)
        
        # SOMETHING GOES WRONG HERE WITH TRANSPOSITIONS OF MATRIX I DO NOT WANT IT TO DO. GRRRR.
        # MIGHT HAVE TO MAKE NODE ESTIMATES A LIST AND REBUILD AROUND THAT.
        
        # Update node estimates with all possible most parsimonious state(s):
        new_node_estimates <- apply(X = cbind(1:nrow(x = permuted_new_lengths), permuted_new_lengths), MARGIN = 1, FUN = function(x) {column_number <- x[1]; x <- x[-1]; possible_states <- names(x = x[x == min(x = x)]); n_states <- length(x = possible_states); new_node_estimates <- matrix(data = rep(x = node_estimates[, column_number], times = n_states), ncol = n_states); new_node_estimates[needle, ] <- possible_states; new_node_estimates})
        
        # If output turned out to be a list coerce back to a matrix:
        if(is.list(x = new_node_estimates)) new_node_estimates <- do.call(what = cbind, args = new_node_estimates)
        
        # Update node estimates:
        node_estimates <- new_node_estimates
        
      }
      
    }

  }
  
  node_estimates <- list(c(unname(obj = tip_states[tree$tip.label]), rep(x = NA, times = node_count - tip_count)))
  
  # Find all possible root states:
  possible_root_states <- colnames(x = node_values)[node_values[tip_count + 1, ] == min(x = node_values[tip_count + 1, ])]
  #possible_root_states <- c("1", "2")
  # Make new node estimates with every possible root state:
  node_estimates <- lapply(X = as.list(x = possible_root_states), FUN = function(x) {new_node_estimates <- node_estimates[[1]]; new_node_estimates[tip_count + 1] <- x; new_node_estimates})
  
  # For each internal node above the root to the tips:
  for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
  
    # Permute new tree lengths given ancestral state(s):
    permuted_new_lengths <- stepmatrix[unlist(x = lapply(X = node_estimates, FUN = function(x) x[tree$edge[tree$edge[, 2] == needle, 1]])), 1] + matrix(data = rep(x = node_values[needle, ], times = length(x = node_estimates)), nrow = length(x = node_estimates), byrow = TRUE)
    
    
    
    
    node_estimates <- do.call(what = rbind, args = lapply(X = apply(X = cbind(1:nrow(x = permuted_new_lengths), permuted_new_lengths), MARGIN = 1, FUN = list), FUN = function(x) {x <- unlist(x = x); column_number <- x[1]; x <- x[-1]; x == min(x = x); possible_states <- colnames(x = stepmatrix)[x == min(x = x)]; n_states <- length(x = possible_states); new_node_estimates <- lapply(X = as.list(x = possible_states), FUN = function(y) {new_node_estimates <- node_estimates[[column_number]]; new_node_estimates[needle] <- y; new_node_estimates}); do.call(what = rbind, args = new_node_estimates)}))
    
    
    node_estimates <- lapply(X = as.list(x = 1:nrow(x = node_estimates)), FUN = function(x) node_estimates[x, ])
    
    
    # SOMEHOW NEED ABOVE NOT TO BE NESTED LISTS!!!!!
 
  }


#any(x = unlist(x = lapply(X = node_estimates, FUN = is.na)))
  
  node_estimates <- do.call(what = cbind, args = node_estimates)
  

par(mfrow = c(1, ncol(node_estimates)))
par(mfrow = c(4, 6))
for(i in 1:ncol(node_estimates)) {
  plot(tree, show.tip.label = FALSE)
  tiplabels(node_estimates[1:tip_count], cex = 2)
  nodelabels(node_estimates[(1 + tip_count):node_count, i], cex = 2)
}


  
  
  
  # TEST BIT BELOW: TEST FAILS!
  
  # Set up list to store estimate(s) for every node:
  #node_estimates <- rep(x = NA, times = node_count)
  
  # Populate tip values immediately:
  #node_estimates[1:tip_count] <- tip_states[tree$tip.label]

  # Get most parsimonious root state(s):
  #root_states <- colnames(node_values)[node_values[tip_count + 1, ] == min(x = node_values[tip_count + 1, ])]
  
  # Set a single root state (picking arbitrarily if multiple options) and store as first node estimate:
  #node_estimates[tip_count + 1] <- ifelse(length(x = root_states) > 0, sample(x = root_states, size = 1), root_states)
  
  # INSERT CONDITIONAL SOMEWHERE TO DO QUICKLY IF NO AMBIGUITIES!
  
  # Second pass (traverse tree from root to tips):
  #for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
    
    # Find ancestor of current node:
    #ancestor <- tree$edge[tree$edge[, 2] == needle, 1]
    
    # Find ancestor state:
    #ancestor_state <- node_estimates[ancestor]
    
    # Find possible states at current node:
    #possible_states <- colnames(node_values)[node_values[needle, ] == min(x = node_values[needle, ])]
    
    # If there are multiple possible states:
    #if(length(x = possible_states) > 1) {
      
      # Find intersecting state(s) of current node and its' ancestor:
      #intersecting_states <- intersect(x = ancestor_state, y = possible_states)
      
      # If a single intersecting state, then store that:
      #if(length(x = intersecting_states) == 1) {
        
        # Simply store that value as the estimate:
        #node_estimates[needle] <- intersecting_states
        
      # If not singe intersecting state:
      #} else {
        
        # DEFO BREAKS HERE!
        
        #ancestor_distances <- stepmatrix[ancestor_state, possible_states]
        
        # DELTRAN SO PICK CLOSEST VALUE:
        #node_estimates[needle] <- names(ancestor_distances)[ancestor_distances == min(x = ancestor_distances)]
        
        # ACCTRAN SO PICK LARGEST VALUE
        #node_estimates[needle] <- names(ancestor_distances)[ancestor_distances == max(x = ancestor_distances)]
        
        #}
      
    # If only a single possible state:
    #} else {
      
      # Assign that state as the estimate:
      #node_estimates[needle] <- possible_states
      
      #}
    
    #}
    
    
    
  


  # OUTPUT:
  # - Ancestral states
  # - Ambiguities?
  # - Stepmatrix used
  # - Tree length
  # - How root value was chosen (arbitraty forced, unambiguous?)
  # - Algorithm (parsimony, ML, etc.)
  # - CI?
  # - RI?
  # - N reversals
  # - N parallelisms
  # - Plot option (one panel for each possible MPR?)
  # - MAKE THIS A CLASS TOO!
  
}

tree <- ape::read.tree(text = "((C,(A,B)),(D,E));")
tip_states <- c(A = 1, B = 2, C = 2, D = 1, E = 0)

# Ordered:
plot(tree, show.tip.label = FALSE); tiplabels(tip_states[tree$tip.label]); nodelabels(unlist(lapply(apply(parsimony_ancestors(tree, tip_states = c(A = 1, B = 2, C = 2, D = 1, E = 0), stepmatrix = "ordered")[(ape::Ntip(tree) + 1):(ape::Ntip(tree) + tree$Nnode), , drop = FALSE], 1, list), function(x) {x <- unlist(x); paste(names(x)[x == min(x)], collapse = ",")})))

# Unordered:
plot(tree, show.tip.label = FALSE); tiplabels(tip_states[tree$tip.label]); nodelabels(unlist(lapply(apply(parsimony_ancestors(tree, tip_states = c(A = 1, B = 2, C = 2, D = 1, E = 0), stepmatrix = "unordered")[(ape::Ntip(tree) + 1):(ape::Ntip(tree) + tree$Nnode), , drop = FALSE], 1, list), function(x) {x <- unlist(x); paste(names(x)[x == min(x)], collapse = ",")})))

# Penalise gains:
plot(tree, show.tip.label = FALSE); tiplabels(tip_states[tree$tip.label]); nodelabels(unlist(lapply(apply(parsimony_ancestors(tree, tip_states = c(A = 1, B = 2, C = 2, D = 1, E = 0), stepmatrix = matrix(data = c(0,1,2,999,0,1,999,999,0), nrow = 3, ncol = 3, dimnames = list(0:2, 0:2)))[(ape::Ntip(tree) + 1):(ape::Ntip(tree) + tree$Nnode), , drop = FALSE], 1, list), function(x) {x <- unlist(x); paste(names(x)[x == min(x)], collapse = ",")})))

# Penalise losses:
plot(tree, show.tip.label = FALSE); tiplabels(tip_states[tree$tip.label]); nodelabels(unlist(lapply(apply(parsimony_ancestors(tree, tip_states = c(A = 1, B = 2, C = 2, D = 1, E = 0), stepmatrix = matrix(data = c(0,999,999,1,0,999,2,1,0), nrow = 3, ncol = 3, dimnames = list(0:2, 0:2)))[(ape::Ntip(tree) + 1):(ape::Ntip(tree) + tree$Nnode), , drop = FALSE], 1, list), function(x) {x <- unlist(x); paste(names(x)[x == min(x)], collapse = ",")})))

# Ordered with polytomies:
tree <- ape::read.tree(text = "((C,A,B),(D,E));") # Two non-root internal
plot(tree, show.tip.label = FALSE); tiplabels(tip_states[tree$tip.label]); nodelabels(unlist(lapply(apply(parsimony_ancestors(tree, tip_states = c(A = 1, B = 2, C = 2, D = 1, E = 0), stepmatrix = "ordered")[(ape::Ntip(tree) + 1):(ape::Ntip(tree) + tree$Nnode), , drop = FALSE], 1, list), function(x) {x <- unlist(x); paste(names(x)[x == min(x)], collapse = ",")})))
tree <- ape::read.tree(text = "((C,A,B),D,E);") # One non-root internal
plot(tree, show.tip.label = FALSE); tiplabels(tip_states[tree$tip.label]); nodelabels(unlist(lapply(apply(parsimony_ancestors(tree, tip_states = c(A = 1, B = 2, C = 2, D = 1, E = 0), stepmatrix = "ordered")[(ape::Ntip(tree) + 1):(ape::Ntip(tree) + tree$Nnode), , drop = FALSE], 1, list), function(x) {x <- unlist(x); paste(names(x)[x == min(x)], collapse = ",")})))
tree <- ape::read.tree(text = "(C,A,B,D,E);") # Star tree
plot(tree, show.tip.label = FALSE); tiplabels(tip_states[tree$tip.label]); nodelabels(unlist(lapply(apply(parsimony_ancestors(tree, tip_states = c(A = 1, B = 2, C = 2, D = 1, E = 0), stepmatrix = "ordered")[(ape::Ntip(tree) + 1):(ape::Ntip(tree) + tree$Nnode), , drop = FALSE], 1, list), function(x) {x <- unlist(x); paste(names(x)[x == min(x)], collapse = ",")})))

