#' Reconstruct ancestral states using parsimony
#'
#' @description
#'
#' Given a tree, discrete states for each tip, and either a specific stepmatrix or character type, returns the most parsimonious ancestral reconstruction(s).
#'
#' @param tree A tree (phylo object).
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}. These should be discrete with labels matching the tip labels of \code{tree}.
#' @param inapplicables_as_missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default).
#' @param include_polymorphisms Argument passed to \link{make_stepmatrix}.
#' @param polymorphism_shape Argument passed to \link{make_stepmatrix}.
#' @param polymorphism_distance Argument passed to \link{make_stepmatrix}.
#' @param state_ages Argument passed to \link{make_stepmatrix}.
#' @param dollo_penalty Argument passed to \link{make_stepmatrix}.
#'
#' @details
#'


# MAKE A TREE-LENGTH FUNCTION (THAT HAPPENS TO EXPORT NODE VALUES)
# THEN MAKE SEPARATE ANCESTRAL STATE RECONSTRUCTION FUNCTION THAT CALLS THIS OTHER FUNCTION



#' Can be used to infer ancestral states or get tree lengths (could theoretically buol a tree inference procedure around this by proposing different topologies, but ikely to be much slower than compiled software). Also a first step to generating character maps that can be used for plotting changes on a tree, performing rate analysis (see \link{test_rates}), or generating phylomorphospaces (not implemented yet).
#'
#' \bold{Input}
#'
#' To use this function data must first be in the correct format.
#'
#' - Tree
#' - Cladistic matrix (currently tip_states)
#'
#' \bold{Options}
#'
#' Come in two forms, those implicit (as decided by values in cladisticmatrix) and explicit to this function.
#'
#' \emph{Implicit options}
#'
#' - Stepmatrix
#' - Weight (only affects tree length) - defaults to one as this is typical value used. This can disappear once cladistic matrix is input.
#'
#' \emph{Explicit options}
#'
#' - inapplicables_as_missing = FALSE, what to do with inapplicable characters. Currently treat as missing is only real option, but can be differet state that is also inferred as an ancestor.
#' - More for uncertianties and polymorphisms (see below).
#' - Add option to predict missing tip values (but definitely not recommended! - see Lloyd 2018)
#' - Collpasing to a single output with ACCTRAN etc.?
#'
#' \bold{Implementation}
#'
#' Algorithms for reconstructing ancestral states under parsimony can be specific (e.g., for ordered characters only) or general (allowing for any character type, including highly specific custom-formats). The advantage of the former is usually that they are faster as multiple assumptions are already built in, but the latter offer more flexibility and represent the version applied here.
#'
#' More specifically this function is based on the approach described by Swofford and Maddison (1992) that assumes a rooted tree and uses a stepmatrix to define the costs (in number of steps) of specific transitions, i.e., from state X to state Y. The rooting is key here as it sets the directionality of evolution (i.e., time) and hence allows the use of asymmetric characters (where the cost of going from state X to state Y is different than going from stete Y to state X) as well as other direction-based optimisations (e.g., ACCTRAN and DELTRAN).
#'
#'
#' \bold{Special character types}
#' \emph{Missing values}
#' \emph{Inapplicables}
#' \emph{Uncertainties}
#'
#' In Claddis, uncertainties are coded by separating states with the slash character, e.g., 0/1. The Swofford and Maddison (1992) algorithm - and the implementation used here - treats uncertainties (automatically) by allowing any of the possible states at the tip (here 0 or 1, but not, say, 2 or 3). In other words, it will opt for solutions that effectively prefer whichever state is easiest (fewest steps) to reach from a proposed ancestral value. This could therefore be used to predict what the "true" value may be, at least under the assumption of parsimony. However, such phylogenetic "prediction" [OPTION?] is generally advised against (Lloyd 2018).
#'
#' \emph{Polymorphisms}
#'
#' TEXT?
#'
#' \bold{Polytomies}
#'
#' Another advantage of the Swofford and Maddison (1992) approach is that it does not require fully-bifurcating topologies, hence these are "allowed" as input. However, it is important to note that the implementation here will treat any such polytomies (also known as multifurcations) as "hard". I.e., that they genuinely represent one species splitting into three or more descendant species. As always, the Ian Malcolm caveat applies - "your scientists were so preoccupied with whether or not they sould, they didn’t stop to think if they should."


#'
#'
#' # Option to fix root state (polarize) somehow
#' # From is row, to is column (directional hence requires/assumes rooted trees)
#' # Record in manual that weighted parsimony is folly, but could use to chose between MPRs (maybe?)

#' \bold{Polymorphic characters}
#'

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
#' Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. \emph{Palaeontology}, \bold{61}, 637-645.
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
#' # Calculate maximum parsimony tree length for Gauthier 1986:
#' cladistic_matrix <- Claddis::gauthier_1986
#' tree <- ape::read.tree(
#'   text = "(Outgroup,(Ornithischia,(Sauropodomorpha,(Ceratosauria,Procompsognathus,
#'   Liliensternus,(Carnosauria,(Ornithmimidae,Saurornitholestes,Hulsanpes,(Coelurus,
#'   Elmisauridae,(Compsognathus,(Ornitholestes,Microvenator,Caenagnathidae,
#'   (Deinonychosauria,Avialae))))))))));"
#' )
#' inapplicables_as_missing = FALSE
#' calculate_tree_length(
#'   tree = tree,
#'   cladistic_matrix,
#'   inapplicables_as_missing = FALSE
#' )
#'
#' @export calculate_tree_length
calculate_tree_length <- function(tree, cladistic_matrix, inapplicables_as_missing = FALSE, include_polymorphisms, polymorphism_shape, polymorphism_distance, state_ages, dollo_penalty) {
  
  # WHOLE THING SHOULD BE REFACTORED TO DEAL WITH ENTIRE MATRIX (THE OBVIOUS END USE) PLUS BRINGS INTO LINE WITH OTHER FUNCTIONS AND TRANSFERS MANY OPTIONS TO MATRIX STRUCTURE (ORDERING AND WEIGHTS PRIMARILY)
  # ALSO ALLOWS MOVING MANY CHECKS TO A SINGLE PASS OF MATRIX AND HENCE FASTER
  # ALLOW FOR MULTIPLE TREES?
  
  # TO DO
  # - Conditional for invariant character plus autapomorphic ones (fast solution that escapes estimates of other values)
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
  
  
  
  # Combine all blocks into a single input matrix:
  single_input_matrix <- list(
    matrix = do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$matrix)),
    ordering = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$ordering))),
    weights = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$character_weights))),
    minima = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$minimum_values))),
    maxima = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$maximum_values)))
  )
  
  # Build character list:
  character_list <- lapply(
    X = as.list(x = 1:length(x = single_input_matrix$ordering)),
    FUN = function(x) list(
      tip_states = single_input_matrix$matrix[, x],
      stepmatrix = make_stepmatrix(
        min_state = single_input_matrix$minima[x],
        max_state = single_input_matrix$maxima[x],
        character_type = single_input_matrix$ordering[x],
        include_polymorphisms = include_polymorphisms,
        polymorphism_shape = polymorphism_shape,
        polymorphism_distance = polymorphism_distance,
        state_ages = state_ages,
        dollo_penalty = dollo_penalty
      ),
      weight = single_input_matrix$weights[x],
      inapplicables_as_missing = inapplicables_as_missing
    )
  )
  
  
  
  
  
  tree <- ape::read.tree(text = "((A,B),(C,(D,E)));")
  tip_states <- c(A = 0, B = 1, C = 0, D = 2, E = 1)
  stepmatrix <- make_stepmatrix(min_state = 0, max_state = 2, character_type = "ordered", include_polymorphisms = FALSE, polymorphism_shape = "hypersphere", polymorphism_distance = "great_circle")
  weight <- 1
  inapplicables_as_missing = FALSE




  reconstruct_ancestral_states <- function(tree, tip_states, stepmatrix, weight, inapplicables_as_missing = FALSE) {
    
    # Reorder tip states (1 to N) and store an unmodified version:
    pristine_tip_states <- tip_states <- tip_states[tree$tip.label]
    
    # Store pristine input tree:
    input_tree <- tree
    
    # If there are no branch durations set these as all one:
    if (is.null(x = tree$edge.length[1])) tree$edge.length <- rep(x = 1, length.out = nrow(x = tree$edge))
    
    # If treating inapplicables as missing then replace any inapplicable tip state with NA:
    if (inapplicables_as_missing) tip_states[tip_states == ""] <- NA
    
    # Reformat tip states as character vector:
    tip_states <- as.character(x = tip_states)
    
    # Re-add names to tip states:
    names(tip_states) <- tree$tip.label
    
    # Establish number of tips:
    tip_count <- ape::Ntip(phy = tree)
    
    # Establish total number of nodes (terminal and internal):
    node_count <- tip_count + tree$Nnode
    
    # Initialise node_values matrix:
    node_values <- matrix(data = 0, nrow = node_count, ncol = stepmatrix$size, dimnames = list(c(), colnames(stepmatrix$stepmatrix)))
    
    # Begin by inserting Inf as default tip value:
    node_values[1:tip_count, ] <- Inf
    
    # Populate tip values (zeroes for assigned state(s):
    for(i in 1:tip_count) node_values[i, match(x = strsplit(x = tip_states[i], split = "/")[[1]], table = colnames(x = node_values))] <- 0
    
    # If any tips are missing then set all states for those as zero:
    if (any(x = is.na(x = tip_states))) node_values[which(x = is.na(x = tip_states)), ] <- 0
    
    ### FORCE ROOT STATES IF ASKED TO:
    
    # First pass (traverse tree from tips to root):
    for(needle in (tip_count + node_count - tip_count):(tip_count + 1)) {
      
      # Find decsendants of current node:
      descendants <- tree$edge[tree$edge[, 1] == needle, 2]
      
      # Calculate and store new node values:
      node_values[needle, ] <- apply(
        X = rbind(
          node_values[needle, ],
          unlist(x = lapply(
            X = as.list(x = colnames(x = node_values)),
            FUN = function(fromstate) {
              sum(x = unlist(x = lapply(
                X = as.list(x = descendants),
                FUN = function(descendant) min(x = node_values[descendant, ] + stepmatrix$stepmatrix[fromstate, ]))))
            })
          )
        ),
        MARGIN = 2,
        FUN = max
      )
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
    if (any(x = is.na(x = node_estimates[(tip_count + 1):node_count]))) {
      
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
          updated_tree_lengths <- stepmatrix$stepmatrix[x[ancestor_node], ] + node_values[needle, ]
          
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
    
    # Return compiled output:
    list(length = tree_length, most_parsimonious_reconstructions = node_estimates, input_tree = input_tree)
    
  }
  
  
  
  ancestral_values <- lapply(
    X = character_list,
    function(x) reconstruct_ancestral_states(
      tree = tree,
      tip_states = x$tip_states,
      stepmatrix = x$stepmatrix,
      weight = x$weight,
      inapplicables_as_missing = x$inapplicables_as_missing
    )
  )
  
  
  
  
  
  # TO ADD INTO FUNCTION:
  estimate_all_nodes <- FALSE
  estimate_tip_values <- FALSE # Should not happen for inapplicables when they are treated as inapplicables
  
  
  
  
  
  
  




  
  # KEY THING HERE IS ALLOWS POLYMORPHISMS AT INTERNAL NODES!
  # KEY QUESTION IS ARE THERE EVER CASES WHERE //ALL// POSSIBLE MPRS INCLUDE A POLYMORPHIC ANCESTOR?
  
  
  
  # NEED TO COUNT DOLLO CHARACTERS DIFFERENTLY FOR PARSIMONY OR WILL OVERWHELM BROADER CHARACTER TRENDS.
  
  # ADD SWOFFORD REF TO DESCRIPTION FILE (PROLLY NO DOI!)
  
  
  
  
  # Output:
  list(tree_length = tree_length)
  
  
  

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
# TEST WEIGHTING OF STRATOCLADISTICS BY USING A STRATIGRAPHIC CHARACTER AND WEIGHTING IT MULTIPLE WAYS (AS A SLIDER) AND SHOW HOW IT EFFECTS PARSIMONY VS STRAT CONGRUENCE TREE LANDSCAPES

# ADD POLLY 1997 TO PACKAGE AS EXAMPLE OF STRATOCLADISTICS FUNCTION?
# Fisher argues (correctly) that treelength is really two components, minimum change based on the matrix plus additional changes implied by the tree (homoplas), the latter of which he terms parsimony debt. Something to think about for homoplasy metrics.

# ACCTRAN
# DELTRAN
# MINF (Swofford and Maddison 1987) - I THINK NOT!
# MINSTATE *shrug emoi*
# MAXSTATE *shrug emoi*
# WEIGHTED BY BRANCH DURATION?


