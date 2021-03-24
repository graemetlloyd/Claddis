#' Reconstruct ancestral states using parsimony
#'
#' @description
#'
#' Given a tree, discrete states for each tip, and either a specific stepmatrix or character type, returns the most parsimonious ancestral reconstruction(s).
#'
#' @param tree A tree (phylo object).
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}. These should be discrete with labels matching the tip labels of \code{tree}.
#' @param inapplicables_as_missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default).
#'
#' @details
#'
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
#' \emph{Stepmatrices}
#'
#' Although rare as explicit statements (see Hooker 2014 for a rare example), any character type (e.g., ordered, unordered, Dollo) can be expressed as a step matrix. These are always square, with rows and columns representing the same set of states that will typically reflect (but are not restricted to) the range of values observed in the data. Individual values represent the cost (in steps) of each transition, with the row value being the "from" state and the column value being the "to" state. By convention the cost of going from any character to itself is zero and hence so is the diagonal of the matrix. An example stepmatrix for a three state unordered character would thus look like this:
#'
#' \preformatted{   -------------
#'    | 0 | 1 | 2 |
#' ----------------
#'  0 | 0 | 1 | 1 |
#' ----------------
#'  1 | 1 | 0 | 1 |
#' ----------------
#'  2 | 1 | 1 | 0 |
#' ----------------}
#'
#' Hence going from state 2 (third row) to state 1 (second column) reveals a cost of 1. Indeed, as this is an unordered character \emph{any} off-diagonal value is the same (here 1). By contrast an ordered matrix might look like this:
#'
#' \preformatted{   -------------
#'    | 0 | 1 | 2 |
#' ----------------
#'  0 | 0 | 1 | 2 |
#' ----------------
#'  1 | 1 | 0 | 1 |
#' ----------------
#'  2 | 2 | 1 | 0 |
#' ----------------}
#'
#' Now going from state 0 (first row) to state 2 (third column) costs 2 steps as by implication this transition must pass through the intermediate state 1.
#'
#' So far these examples are symmetric - i.e., you can imagine the diagonal as a line of reflection, or alternatively that going from state X to state Y \emph{always} costs the same as going from state Y to state X. However, asymmetric stepmatrices are also possible and there are multiple such character=types that may be relevant here, e.g., Dollo, Camin-Sokal, and stratigraphy.
#'

#' Dollo: (bird teeth example Field and Brocklehurst) M is large but not infinite as an acquisition is still expected, but should still outweigh any second acqusition such that when optimised multiple losses are allowed.
#'
#' \preformatted{   -------------------
#'    |  0  |  1  |  2  |
#' ----------------------
#'  0 |  0  |  1  |  2  |
#' ----------------------
#'  1 |  M  |  0  |  1  |
#' ----------------------
#'  2 |  2M |  M  |  0  |
#' ----------------------}

#' Camin-Sokal (irreversible): Less obvious what an example may be, perhaps organelles in eukaryotes?
#'
#' \preformatted{   -------------------
#'    |  0  |  1  |  2  |
#' ----------------------
#'  0 |  0  | Inf | Inf |
#' ----------------------
#'  1 |  1  |  0  | Inf |
#' ----------------------
#'  2 |  2  |  1  |  0  |
#' ----------------------}

#' Stratigraphy (irreversible, but costs reflect time differences, i.e., in millions of years).
#'
#' 0 - Santonian (85.8 - 83.5 Ma)
#' 1 - Campanian (83.5 - 70.6 Ma)
#' 2 - Maastrichtian (70.6 - 65.5 Ma)
#'
#' Midpoints: 0 (84.65 Ma), 1 (77.05 Ma), and 2 (68.05 Ma)
#'
#' \preformatted{    --------------------
#'     |   0  |  1  |  2  |
#' ------------------------
#' | 0 |   0  | Inf | Inf |
#' ------------------------
#' | 1 |  7.6 |  0  | Inf |
#' ------------------------
#' | 2 | 16.6 |  9  |  0  |
#' ------------------------}




#' [MORE HERE]
#'
#' These are (in a very general way) related to the Q-matrices used in likelihood based approaches, but note that their use has multiple key differences and a Q-matrix cannot be directly interpretetd as a stepmatrix or vice versa.
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
#' Polymorphisms represent perhaps the most complex problem for phylogenetic inference and ancestral state estimation (reconstruction), at least in part due to the varied nature their meaning can have (Swofford and Maddison 1992). For example, they may represent true variance within a species or composite variance in some higher taxon represented by a single OTU (something to be avoided if at all possible). Similarly, as different cladistic matrix formats do not necessarily allow the distinction between uncertainties and polymorphisms they can sometimes represent the former even though their encoding suggest the latter. (For example, TNT (Goloboff et al. 2008; Goloboff and Catalano 2016) conflates the two with a single representation within square brackets, [].)
#'
#' Misrepresenting true polymorphisms as uncertainty leads to incorrect estimates of the amount of evolution (tree length under parsimony; Nixon and Davis 1991) and disallows, potentially inapproriately, polymorphic values as ancestral state estimates (reconstructions).

#' Really need explicit options before talking about this too much further.
#' Maddison and Maddison (1987) suggested stepmatrix solutionss that allow for polymorphic ancestors.
#' fitpolyMk in phytools [LINK?]

#' An unordered stepmatrix using the Maddison and Maddison (1987) approach (NB: modified such that the single state transitions are scaled to one):
#'
#' \preformatted{        ---------------------------------------------
#'         |  0  |  1  |  2  | 0&1 | 0&2 | 1&2 | 0&1&2 |
#' -----------------------------------------------------
#' |   0   |  0  |  1  |  1  | 0.5 | 0.5 | 1.5 |   1   |
#' -----------------------------------------------------
#' |   1   |  1  |  0  |  1  | 0.5 | 1.5 | 0.5 |   1   |
#' -----------------------------------------------------
#' |   2   |  1  |  1  |  0  | 1.5 | 0.5 | 0.5 |   1   |
#' -----------------------------------------------------
#' |  0&1  | 0.5 | 0.5 | 1.5 |  0  |  1  |  1  |  0.5  |
#' -----------------------------------------------------
#' |  0&2  | 0.5 | 1.5 | 0.5 |  1  |  0  |  1  |  0.5  |
#' -----------------------------------------------------
#' |  1&2  | 1.5 | 0.5 | 0.5 |  1  |  1  |  0  |  0.5  |
#' -----------------------------------------------------
#' | 0&1&2 |  1  |  1  |  1  | 0.5 | 0.5 | 0.5 |   0   |
#' -----------------------------------------------------}
#'
#' These can also be represented as Manhattan distance between the vertices of a hypercube plotted in a space with axes representing the presence of each individual state (0, 1, or 2). Hence the value 0&1 would have the coordinates 1,1,0. For this three state example this will be a cube with every vertice except the origin (0,0,0) occupied by a possible state. (NB: This is also the reason for the N^2 - 1 growth of stepmatrix size given below.)
#'
#' However, this approach may not be ideal as it can lead to scenarios where a polymorphism may be seen as non-optimal. For example, the tree:
#'
#' \preformatted{0   1   2
#'  \  |  /
#'   \ | /
#'    \|/
#'     X}
#'
#' If we consider the possibilities for ancestor X being in each state they are:
#'
#' \preformatted{---------------------
#' | Ancestral | Total |
#' |   state   | steps |
#' ---------------------
#' |     0     |   2   |
#' |     1     |   2   |
#' |     2     |   2   |
#' |    0&1    |  2.5  |
#' |    0&2    |  2.5  |
#' |    1&2    |  2.5  |
#' |   0&1&2   |   3   |
#' ---------------------}
#'
#' This seems like a scenario where 0&1&2 ought to be the optimal solution, or at least equally optimal, and yet it comes out as the maximally suboptimal solution. To allow such solutions to be considered optimal this function takes the basic concept of a coordinate space with axes representing each state, but allows ploting of values on either the surface of a hypersphere or a simplex, i.e., the dimensionally-scalable versions of a circle and an equilateral traingle, respectively. (Shapes that place polymorphic values closer to the origin and hence make them more parsimonious.) Similarly, the user can use either the Euclidean [OPTION], Manhattan [OPTION], or Great Circle [OPTION] distance between points. Here it is assumed these will match the shape used, specifically that the Hypercube uses Manhattan distances, a Hypersphere Great Circle distances, and a Simplex Euclidean distances, but no restrictions prevent other combinations being applied. (Caveat emptor, always.)
#'
#' Things become more complex when we consider ordered characters. For "intermediate" polymorphisms they can be simple, i.e.:
#'
#' \preformatted{/---\  0.5  /---\  0.5  /---\  0.5  /---\  0.5  /---\
#' | 0 |<----->|0&1|<----->| 1 |<----->|1&2|<----->| 2 |
#' \---/       \---/       \---/       \---/       \---/}
#'
#' Values between states are steps, leading to the following step matrix:
#'
#' \preformatted{        -------------------------------
#'         |  0  |  1  |  2  | 0&1 | 1&2 |
#' ---------------------------------------
#' |   0   |  0  |  1  |  2  | 0.5 | 1.5 |
#' ---------------------------------------
#' |   1   |  1  |  0  |  1  | 0.5 | 0.5 |
#' ---------------------------------------
#' |   2   |  2  |  1  |  0  | 1.5 | 0.5 |
#' ---------------------------------------
#' |  0&1  | 0.5 | 0.5 | 1.5 |  0  |  1  |
#' ---------------------------------------
#' |  1&2  | 1.5 | 0.5 | 0.5 |  1  |  0  |
#' ---------------------------------------}
#'
#' However, it is less obvious what to do with the "non-intermediate" polymorphisms, 0&2 and 0&1&2. The only suggested solution of which I am aware is that of Maddison and Maddison (2000) where changes (0.5 steps here) represent adding OR losing a state. Thus going from 0 to 1 involves 0.5 steps to gain 1 and then 0.5 steps to lose 0, and going from 0 to 0&1&2 involves (first) adding 1 (as this is already adjacent, i.e., accessible, to zero) then another 0.5 steps to add state 2 (now accessible via the state 1). Hence the extended stepmatrix would be:
#'
#' \preformatted{        ---------------------------------------------
#'         |  0  |  1  |  2  | 0&1 | 0&2 | 1&2 | 0&1&2 |
#' -----------------------------------------------------
#' |   0   |  0  |  1  |  2  | 0.5 | 1.5 | 1.5 |   1   |
#' -----------------------------------------------------
#' |   1   |  1  |  0  |  1  | 0.5 | 1.5 | 0.5 |   1   |
#' -----------------------------------------------------
#' |   2   |  2  |  1  |  0  | 1.5 | 1.5 | 0.5 |   1   |
#' -----------------------------------------------------
#' |  0&1  | 0.5 | 0.5 | 1.5 |  0  |  1  |  1  |  0.5  |
#' -----------------------------------------------------
#' |  0&2  | 1.5 | 1.5 | 1.5 |  1  |  0  |  1  |  0.5  |
#' -----------------------------------------------------
#' |  1&2  | 1.5 | 0.5 | 0.5 |  1  |  1  |  0  |  0.5  |
#' -----------------------------------------------------
#' | 0&1&2 |  1  |  1  |  1  | 0.5 | 0.5 | 0.5 |   0   |
#' -----------------------------------------------------}
#'
#' NOW DEAL WITH DOLLO AND CAMIN-SOKAL AND STRATIGRAPHY VERSIONS (THIS POLYMORPHISM ISSUE IS WHERE STRATIGRAPHY DIFFERS FROM IRRERVERSIBLE CHARACTERS).





#'
#' Forming stepmatrices in this way has a clear size limit. Specifically, it will require 2^N - 1 terms where N is the number of single states. Thus the size of stepmatrix required grows very quickly:
#'
#' \preformatted{------------------------------
#' | N states | Stepmatrix size |
#' ------------------------------
#' |     2    |      3 x 3      |
#' |     3    |      7 x 7      |
#' |     4    |     15 x 15     |
#' |     5    |     31 x 31     |
#' |     6    |     63 x 63     |
#' |     7    |    127 x 127    |
#' |     8    |    255 x 255    |
#' |     9    |    511 x 511    |
#' |    10    |   1023 x 1023   |
#' |    11    |   2047 x 2047   |
#' |    12    |   4095 x 4095   |
#' |    13    |   8191 x 8191   |
#' |    14    |  16383 x 16383  |
#' ------------------------------}
#'
#' Because of this the function will become extremely slow for these higher values and indeed does not allow more than fourteen states.
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
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics/ \emph{Cladistics}, \bold{32}. 221-238
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' Hooker, J. J., 2014. New postcranial bones of the extinct mammalian family Nyctitheriidae (Paleogene, UK): primitive euarchontans with scansorial locomotion. \emph{Palaeontologia Electronica}, \bold{17.3.47A}, 1-82.
#'
#' Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. \emph{Palaeontology}, \bold{61}, 637-645.
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
#' # Calculate maximum parsimony tree length for Gauthier 1986:
#' cladistic_matrix <- Claddis::gauthier_1986
#' tree <- ape::read.tree(
#'   text = "(Outgroup,(Ornithischia,(Sauropodomorpha,(Ceratosauria,Procompsognathus,
#'   Liliensternus,(Carnosauria,(Ornithmimidae,Saurornitholestes,Hulsanpes,(Coelurus,
#'   Elmisauridae,(Compsognathus,(Ornitholestes,Microvenator,Caenagnathidae,
#'   (Deinonychosauria,Avialae))))))))));"
#' )
#' inapplicables_as_missing = FALSE
#' parsimony_function(
#'   tree = tree,
#'   cladistic_matrix,
#'   inapplicables_as_missing = FALSE
#' )
#'
#' @export parsimony_function
parsimony_function <- function(tree, cladistic_matrix, inapplicables_as_missing = FALSE) {
  
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

  reconstruct_ancestral_states <- function(tree, tip_states, stepmatrix, weight = 1, inapplicables_as_missing = FALSE) {
    
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
    
    # Return compiled output:
    list(length = tree_length, most_parsimonious_reconstructions = node_estimates, input_tree = input_tree)
    
  }
  
  
  # TO ADD INTO FUNCTION:
  estimate_all_nodes <- FALSE
  estimate_tip_values <- FALSE # Should not happen for inapplicables when they are treated as inapplicables
  
  # Combine all blocks into a single input matrix:
  single_input_matrix <- list(
    matrix = do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$matrix)),
    ordering = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$ordering))),
    weights = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$character_weights))),
    minima = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$minimum_values))),
    maxima = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$maximum_values)))
  )
  
  # Above needs to deal with actual stepmatrices:
  
  full_output <- lapply(X = as.list(1:ncol(x = single_input_matrix$matrix)), FUN = function(x) reconstruct_ancestral_states(tree = tree, tip_states = single_input_matrix$matrix[, x], stepmatrix = "unordered", weight = single_input_matrix$weights[x], inapplicables_as_missing = inapplicables_as_missing))
  
  
  tree_length <- do.call(what = sum, args = lapply(X = full_output, FUN = function(x) x$length))
  
  


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


