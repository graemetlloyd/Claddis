#' Converts a character state tree to an adjacency matrix
#'
#' @description
#'
#' Takes a character state tree as input and returns the corresponding adjacency matrix.
#'
#' @param state_tree A text string describing a character state tree.
#'
#' @details
#'
#' There are multiple ways to define discrete character states for use in phylogenetic analysis or ancestral state reconstruction in terms of transitions between states. The two most common being "ordered" ((Wagner optimization; Farris, 1970) or "unordered" (Fitch optimization; Fitch, 1971). These can be sonsidered as representing simple Markov models with symmetric transition probabilities, with the only difference being whether any direct state-tostae transition is possible (unordered) or only immediately "adjacent" states can be transitioned too.
#'
#' In practice, however, there may be reason to believe in more complex transition models. For an example, see various characters in Hooker (2014). If these are still symmetric (i.e., going from state X to state Y is just as probable as going from state Y to state X) then they can be expressed using a "character state tree". For example, the model:
#'
#' \preformatted{  0-1-2
#'   |
#'   3}
#'
#' This model has four states (labelled 0-3). However, it is clearly not an ordered model:
#'
#' \preformatted{  0-1-2-3}
#'
#' Where only adjacent state transitions are possible, e.g., to get from state 0 to state 2 requires passing through state 1.
#'
#' Nor is it an unordered model:
#'
#' \preformatted{  0---1
#'   |\ /|
#'   | X |
#'   |/ \|
#'   2---3}
#'
#' Where any state-to-state transition is possible.
#'
#' In other words we cannot simply label this model as "ordered" or "unordered". It can, however, be expressed directly as character state tree:
#'
#' \preformatted{  (3,(2)1)0}
#'
#' Superficially these appear very similar to Newick strings (ways or representing bifurcating or multifurcating trees of species). However, the major difference is that with Newick strings we often only label the "tips" (terminal nodes, or vertices of degree one). Here there are no unlabelled "internal nodes", also known as "hypothetical ancestors". Instead, numbers after parentheses indicate the label of that node (vertex).
#'
#' Thus we can read the above tree as node 2 is only directly connected to node 1, 3 and 1 are directly connected to node zero.
#'
#' This is a very compact way of representing a Markov model, but unfortunately it is hard to read and difficult to do anything quantitative with. Thus this function takes such text strings and converts them to a more easily usable graph representation: the adjacency matrix.
#'
#' Adjacency matrices represent graphs by showing links between vertices as values of one in symmetric square matrices. These always have a diagonal of zero (you cannot be adjacent to yourself) and any other vertices that are not directly linked are also coded as zero.
#'
#' Note that here such matrices are also considered symmetric and hence the diagonal is also a line of reflection.
#'
#' In practice more complex relationships may be desired, including differential weighting of specific transitions (without additional intermediate states; where going from state X to state Y has a cost other than zero or one), or asymmetric models (where going from state X to state Y has a different "cost" than going from state Y to state X). These are still possible in Claddis, and phylogenetic analysis in general, but require "costmatrices" to define them. (Character state trees are insufficent.)
#'
#' Thus, costmatrices (or their probabilistic equivalent, Q-matrices) are the only generalisable form of defining \emph{any} Markov model.
#'
#' The output from this function can also be represented as a costmatrix with \link{convert_adjacency_matrix_to_costmatrix}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Farris, J. S., 1970. Methods for computing Wagner trees. \emph{Systematic Zoology}, \bold{19}, 83-92.
#'
#' Fitch, W. M., 1971. Towards defining the course of evolution: minimum change for a specified tree topology. \emph{Systematic Zoology}, \bold{20}, 406-416.
#'
#' Hooker, J. J., 2014. New postcranial bones of the extinct mammalian family Nyctitheriidae (Paleogene, UK): primitive euarchontans with scansorial locomotion. \emph{Palaeontologia Electronica}, \bold{17.3.47A}, 1-82.
#'
#' @return
#'
#' A symmetric adjacency matrix indicating which states are linked (value 1) by an edge.
#'
#' @seealso
#'
#' \link{convert_adjacency_matrix_to_costmatrix}, \link{locate_bracket_positions}
#'
#' @examples
#'
#' # Convert a simple state tree to an adjacency matrix:
#' convert_state_tree_to_adjacency_matrix(
#'  state_tree = "(3,(2)1)0"
#' )
#'
#' # Convert a more complex state tree to an adjacency matrix:
#' convert_state_tree_to_adjacency_matrix(
#'   state_tree = "(((5)4)3,(2)1)0"
#' )
#'
#' @export convert_state_tree_to_adjacency_matrix
convert_state_tree_to_adjacency_matrix <- function(state_tree) {
  
  # TO DO:
  #
  # - Make some kind of check that state_tree is formatted correctly
  # - Possible make a state tree class
  # = Remove any whitespace in input
  # - RULE: state labels must be nchar == 1 (no "AA" or similar)
  
  # Begin process of isolating state labels:
  state_labels <- strsplit(x = state_tree, split = "\\(|)|,")[[1]]
  
  # Finish isolating state labels:
  state_labels <- state_labels[nchar(state_labels) > 0]
  
  # Build initial adjacency matrix:
  adjacency_matrix <- matrix(data = 0, ncol = length(x = state_labels), nrow = length(x = state_labels), dimnames = list(sort(x = state_labels), sort(x = state_labels)))
  
  # Break apart state tree string:
  state_tree_string <- strsplit(x = state_tree, split = "")[[1]]
  
  # Generate all tree substrings:
  all_subtrees <- apply(X = locate_bracket_positions(input_string = state_tree), MARGIN = 1, FUN = function(x) paste(state_tree_string[x[1]:(x[2] + 1)], collapse = ""))
  
  # Subfunction to prune internal parentheses from a subtree:
  prune_internal_parentheses <- function(tree_string) {
    
    # Locate bracket positions:
    bracket_positions <- locate_bracket_positions(input_string = tree_string)
    
    # If there are any internal parentheses:
    if (nrow(x = bracket_positions) > 1) {
      
      # Split string:
      split_string <- strsplit(x = tree_string, split = "")[[1]]
      
      # For every internal parentheses set values to NA:
      for(i in 2:nrow(x = bracket_positions)) split_string[bracket_positions[i, 1]:bracket_positions[i, 2]] <- NA
      
      # Reform tree string with NAs excluded:
      tree_string <- paste(x = split_string[!is.na(x = split_string)], collapse = "")
      
    }
    
    # Return tree string without internal parentheses:
    tree_string
    
  }
  
  # Reformat subtrees without internal parentheses:
  all_subtrees <- lapply(X = all_subtrees, FUN = prune_internal_parentheses)
  
  # Subfunction to convert a split string into edge(s):
  convert_string_to_edges <- function(subtree_string) {
    
    # Get bracket positions:
    bracket_positions <- locate_bracket_positions(input_string = subtree_string)
    
    # Isolate link beginnings:
    link_starts <- strsplit(x = paste(x = strsplit(x = subtree_string, split = "")[[1]][(bracket_positions[1, 1] + 1):(bracket_positions[1, 2] - 1)], collapse = ""), split = ",")[[1]]
    
    # Isolate link end:
    link_end <- strsplit(x = subtree_string, split = "")[[1]][(bracket_positions[1, 2] + 1)]
    
    # Return edges:
    cbind(link_starts, link_end)
    
  }
  
  # Get matrix of all edges:
  edges_matrix <- do.call(what = rbind, args = lapply(X = all_subtrees, FUN = convert_string_to_edges))
  
  # For each edge, populate the adjacnecy matrix accordingly:
  for(i in 1:nrow(edges_matrix)) adjacency_matrix[edges_matrix[i, 1], edges_matrix[i, 2]] <- adjacency_matrix[edges_matrix[i, 2], edges_matrix[i, 1]] <- 1
  
  # Return adjacency matrix:
  adjacency_matrix
  
}
