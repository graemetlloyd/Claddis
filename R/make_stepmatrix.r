#' Make a stepmatrix for a given set of states
#'
#' @description
#'
#' Given a set of discrete states and a character type will make the approriate stepmatrix for them.
#'
#' @param min_state The minimum character state (defaults to \code{0}).
#' @param max_state The maximum character state. Must be \code{1} or greater.
#' @param character_type The type of character desired. Must be one of: \code{"ordered"}, \code{"unordered"}, \code{"dollo"}, \code{"irreversible"}, or \code{"stratigraphy"}.
#' @param include_polymorphisms Logical indicating whether or not to include polymorphic state combinations (defaults to \code{FALSE}).
#' @param polymorphism_shape The shape to use for assigning polymorphism coordinates. Must be one of: \code{"hypercube"}, \code{"hypersphere"}, or \code{"simplex"}. (Only relevant if \code{character_type = "unordered"} and \code{include_polymorphisms = TRUE}.)
#' @param polymorphism_distance The distance to use to set transformation costs between polymorphic states. Must be one of: \code{"euclidean"}, \code{"great_circle"}, or \code{"manhattan"}.
#' @param state_ages A vector of ages assigned to each state (only relevant if \code{character_type = stratigraphy}).
#' @param dollo_penalty The size of the cost penalty for the acquisition of a Dollo character (defaults to \code{999}). Note: this should always be a positive real value greater than one, and never infinity (\code{Inf}), as at least one acquisition is expected.
#'
#' @details
#'
#' Stepmatrices encode the parsimony cost (number of evolutionary "steps") for each possible state-to-state transition and can be used to estimate total lengths and make estimates (often also called reconstructions) for ancestral states at branching points on a given topology (tree). This function automates the generation of some common stepmatrix types for use elsewhere in Claddis and hence is primarily an internal function. However, as stepmatrices are fundamental to various key analyses their operation is detailed extensively below.
#'
#' \bold{Stepmatrix basics}
#'
#' Although their usage is rare as explicit statements (see Hooker 2014 for an example), any character type (e.g., ordered, unordered, Dollo) can be expressed as a step matrix. These are always square, with rows and columns representing the same set of states that will typically reflect (but are not restricted to) the range of values observed in the data. Individual values represent the cost (in steps) of each transition, with the row value being the "from" state and the column value being the "to" state. By convention, the cost of going from any character to itself is zero and hence so is the diagonal of the matrix. An example stepmatrix for a three state unordered character would thus look like this:
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
#' Hence going from state 2 (third row) to state 1 (second column) reveals a cost of 1. Indeed, as this is an unordered character, \emph{any} off-diagonal value has the same cost of 1. By contrast an ordered matrix might look like this:
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
#' Now going from state 0 (first row) to state 2 (third column) costs 2 steps, as by implication this transition must pass through the intermediate state (1).
#'
#' So far these examples are symmetric - i.e., you can imagine the diagonal as a line of reflection, or, alternatively, that going from state X to state Y \emph{always} costs the same as going from state Y to state X. However, asymmetric stepmatrices are also possible and there are multiple such character-types that may be relevant here. (Specifically, Dollo, irreversible (also known as Camin-Sokal), and stratigraphic.)
#'
#' \bold{Character types}
#'
#' This function will generate a stepmatrix for every possible transition (between \code{min_state} and \code{max_state}, inclusive) for a requested \code{character_type}, each of which is detailed further below.
#'
#' \emph{Ordered} - \code{character_type = "ordered"}
#'
#' An ordered character, also known as a Wagner character (Wagner 1961) and formalised by Farris (1970), is one where the order of states matters, and that state-to-state transitions must occur through a linear series. For example, for a three state character the states must be in order of transition, 0 <-> 1 <-> 2. Additionally, all transitions are symmetric. 0 -> 1 = 1 -> 0.
#'
#' \emph{Unordered} - \code{character_type = "unordered"}
#'
#' An unordered character, also known as a Fitch character (Fitch 1971), is one where the order of states does not matter, and that state-to-state transitions are all direct and symmetric. For example, for a five-state character the transition 0 -> 4 is direct and costs the same as 3 -> 1, or any other transition, excepting those from any state to itself.
#'
#' \emph{Dollo} - \code{character_type = "dollo"}
#'
#' Under a Dollo assumption the acquisition of a \emph{derived} character state is considered to be sufficiently complex (biologically difficult) that in all probability it only occurred once. Even if lost later in evolution, it is considered that it cannot be reacquired. An example of this is the frequent loss of teeth in dinosaurs (Brocklehurst and Field 2021). Teeth were never reacquired in this group, with tooth-like serrations forming in many birds instead. Modelling such characters with a stepmatrix is challenging (indeed, it cannot truly been done with \emph{just} a stepmatrix - see below). The way it is typically done (Kitching et al. 1998) is to form the stepmatrix like this:
#'
#' \preformatted{   -------------------
#'    |  0  |  1  |  2  |
#' ----------------------
#'  0 |  0  |  1M |  2M |
#' ----------------------
#'  1 |  1  |  0  |  1M |
#' ----------------------
#'  2 |  2  |  1  |  0  |
#' ----------------------}
#'
#' Where \code{M} is some arbitrary large number, set here using the \code{dollo_penalty} option. Note that \code{M} should still be finite as the purpose is to weight acquisitions such that they are sufficiently expensive that any most parsimonious reconstruction will favour only one, but not make them so expensive that they are not favoured at all. Furthermore, in this example there are twp derived states (1 and 2) and hence logically the weighting should be similar to an ordered character. Specifically, that there should be a single acquisition of state 1 (from state 0) and that this should precede a single acquistion of state 2 (from state 1). Most of the time, though, a Dollo character will be binary (e.g., see \link{map_dollo_changes}).
#'
#' Importantly, and as stated above, a Dollo stepmatrix, unlike the ordered and unordered stepmatrices, cannot be used without further analytical restrictions. Specifically, because it assumes an asymmetric acquisition of \emph{derived} states the root or "primitive" value must be 0. (It would not be logical to apply such a steoatrix where the root state is 1 or 2.) Thus use of a Dollo character requires the additional assumption that the primitive (root) state is known \emph{a priori}. As always, it is up to the user to know that this is a valid assumption for their data, and to set a value for the \code{dollo_penalty} accordingly.
#'
#' \emph{Irreversible (Camin-Sokal)} - \code{character_type = "irreversible"}
#'
#' Sometimes confused with a Dollo character, irreversible or Camin-Sokal (Camin and Sokal 1965) characters, do not allow for any reversals to a prior state once the derived state has been acquired. An irreversible stematrix might look like this:
#'
#' \preformatted{   -------------------
#'    |  0  |  1  |  2  |
#' ----------------------
#'  0 |  0  |  1  |  2  |
#' ----------------------
#'  1 | Inf |  0  |  1  |
#' ----------------------
#'  2 | Inf | Inf |  0  |
#' ----------------------}
#'
#' The upper triangle (representing gains) is the same as for an ordered character, but the lower triangle (representing losses, or reversals) will always be comprised of infinite costs (\code{Inf}), precluding their possibility in a most parsimonious reconstruction. Thus, in practice, any number of gains is allowed, but losses never are.
#'
#' Like a Dollo character, this approach assumes the direction of evolution is known \emph{a priori} and that 0 is always the primitive, or root, state.
#'
#' \emph{Stratigraphic} - \code{character_type = "stratigraphy"}
#'
#' Stratigraphic characters, where states represent geologic time units in which an OTU is present, are also logically irreversible - as a younger unit cannot exist before (be "ancestral" to) an older one. Thus here, too, the lower triangle of the stepmatrix is always populated by infinities to preclude their possibility. However, this time the upper triangle requires additional information to populate its' values, i.e., the gap between states in units of time (typically millions of years).
#'
#' This is best illustrated with an example. Let's assume we have three states and they represent three consecutive geologic stages: 0 (Santonian, 85.8-83.5 Ma), 1 (Campanian, 83.5-70.6 Ma), and 2 (Maastrichtian, 70.6-65.5 Ma). Now the cost of a transition between states should be expressed in some difference between the ages of each state. In this example we will use the midpoints as our point estimate: 0 (84.65 Ma), 1 (77.05 Ma), and 2 (68.05 Ma). These must also be supplied to the function using the \code{state_ages} option: \code{state_ages = c(84.65, 77.05, 68.05)}. The resulting stepmatrix would look like this:
#'
#' \preformatted{    --------------------
#'     |  0  |  1  |   2  |
#' ------------------------
#' | 0 |  0  | 7.6 | 16.6 |
#' ------------------------
#' | 1 | Inf |  0  |   9  |
#' ------------------------
#' | 2 | Inf | Inf |   0  |
#' ------------------------}
#'
#' Note that the user needn't use midpoint (a first or last appearance date may be preferable). Indeeed, a simple unit count could be applied instead but this would then be identical to an irreversible character (except, perhaps, in the case of polymorphisms - see below).
#'
#' The use of stratigraphy as a character is probably not generally recommended, but is offered here for those that wish to explore stratocladistic methods that aim to incorporate stratigraphic information in phylogenetic inference and identify ancestor-descendant relationships (Fisher 1994).
#'
#' As with the other asymmetric characters it is an additional assumption that the root represent 0 (the oldest sampled value).
#'
#' \emph{Other character types}
#'
#' Although no other character types are supported here, other forms of stepmatrix can of course exist (e.g., see those in Hooker 2014). These can still be used elsewhere in Claddis, but users must either generate them themselves or have them be specified within a \#NEXUS file imported into Claddis. Note that general advice is to always have the diagonal be zero and the off-diagonal values be positive, otherwise downstream use will typically be confounded.
#'
#' \bold{Polymorphic characters}
#'
#' True polymorphisms represent perhaps the most complex problem for phylogenetic inference and ancestral state estimation (reconstruction), at least in part due to the varied nature their meaning can have (Swofford and Maddison 1992). For example, they may represent true variance within a species or composite variance in some higher taxon represented by a single OTU (something to be avoided if at all possible). Similarly, as different cladistic matrix formats do not necessarily allow the distinction between uncertainties and polymorphisms they can sometimes represent the former even though their encoding suggest the latter. For example, TNT (Goloboff et al. 2008; Goloboff and Catalano 2016) conflates the two with a single representation within square brackets, [].
#'
#' Misrepresenting true polymorphisms as uncertainties leads to incorrect estimates of the amount of evolution (tree length under parsimony; Nixon and Davis 1991) and disallows, potentially inappropriately, polymorphic values as ancestral state estimates (reconstructions). Hence, it can be argued polymorphic combinations should also be considered in stepmatrices, too (e.g., Maddison and Maddison 1987). This is allowed here (\code{include_polymorphisms = TRUE}), but is only approprate for some specific character types which we can consider in turn below.
#'
#' \emph{Unordered polymorphisms}
#'
#' Maddison and Maddison (1987) proposed the following approach for a three-state character (NB: this is modified slightly such that the single state transitions are scaled to one):
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
#' Here we see every possible polymorphism is accounted for in the rows and columns and a cost associated with each transition to, from, or between polymorphic states. The assignment of costs can be thought of as the presence-absence "switching" to go from one state to another, with each switch costing 0.5 steps. Thus to go from 0 to 0&1 we are making one "switch", turning state 1 "on", hence the total cost is just 0.5. Similarly, to go from 0 to state 1 we turn "on" state 1 and turn "off" state 0, costing us 0.5 + 0.5 = 1 steps.
#'
#' Assigning costs in this way can become tedious as the total number of states grows and fortunately there is another way to conceptualize this data. Specifically, a coordinate space with dimensionality equivalent to the number of single states which themselves represent orthogonal axes. So for our 0, 1, and 2 example we would have three-axes (three-dimensions): the 0-axis, the 1-axis, and the 2-axis. Then each state (single or polymorphic) can be plotted into the space as being "present" ("on" in the above analogy) or "absent" ("off"). So the state 0&1 would be present-present-absemt, or have the coordinates 1, 1, 0 (with values given in 0, 1, 2 order as a convention). In such a space, which would really form a cube, the states occupy each vertex (corner) and the Manhattan distances (shortest path when confined to traversing the edges of the cube) between them would represent the transition cost (\code{polymorphism_shape = "hypercube"} and \code{polymorphism_distance = "manhattan"}). Here each edge-length would be 0.5, but if you prefer you can consider the unit cube and then rescale at the end (which is what the function does in practice). (Note that the origin, i.e., the point 0, 0, 0, will always be unoccupied.) Of course, the use of such a coordinate space makes sense for an unordered character, but also suggests other possibilties for where polymorphic states should be positioned and how distances should be calculated. Although, importantly, any other candidate shapes or distances should make logical sense and be freely extensible to higher dimensions (larger numbers of states), just as the cube generalizes to the hypercube (e.g., in two-dimensions it is the square, and in four the tesseract). However, first I offer a reason to consider such alternatives.
#'
#' Let's imagine a "good" polymorphic stepmatrix is one that assigns costs such that polymorphisms are truly considered optimal solutions when the situation warrants it. We can use a contrived example to show this - the star phylogeny for three tips each with a different state:
#'
#' \preformatted{0   1   2
#'  \  |  /
#'   \ | /
#'    \|/
#'     X}
#'
#' We might imagine that the polymorphism 0&1&2 ought to be optimal (lowest cost) here as our value for the ancestor X. But if we use the stepmatrix from above we find the following costs:
#'
#' \preformatted{---------------------
#' | Ancestral | Total |
#' |   state   |  cost |
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
#' Thus, 0&1&2 is neither uniquely optimal, nor amongst the optimal solutions. In fact it is uniquely maximally sub-optimal. We could consider the same coordinate space as before, but assign transition costs based on Euclidean (straight-line) distances instead, but we find this doesn't change the rankings (distances are given to 1dp):
#'
#' \preformatted{---------------------
#' | Ancestral | Total |
#' |   state   |  cost |
#' ---------------------
#' |     0     |   2   |
#' |     1     |   2   |
#' |     2     |   2   |
#' |    0&1    |  2.6  |
#' |    0&2    |  2.6  |
#' |    1&2    |  2.6  |
#' |   0&1&2   |   3   |
#' ---------------------}
#'
#' Indeed, we would never favour any polymorphism here. Thus to allow for 0&1&2 to be optimal we must consider other "shapes" for our polymorphism space. This function offers two alternatives: 1. the simplex (the hyperdimensional equivalent of an equilateral triangle), and 2. the hypersphere (the hyperdimensional equivalent of the circle). These two shapes also imply their own distance measures (Euclidean and Great Circle, respectively). Both extend to any dimension as required and both have the critical property of bringing the coordinates of polymorphic combinations closer to the origin and hence making them more plausible candidates under the parsimony criterion.
#'
#' For the simplex (\code{polymorphism_shape = "simplex"}) the single states represent the vertices, with double states (e.g., 0&1) plotting midway along edges. More specifically, coordinates will always sum to one, e.g., for our three-state scenario, state 1 (0 + 1 + 0 = 1), state 1&2 (0 + 1/2 + 1/2 = 1), and state 0&1&2 (1/3 + 1/3 + 1/3 = 1), hence this is how they are assigned. Transition costs are then best considered as the Euclidean distances (\code{polymorphism_distance = "euclidean"}) between points. Repeating our scenario above this would lead to the following costs for each candidate for ancestor X (distances given to 1dp):
#'
#' \preformatted{---------------------
#' | Ancestral | Total |
#' |   state   |  cost |
#' ---------------------
#' |     0     |   2   |
#' |     1     |   2   |
#' |     2     |   2   |
#' |    0&1    |  1.8  |
#' |    0&2    |  1.8  |
#' |    1&2    |  1.8  |
#' |   0&1&2   |  1.7  |
#' ---------------------}
#'
#' For the hypersphere (\code{polymorphism_shape = "hypersphere"}) the single states represent the vertices, with double states (e.g., 0&1) plotting midway along the great circle linking them (i.e., the surface of a hypersphere). Now coordinates are given as the square root of 1/N along each axis the state is present on (0 if absent), and where N is the number of states present in the polymorphism. Thus, 0&2 would be given as sqrt(1/2), 0, sqrt(1/2). Transition costs are now best considered as the equivalent of a Great Circle distance - the shortest distance across the surface of the hypersphere (\code{polymorphism_distance = "great_circle"}) between two points. Repeating our scenario above again this would lead to the following costs for each candidate for ancestor X (distances given to 1dp):
#'
#' \preformatted{---------------------
#' | Ancestral | Total |
#' |   state   |  cost |
#' ---------------------
#' |     0     |   2   |
#' |     1     |   2   |
#' |     2     |   2   |
#' |    0&1    |   2   |
#' |    0&2    |   2   |
#' |    1&2    |   2   |
#' |   0&1&2   |  1.8  |
#' ---------------------}
#'
#' Note that this makes every candidate equally plausible, save 0&1&2 which is the uniquely maximally optimal solution.
#'
#' No recommendation is made here on which of these shape-distance combinations is "best", rather they are offered as potentially superior alternatives for a given situation. Note also, that although it is assumed the user will want to match the shape and distance as suggested above there are no restrictions preventing other combinations.
#'
#' \emph{Ordered polymorphisms}
#'
#' Initially the case for ordered characters may seem simple as they normally naturally form a linear series (one-dimensional space) and hence we can imagine inserting polymorphisms as intermediates along this same axis, e.g.:
#'
#' \preformatted{/---\  0.5  /---\  0.5  /---\  0.5  /---\  0.5  /---\
#' | 0 |<----->|0&1|<----->| 1 |<----->|1&2|<----->| 2 |
#' \---/       \---/       \---/       \---/       \---/}
#'
#' Where, again, costs are scaled such that those between single states are the same for non-polymorphic characters. This gives us the following step matrix:
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
#' However, it is less obvious what to do with the "non-intermediate" polymorphisms, 0&2 and 0&1&2, and here the one-dimensional space analogy breaks down. The only suggested solution of which I am aware is that of Maddison and Maddison (2000), where changes (0.5 steps here) represent adding OR losing a state (i.e.,  the on-off switch analogy given above). The difference this time is that the only switches that are "reachable" are those adjacent to ones already turned on. I.e., you cannot turn state 2 "on" if you currently only have state 0 turned on as the only reachable switch is state 1. Otherwise, the cost is like before, e.g., going from 0 to 0&1&2 involves (first) turning on state 1 (as this is already adjacent, i.e., accessible, to zero) then another 0.5 steps to add state 2 (now accessible via state 1). Extrapolating this to the complete stepmatrix gives us:
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
#' Note that in practical terms the function generates such stepmatrices by first establishing an adjacency matrix (which siwtches are reachable from the current state?) and then converting these to stepmatrices using the \link{convert_adjacency_matrix_to_stepmatrix} function.
#'
#' \emph{Dollo and irreversible polymorphisms}
#'
#' Dollo and irreversible (Camin-Sokal) characters suffer from the same problem when it comes to polymorphisms, in that polymorphic states confound stepmatrix construction. More specifically, for a Dollo character, the state 0&1 would suggest state 1 was acquired previously, but lost in some member(s) of the OTU. For an irreversible character the opposite would be true, and we would assume that state 1 was acquired independently by some member(s) of the OTU and state 0 was acquired previously. This nuance is not really codable in a stepmatrix as their construction doesn't allow for the correct information to be conveyed. Thus in practice these options are not allowed by the function. Instead a user may wish to encode the more important part of the polymorphism as a single state (e.g., 0&1 would become 1 for a Dollo character as this would better allow for correct placement of the acqusition of state 1 on the tree). However, this is non-ideal too as it would still lead to the problem of undercounting true changes that ignoring polymorphisms generates in the first place (see Nixon and Davis 1991, who also discuss alternative approaches to dealing with them).
#'
#' \emph{Stratigraphic polymorphisms}
#'
#' The stratigraphic case may seem like it ought to follow the same rule as irreversible characters, but this time polymorphisms have a different logical meaning. Specifically they would encode the presence of a taxon in multiple time units. However, this means that in practice (for estimating tree lengths or ancestral states) the taxon should simply be recoded as the oldest unit. E.g., 1&2&3 should just be recoded as 1. Thus, again, polymorphisms do not belong in a stratigraphc stepmatrix, however nothing is really lost by recoding them as a single state.
#'
#' \emph{Polymorphism limits}
#'
#' Even where polymorphisms may be appropriate (i.e., for some ordered and unordered characters) they still represent a major increase to stepmatrix size that can cause memory limits to be hit far in advance of stepmatrices that exclude them. Specifically, and as shown in our hypercube example, there will be as many possible states as N^2 - 1, where N is the number of single states and the minus one indicates the exclusion of the origin value. Thus the size of any stepmatrix required grows very quickly:
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
#' Because of this the function will become extremely slow for these higher values and indeed here is hard-capped at no more than fourteen states. Consequently any gap-weighted characters are also most likely inappropriate for polymorphism use.
#'
#' \bold{Character weights}
#'
#' Note that here here any stepmatrix (save the stratigraphic kind) is rescaled such that the "base weight" - which can be thought of as the lowest single state to single state transition possible (excluding the diagonal) - is set to 1. However, users may also wish to weight individual characters differently. It may be tempting, then, to simply multiply the values of each stepmatrix by that weight for use elsewhere. This is discouraged, however, as: 1. Other Claddis functions and data matrix data structures allow for weights to be encoded elsewhere (this would be at best redundant and at worst incorrect, if weights get applied twice), and 2. This can lead to confounded results in some cases. For example, if weighting a character zero to exclude it from an analysis any ancestral state estimation would become equally likely, meaning (inadvertently) large numbers of redundant calculations get made. In short, do not weight steomatrices direcrtly, but use the weights component of each matrix block in a \code{cladisticMatrix} object instead.
#'
#' \bold{Relation to Q-matrices}
#'
#' Q-matrices are used in likelihood and Bayesian approaches and also deal with state-to-state transitions, but are fundamentally distinct from stepmatrices in multiple ways. For example, Q-matrices encode transition rates, not costs, and their values are typically parameters that are estimated from the data not prior statements about the data. However, both can be rendered as Markov chain diagrams. Put simply though, they are not equivalent or interchangable.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Brocklehurst, N. and Field, D. J., 2021. Macroevolutionary dynamics of dentition in Mesozoic birds reveal no long-term selection towards tooth loss. \emph{iScience}, \bold{24}, 102243.
#'
#' Camin, J. H. and Sokal, R. R., 1965. A method for deducing branching sequences in phylogeny. \emph{Evolution}, \bold{19}, 311-26.
#'
#' Farris, J. S., 1970. Methods for computing Wagner trees. \emph{Systematic Zoology}, \bold{19}, 83-92.
#'
#' Fisher, D. C., 1994. Stratocladistics: morphological and temporal patterns and their relation to phylogenetic process. In L. Grande and O. Rieppel (eds.), \emph{Interpreting the Hierarchy of Nature}. Academic Press, San Diego. pp133–171.
#'
#' Fitch, W. M., 1971. Towards defining the course of evolution: minimum change for a specified tree topology. \emph{Systematic Zoology}, \bold{20}, 406-416.
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics/ \emph{Cladistics}, \bold{32}. 221-238
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' Hooker, J. J., 2014. New postcranial bones of the extinct mammalian family Nyctitheriidae (Paleogene, UK): primitive euarchontans with scansorial locomotion. \emph{Palaeontologia Electronica}, \bold{17.3.47A}, 1-82.
#'
#' Kitching, I. J., Forey, P. L., Humphries, C. J. and Williams, D. M., 1998. \emph{Cladistics: The Theory and Practice of Parsimony Analysis (Second Edition)}. Oxford University Press, Oxford. 229pp.
#'
#' Nixon, K. C. and Davis, J. I., 1991. Polymorphic taxa, missing values and cladistic analysis. \emph{Cladistics}, \bold{7}, 233-241.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. In R. L. Mayden (ed.), \emph{Systematics, Historical Ecology, and North American Freshwater Fishes}. Stanford University Press, Stanford. pp187-223.
#'
#' Wagner, W. H., 1961. Problems in the classification of ferns. \emph{Recent Advances in Botany}, \bold{1}, 841-844.
#'
#' @return
#'
#' An object of class \code{stepMatrix} - a square stepmatrix with rows representing "from" states, columns "to" states, and individual cells the cost in steps of that transition.
#'
#' @seealso
#'
#' \link{make_all_polymorphisms}
#'
#' @examples
#'
#' # Make an unordered stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make an unordered stepmatrix including polymorphisms:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = TRUE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make an ordered stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make an ordered stepmatrix including polymorphisms:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered",
#'   include_polymorphisms = TRUE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make a Dollo stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "dollo",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle",
#'   dollo_penalty = 100
#' )
#'
#' # Make an irreversible stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "irreversible",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make a stratigraphic stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "stratigraphy",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle",
#'   state_ages = c(52, 34, 12)
#' )
#'
#' @export make_stepmatrix
make_stepmatrix <- function(min_state = 0, max_state, character_type, include_polymorphisms = FALSE, polymorphism_shape, polymorphism_distance, state_ages, dollo_penalty = 999) {
  
  # - Need more incoming checks probably?
  # - Option to disallow polymorphisms as ancestors (still allows costs for transitions to polymorphic states to be taken into account)
  
  # Check character_type is a valid value and stop and warn user if not:
  if (length(x = setdiff(x = character_type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy"))) > 0) stop("character_type must be one of \"ordered\", \"unordered\", \"dollo\", \"irreversible\", or \"stratigraphy\".")
  
  # Check polymorphism_shape is a valid value and stop and warn user if not:
  if (length(x = setdiff(x = polymorphism_shape, y = c("hypercube", "hypersphere", "simplex"))) > 0) stop("character_type must be one of \"hypercube\", \"hypersphere\", \"simplex\".")
  
  # Check polymorphism_distance is a valid value and stop and warn user if not:
  if (length(x = setdiff(x = polymorphism_distance, y = c("manhattan", "euclidean", "great_circle"))) > 0) stop("character_type must be one of \"manhattan\", \"euclidean\", \"great_circle\".")
  
  # Get single states:
  single_states <- min_state:max_state
  
  # Case if character is ordered:
  if (character_type == "ordered") {
    
    # Store symmetry of stepmatrix as symmetric:
    symmetry <- "Symmetric"
    
    # If including polymorphisms:
    if (include_polymorphisms) {
      
      # Generate and store all possible states:
      all_states <- make_all_polymorphisms(single_states = single_states)
      
      # Check data are not too big (>= 2^14 states) and stop and warn user if so:
      if (length(x = all_states) >= 16384) stop("Stepmatrix would be too large. Use fewer states.")
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = all_states)
      
      # Initialise adjacency atrix with zeroes:
      adjacency_matrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(all_states, all_states))
      
      # For each state save the last one:
      for(state in all_states[-length(x = all_states)]) {
        
        # Find matching states for current state:
        state_matches <- apply(X = do.call(what = rbind, args = lapply(X = as.list(x = strsplit(x = state, split = "&")[[1]]), FUN = function(x) single_states == x)), MARGIN = 2, FUN = any)
        
        # Find states that are addable (reachable) from current state(s):
        addable_states <- single_states[apply(X = rbind(!state_matches, unlist(lapply(X = as.list(x = 1:length(x = state_matches)), FUN = function(x) any(x = c(state_matches[x - 1], state_matches[x + 1]), na.rm = TRUE)))), MARGIN = 2, FUN = all)]
        
        # Compile state(s) adjacent to current state:
        adjacent_states <- unlist(x = lapply(X = as.list(x = addable_states), FUN = function(new_state) paste(x = sort(x = c(new_state, strsplit(x = state, split = "&")[[1]])), collapse = "&")))
        
        # Populate adjaceney matrix with these states:
        adjacency_matrix[state, adjacent_states] <- adjacency_matrix[adjacent_states, state] <- 1
        
      }
      
      # Convert adjacency matrix to step matrix:
      stepmatrix <- convert_adjacency_matrix_to_stepmatrix(adjacency_matrix = adjacency_matrix)
      
      # Rescale such that adjacent single states are one step:
      stepmatrix <- stepmatrix / stepmatrix[1, 2]
    
    # If excluding polymorphisms:
    } else {
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = single_states)
      
      # Initialise stepmatrix with all values set to zero:
      stepmatrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(single_states, single_states))
      
      # Store ordered values in upper and lower triangles:
      stepmatrix[upper.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1))
      stepmatrix[lower.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = (stepmatrix_size - 1):1), FUN = function(steps) 1:steps))

    }
  
  }

  # Case if character is unordered:
  if (character_type == "unordered") {
    
    # Store symmetry of stepmatrix as symmetric:
    symmetry <- "Symmetric"

    # If including polymorphisms:
    if (include_polymorphisms) {
      
      # Generate and store all possible states:
      all_states <- make_all_polymorphisms(single_states = single_states)
      
      # Check data are not too big (>= 2^14 states) and stop and warn user if so:
      if (length(x = all_states) >= 16384) stop("Stepmatrix would be too large. Use fewer states.")
      
      # Create coordinate matrix and initialise with zeroes:
      state_presence_matrix <- matrix(data = 0, nrow = length(x = single_states), ncol = length(x = all_states), dimnames = list(single_states, all_states))
      
      # If using the hypercube coordinate space assign coordinates accordingly:
      if (polymorphism_shape == "hypercube") for(i in 1:ncol(x = state_presence_matrix)) state_presence_matrix[strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]], i] <- 1
      
      # If using the hypersphere or simplex coordinate space:
      if (polymorphism_shape == "hypersphere" || polymorphism_shape == "simplex") {
        
        # For each coding:
        for(i in 1:ncol(x = state_presence_matrix)) {
          
          # Isolate components of polymorphism:
          components <- strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]]
          
          # If using the hypersphere coordinate space then apply coordinates accordingly (the square root of 1/N states in polymorphism on each axis state is present):
          if (polymorphism_shape == "hypersphere") state_presence_matrix[components, i] <- sqrt(x = 1 / length(x = components))
          
          # If using the simplex coordinate space then apply coordinates accordingly (1/N states in polymorphism on each axis state is present):
          if (polymorphism_shape == "simplex") state_presence_matrix[components, i] <- 1 / length(x = components)
          
        }
        
      }
      
      # If using a manhattan or euclidean distance, calculate distance directly from coordinate-space:
      if (polymorphism_distance == "euclidean" || polymorphism_distance == "manhattan") stepmatrix <- as.matrix(x = dist(x = t(x = state_presence_matrix), method = polymorphism_distance, diag = TRUE, upper = TRUE))
      
      # If using a great circle distance:
      if (polymorphism_distance == "great_circle") {
        
        # Start by calculating the euclidean distances between each point:
        stepmatrix <- as.matrix(x = dist(x = t(x = state_presence_matrix), method = "euclidean", diag = TRUE, upper = TRUE))
        
        # Treating distances as the chord length of a circle of radius one transform them to the arc length for the same:
        stepmatrix <- 2 * asin(x = stepmatrix / 2)
        
      }
      
      # Rescale such that single state distances (e.g., 0 to 1) are one (i.e., a normal unordered character):
      stepmatrix <- stepmatrix / stepmatrix[1, 2]
      
    # If excluding polymorphisms:
    } else {
      
      # Initialise stepmatrix with all values set to one:
      stepmatrix <- matrix(data = 1, nrow = length(x = single_states), ncol = length(x = single_states), dimnames = list(single_states, single_states))
      
      # Make diagional of stepmatrix zero:
      diag(x = stepmatrix) <- 0
      
    }
    
  }
  
  # Case if character is Dollo:
  if (character_type == "dollo") {
    
    # Store symmetry of stepmatrix as asymmetric:
    symmetry <- "Asymmetric"

    # If including polymorphisms:
    if (include_polymorphisms) {
      
      # Stop and warn user:
      stop("If character_type is \"dollo\" then include_polymorphisms cannot be TRUE. This is because there is no safe way to record acquisitions and losses of single states where polymorphisms exist using just a stepmatrix approach. Consider recoding as a single state, although note that this will necessarily mean some change(s) are not accounted for.")
      
    # If excluding polymorphisms:
    } else {
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = single_states)
      
      # Initialise stepmatrix with all values set to zero:
      stepmatrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(single_states, single_states))
      
      # Set upper triangle as an ordered character times the dollo_penalty
      stepmatrix[upper.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1)) * dollo_penalty
      
      # Store regular ordered values in lower triangles::
      stepmatrix[lower.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = (stepmatrix_size - 1):1), FUN = function(steps) 1:steps))
      
    }
    
  }
  
  # Case if character is irreversible:
  if (character_type == "irreversible") {
    
    # Store symmetry of stepmatrix as asymmetric:
    symmetry <- "Asymmetric"

    # If including polymorphisms:
    if (include_polymorphisms) {
      
      # Stop and warn user:
      stop("If character_type is \"irreversible\" then include_polymorphisms cannot be TRUE. This is because there is no safe way to record acquisitions and losses of single states where polymorphisms exist using just a stepmatrix approach. Consider recoding as a single state, although note that this will necessarily mean some change(s) are not accounted for.")
      
    # If excluding polymorphisms:
    } else {
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = single_states)
      
      # Initialise stepmatrix with all values set to zero:
      stepmatrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(single_states, single_states))
      
      # Exclude losses by setting lower triangle values to infinity:
      stepmatrix[lower.tri(x = stepmatrix)] <- Inf
      
      # Set upper triangle as an ordered character:
      stepmatrix[upper.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1))
      
    }
    
  }
  
  # Case if character is stratigraphy:
  if (character_type == "stratigraphy") {
    
    # Store symmetry of stepmatrix as asymmetric:
    symmetry <- "Asymmetric"

    # Add state names to ages:
    names(x = state_ages) <- min_state:max_state
    
    # If including polymorphisms:
    if (include_polymorphisms) stop("If character_type is \"stratigraphy\" then include_polymorphisms cannot be TRUE. If the age of the OTU is uncertain then code it as such (use / instead of & between states). If the OTU is truly present in multiple units then only code it as the oldest one, or, alternatively, break the OTU up into multiple OTUs.")
    
    # Generate initial stepmatrix using temporal distances:
    stepmatrix <- as.matrix(x = dist(x = state_ages, diag = TRUE, upper = TRUE))
    
    # Set infinite cost to all reversals:
    stepmatrix[lower.tri(x = stepmatrix)] <- Inf
    
  }
  
  # Create full stepmatrix object:
  stepmatrix <- list(size = ncol(x = stepmatrix), stepmatrix = stepmatrix, symmetry = symmetry, includes_polymorphisms = include_polymorphisms)
  
  # Set class of output as stepMatrix:
  class(stepmatrix) <- "stepMatrix"
  
  # Return stepmatrix:
  stepmatrix
  
}
