GetAllStateChanges <- function(clad.matrix, tree, time.bins, Nsim=10) {

# The above could maybe be modified to allow other options to be specified and handed to make.simmap (pi, Q etc.)
# ADD CONDITIONALS HERE TO CHECK DATA ARE VALID
# NEED TO TIME BIN INTERNAL VERSUS TERMINAL BRANCHES; ALSO MAYBE BY CHARACTER?
# MAKE LESS VERBOSE OUTPUT?
	
# NEED TO SPLIT APART TERMINAL AND INTERNAL EDGES IN BINS:
#terminal.edges <- match(1:Ntip(tree), tree$edge[, 2])
#internal.edges <- setdiff(1:nrow(tree$edge), terminal.edges)
	
	# Create strings that define character models (non-missing states plus ordering):
	char.model.sets <- paste(unlist(lapply(lapply(lapply(lapply(apply(clad.matrix$matrix, 2, strsplit, split="&"), unlist), sort), unique), paste, collapse="")), clad.matrix$ordering, sep="")
	
	# Create empty list to store character models:
	char.models <- list()
	
	# For each unique character model:
	for(i in unique(char.model.sets)) {
		
		# Get the ordering type:
		ordering.type <- gsub("[0-9]", "", i)
		
		# Get the unique non-missing states:
		states <- as.numeric(strsplit(strsplit(i, "[A-z:a-z]")[[1]][1], "")[[1]])
		
		# Create an empty square matrix of zeroes that will serve as the base for the character model:
		char.model <- matrix(0, nrow=length(states), ncol=length(states))

		# Only continue if there are enough states to alter the character model:
		if(length(states) > 1) {
			
			# Check for confounding characters (those which are unordered, but not all states between the minimum and the maximum are sampled):
			if((diff(range(states)) + 1) != length(states) && ordering.type == "ord") {
				
				# Stop and give the user an error message:
				stop(paste("ERROR: cannot currently deal with the following ordered character(s) where not all states are sampled amongst tips: ", paste(which(char.model.sets == i), collapse=" "), sep=""))
				
			}
			
			# If character is unordered:
			if(ordering.type == "unord") {
				
				# Make all off-diagonal elements equal to 1:
				for(j in 1:nrow(char.model)) char.model[setdiff(1:nrow(char.model), j), j] <- 1
				
			}
			
			# If character is ordered:
			if(ordering.type == "ord") {
				
				# Make only immediate off-diagonal elements equal to 1:
				for(j in 1:nrow(char.model)) char.model[setdiff(1:nrow(char.model), j)[abs(setdiff(1:nrow(char.model), j) - j) == 1], j] <- 1
				
			}
			
			# If character is neither ordered or unordered:
			if(ordering.type != "unord" && ordering.type != "ord") {
				
				# Simply send an error message for now:
				stop(paste("ERROR: the following character(s) are neither ordered nor unordered and no model can be specified: ", paste(which(char.model.sets == i), collapse=" "), sep=""))
				
# NEED MORE HERE FOR CONTINUOUS AND STEP MATRIX CHARACTERS; CAN THEORETICALLY DO LATTER, BUT NOT SURE HOW CAN DO FORMER YET - IS THERE A WAY TO POINT SAMPLE A BROWNIAN MOTION MODEL? EVEN IF YES, WHEN IS IT SAMPLED? HAS DOWNSTREAM EFFECTS ON TIME BIN CHOICES!
				
			}
			
			# Store character model for each character that has that model:
			for(j in which(char.model.sets == i)) char.models[[j]] <- char.model
			
		}
		
	}
	
	# Create empty matrix to store all state changes:
	allstatechanges <- matrix(nrow=0, ncol=6, dimnames=list(c(), c("Character", "Nsim", "Edge", "Age", "FromState", "ToState")))

	# Vector to store time (root denominator) for each character and simulation:
	chartime <- rep(0, ncol(clad.matrix$matrix))
	
	# Matrix to store time for each character:
	terminal.edge.lengths.per.bin <- internal.edge.lengths.per.bin <- edge.lengths.per.bin <- matrix(0, ncol=length(time.bins) - 1, nrow=ncol(clad.matrix$matrix))
	
	# Get node ages for original tree:
	tree.nodeages <- GetNodeAges(tree)
	
	# Get binary strings for missing data distributions for each character:
	missing.strings <- apply(matrix(as.numeric(is.na(clad.matrix$matrix)), nrow=nrow(clad.matrix$matrix)), 2, paste, collapse="")
	
	# Create empty list to store characters that share the saem distribution of missing data:
	MissingStringCharacters <- list()
	
	# Find the characters that share the same distribution of missing data:
	for(i in 1:length(unique(missing.strings))) MissingStringCharacters[[i]] <- which(missing.strings == unique(missing.strings)[i])
	
	# For each character that has the same distribution of missing data:
	for(i in 1:length(unique(missing.strings))) {
		
		# Find taxa to drop to make pruned tree:
		taxa.to.drop <- rownames(clad.matrix$matrix)[which(strsplit(unique(missing.strings)[i], "")[[1]] == "1")]
		
		# Case if there are too few taxa to return a meaningful pruned tree:
		if((length(tree$tip.label) - length(taxa.to.drop)) < 3) {
			
			# If no taxa are scored for the character(s):
			if((length(tree$tip.label) - length(taxa.to.drop)) == 0) {
			
				# Warn user about zero taxon scoring:
				cat(paste("Warning: the following character(s) are not scored for any taxa and hence no changes are recorded: ", paste(MissingStringCharacters[[i]], collapse=" "), "\n", sep=""))

				# Add information on root state (NA) to allstatechanges (no need to update time fields as these are simply the default of zero):
				allstatechanges <- rbind(allstatechanges, cbind(rep(MissingStringCharacters[[i]], Nsim), sort(rep(1:Nsim, length(MissingStringCharacters[[i]]))), rep(0, Nsim * length(MissingStringCharacters[[i]])), rep(tree$root.time, Nsim * length(MissingStringCharacters[[i]])), rep(NA, Nsim * length(MissingStringCharacters[[i]])), rep(NA, Nsim * length(MissingStringCharacters[[i]]))))
			
			}
			
			# Case if only one taxon is scored for the character(s):
			if((length(tree$tip.label) - length(taxa.to.drop)) == 1) {

				# Warn user about single taxon scoring:
				cat(paste("Warning: the following character(s) are only scored for one taxon and hence no changes are recorded: ", paste(MissingStringCharacters[[i]], collapse=" "), "\n", sep=""))
				
				# Get terminal node number:
				tip.number <- match(setdiff(tree$tip.label, taxa.to.drop), tree$tip.label)
				
				# Get edge number:
				edge.number <- match(tip.number, tree$edge[, 2])
				
				# Get beginning and end of branch:
				branch.end.dates <- tree.nodeages[tree$edge[edge.number, ]]
				
				# Get change times:
				change.times <- runif(length(MissingStringCharacters[[i]]) * length(MissingStringCharacters[[i]]), min=branch.end.dates[2], max=branch.end.dates[1])
				
				# For each character:
				for(j in MissingStringCharacters[[i]]) {
					
					# Get states for terminal taxa:
					terminal.states <- clad.matrix$matrix[setdiff(tree$tip.label, taxa.to.drop), MissingStringCharacters[[i]]]
				
					# Expand by number of simulations requested:
					terminal.states <- rep(terminal.states, Nsim)
				
					# Conditional if there are polymorphisms:
					if(length(grep("&", terminal.states)) > 0) {
				
						# Get polymorhisms indices:
						polymorphisms <- grep("&", terminal.states)
					
						# Replace polymorphisms with any included state at random:
						for(k in polymorphisms) terminal.states[k] <- sample(strsplit(terminal.states[k], "&")[[1]])[1]
				
					}
				
					# Convert to numerics for later storing:
					terminal.states <- as.numeric(terminal.states)
				
					# Update all state changes:
					allstatechanges <- rbind(allstatechanges, cbind(rep(j, Nsim), c(1:Nsim), rep(edge.number, Nsim), rep(change.times[match(j, MissingStringCharacters[[i]])], Nsim), rep(NA, Nsim), terminal.states))
					
					# Create pruned tree for character:
					pruned.tree <- tree
				
					# Set all branch lengths to zero:
					pruned.tree$edge.length <- rep(0, nrow(pruned.tree$edge))
				
					# Update root age::
					pruned.tree$root.time <- change.times[match(j, MissingStringCharacters[[i]])]
					
					# Add edge length for character:
					pruned.tree$edge.length[edge.number] <- pruned.tree$root.time - branch.end.dates[2]
					
					# Get edge length per bin:
					edge.length.per.bin <- EdgeLengthsInBins(pruned.tree, time.bins)
					
					# Store edge length per bin for character:
					edge.lengths.per.bin[j, ] <- edge.length.per.bin$edge.length.in.bin
					
					# Store terminal edge length per bin for character:
					terminal.edge.lengths.per.bin[j, ] <- edge.length.per.bin$terminal.edge.length.in.bin
					
					# Store internal edge length per bin for character:
					internal.edge.lengths.per.bin[j, ] <- edge.length.per.bin$internal.edge.length.in.bin
					
					# Update character time:
					chartime[j] <- sum(pruned.tree$edge.length)
					
				}
			
			}
			
			# Case if only two taxa are scored for the character(s):
			if((length(tree$tip.label) - length(taxa.to.drop)) == 2) {
				
				# Make miniature matrix for just the taxa and characters involved:
				scored.matrix <- matrix(clad.matrix$matrix[-match(taxa.to.drop, rownames(clad.matrix$matrix)), MissingStringCharacters[[i]]], nrow=2, dimnames=list(setdiff(rownames(clad.matrix$matrix), taxa.to.drop), MissingStringCharacters[[i]]))
				
				# Get root node for two taxon tree:
				root.node <- FindAncestor(rownames(scored.matrix), tree)
				
				# Get tip numbers for taxa in tree:
				tip.numbers <- match(rownames(scored.matrix), tree$tip.label)
				
				# Set up intitial list of edge numbers:
				edge.numbers <- match(tip.numbers, tree$edge[, 2])
				
				# While edges are incomplete (do not connect both taxa all the way to the root):
				while(sum(tree$edge[edge.numbers, 1] == root.node) != 2) {
					
					# Add next edge down for al edges not already connected to root:
					edge.numbers <- sort(unique(c(edge.numbers, match(tree$edge[edge.numbers[which(tree$edge[edge.numbers, 1] != root.node)], 1], tree$edge[, 2]))))
					
				}
				
				# Set up pruned tree from full tree:
				pruned.tree <- tree
				
				# make all non-sampled eedges zero in length:
				pruned.tree$edge.length[setdiff(1:length(tree$edge.length), edge.numbers)] <- 0
				
				# Update root time for pruned tree:
				pruned.tree$root.time <- tree.nodeages[root.node]
				
				# Store rate denominator (total edge time) for character(s):
				chartime[MissingStringCharacters[[i]]] <- sum(pruned.tree$edge.length)

				# Get edge length per bin for pruned tree:
				edge.length.per.bin <- EdgeLengthsInBins(pruned.tree, time.bins)
				
				# For each character:
				for(j in MissingStringCharacters[[i]]) {
					
					# Store edge length per bin:
					edge.lengths.per.bin[j, ] <- edge.length.per.bin$edge.length.in.bin
					
					# Store terminal edge length per bin:
					terminal.edge.lengths.per.bin[j, ] <- edge.length.per.bin$terminal.edge.length.in.bin
					
					# Store internal edge length per bin:
					internal.edge.lengths.per.bin[j, ] <- edge.length.per.bin$internal.edge.length.in.bin
				
				}
				
				# For each character:
				for(j in MissingStringCharacters[[i]]) {
				
					# Case if character is variable:
					if(length(unique(unlist(strsplit(scored.matrix[, as.character(j)], "&")))) > 1) {
				
						# Warn user about two taxon scoring:
						cat(paste("Warning: the following character is only scored for two taxa and hence the minimum number of changes are recorded: ", j, "\n", sep=""))
				
					# Case if character is constant:
					} else {
						
						# Warn user about two taxon scoring:
						cat(paste("Warning: the following character is constant and hence no changes are recorded: ", j, "\n", sep=""))
						
						
# Record root value and change(s) from scored to NA (trickier!)
						
						
					}
						
				}
		
# Given it might be easiest to set up all.charmaps as for trees with > 2 taxa maybe this should be remvoed from the current conditional and placed into the one below?
# AH, BUT THIS WILL NOT QUITE WORK AS PRUNED TREE HAS ALL THE EDGES, JUST LENGTHS SET TO ZERO. HMMM.
# And deal with polymorphisms!!!!!!
# root.state <- ??????
# Not sure what changes to record!!!!
# Repeat Nsim times
				
			}
		
		# Can proceed to stochastic character mapping as there are sufficient taxa:
		} else {
			
			# Prune tree to just taxa that can be scored for those character(s):
			pruned.tree <- drop.tip(tree, taxa.to.drop)

			# Ensure pruned trees $root.time value is correct:
			pruned.tree <- CorrectRootTime(tree, pruned.tree)
			
			# Get node ages for pruned tree:
			pruned.nodeages <- GetNodeAges(pruned.tree)
			
			# Store rate denominator (total edge time) for character(s):
			chartime[MissingStringCharacters[[i]]] <- sum(pruned.tree$edge.length)
			
			# Get ages for start of each branch:
			ageatstartofbranch <- pruned.nodeages[pruned.tree$edge[, 1]]

			# Create matrix to store root state for all characters (columns) and simulations (rows) with same pruned.tree:
			root.states <- matrix(NA, ncol=length(MissingStringCharacters[[i]]), nrow=Nsim)
			
			# Create empty list to store character maps for all characters with same pruned.tree:
			all.charmaps <- vector("list", length(MissingStringCharacters[[i]]))
			
			# Add character numbers as names for character maps:
			names(all.charmaps) <- MissingStringCharacters[[i]]
			
			# For each character:
			for(j in MissingStringCharacters[[i]]) {
				
				# Make this into prob matrix to hand to make.simmap (i.e., deal with polymorphisms):
				tip.values <- clad.matrix$matrix[pruned.tree$tip.label, j]
				
				# Specify character model:
				char.model <- char.models[[j]]
				
				# If character varies:
				if(length(sort(unique(unlist(strsplit(tip.values, "&"))))) > 1) {
					
# NEED SOMETHING HERE TO CATCH ISSUE WITH ORDERED CHARACTERS FOR WHICH NOT ALL STATES ARE RECOVERED
# ALSO WHAT IF NEITHER ORDERED NOR UNORDERED?
					
					# Create tip states matrix that will serve as priors for make.simmap (if character is ordered):
					if(clad.matrix$ordering[j] == "ord") tip.states <- matrix(0, nrow=Ntip(pruned.tree), ncol=length(sort(unique(unlist(strsplit(tip.values, "&"))))), dimnames=list(pruned.tree$tip.label, range(as.numeric(sort(unique(unlist(strsplit(tip.values, "&"))))))[1]:range(as.numeric(sort(unique(unlist(strsplit(tip.values, "&"))))))[2]))

					# Create tip states matrix that will serve as priors for make.simmap (if character is unordered):
					if(clad.matrix$ordering[j] == "unord") tip.states <- matrix(0, nrow=Ntip(pruned.tree), ncol=length(sort(unique(unlist(strsplit(tip.values, "&"))))), dimnames=list(pruned.tree$tip.label, sort(unique(unlist(strsplit(tip.values, "&"))))))
					
					# Case if polymorphisms amongst tip states:
					if(length(grep("&", tip.values))) {
					
						# For each polymorphism:
						for(k in grep("&", tip.values)) {
						
							# Update tip states for each state in polymorphism (each considered equally likely):
							tip.states[k, strsplit(tip.values[k], "&")[[1]]] <- 1 / length(strsplit(tip.values[k], "&")[[1]])
						
						}
						
						# Update tip values without polymorphisms:
						tip.values <- tip.values[-grep("&", tip.values)]
					
					}
					
					# Fill in tip states matrix for non-polymorphic codings:
					for(k in colnames(tip.states)) tip.states[names(which(tip.values == k)), k] <- 1
					
					# Do stochastic character mapping Nsim times (using pi="estimated" allows root state to vary, but other options should maybe be allowed):
					trees <- make.simmap(pruned.tree, tip.states, nsim=Nsim, model=char.model, pi="estimated", message=FALSE)
					
					# If Nsim is 1 need to convert to multiPhylo object so handling below works out:
					if(class(trees) != "multiPhylo") {
						
						# Make trees into list:
						trees <- list(trees)
						
						# Make into multiPhylo object:
						class(trees) <- "multiPhylo"
						
					}
					
					# Get just the character maps for each Nsim:
					all.charmaps[[match(j, MissingStringCharacters[[i]])]] <- charmaps <- lapply(trees, '[[', "maps")
					
					# Clean up memory a bit by removing trees (could be large variable if Nsim is large):
					rm(trees); gc()
					
					# For each simulation performed:
					for(k in 1:Nsim) {
					
						# Get root state(s) for pruned.tree:
						root.states[k, match(j, MissingStringCharacters[[i]])] <- root.state <- as.numeric(names(charmaps[[k]][[match(Ntip(pruned.tree) + 1, pruned.tree$edge[, 1])]][1]))
						
						# Find all edges with changes (no point looking at others):
						edgeswithchanges <- which(unlist(lapply(charmaps[[k]], length)) > 1)
						
						# For each branch with changes on it:
						for(l in edgeswithchanges) {
							
							# Get date for each change:
							ageofchanges <- ageatstartofbranch[l] - cumsum(charmaps[[k]][[l]])

							# Get time change happened:
							statechangetime <- ageofchanges[1:(length(ageofchanges) - 1)]
							
							# Get starting state of change:
							startstate <- as.numeric(names(ageofchanges[1:(length(ageofchanges) - 1)]))
							
							# Get finishing state of change:
							endstate <- as.numeric(names(ageofchanges[2:length(ageofchanges)]))
							
							# Add to list of all state changes:
							allstatechanges <- rbind(allstatechanges, cbind(rep(j, length(endstate)), rep(k, length(endstate)), rep(l, length(endstate)), statechangetime, startstate, endstate))
							
						}
						
					}
				
				# If character is constant:
				} else {
					
					# Warn user about constant character:
					cat(paste("Warning: the following character is constant and hence no changes are recorded: ", j, "\n", sep=""))
					
					# Set root states:
					for(k in 1:Nsim) root.states[k, match(j, MissingStringCharacters[[i]])] <- root.state <- as.numeric(unique(tip.values))

					# Create single instance of charmaps (no changes):
					charmaps <- as.list(pruned.tree$edge.length)
					
					# Add name of constant state to all values:
					for(k in 1:length(charmaps)) names(charmaps[[k]]) <- root.state
					
					# Create list to store Nsim repeated charmaps:
					charmaps2 <- list()
					
					# Store charmaps repeated for each Nsim:
					for(k in 1:Nsim) charmaps2[[k]] <- charmaps
					
					# Get just the character maps for each Nsim:
					all.charmaps[[match(j, MissingStringCharacters[[i]])]] <- charmaps <- charmaps2
					
				}
			
			}

			# Match edges between pruned character tree and complete original tree:
			edgematches <- EdgeMatch(tree, pruned.tree)

			# Get rows that correspond to current characters only:
			char.rows <- sort(unlist(lapply(lapply(as.list(MissingStringCharacters[[i]]), '==', allstatechanges[, "Character"]), which)))
			
			# Find edges in character tree that correspond to multiple edges in original tree:
			edgeswithmultiplematches <- as.numeric(names(which(unlist(lapply(edgematches$matching.edges[unique(allstatechanges[char.rows, "Edge"])], length)) > 1)))
			
			# Find edges in character tree that correspond to single edges in original tree:
			edgeswithsinglematches <- as.numeric(names(which(unlist(lapply(edgematches$matching.edges[unique(allstatechanges[char.rows, "Edge"])], length)) == 1)))
			
			# Get single edge match replacemnt edge numbers:
			singlematchreplacements <- unlist(edgematches$matching.edges[edgeswithsinglematches])

			# Create vector of length one to store edge replacements:
			edgereplacements <- NA
			
			# If there is more than one edge to replace create longer vector of NAs to store edge replacements:
			if(length(char.rows) > 1) edgereplacements <- rep(NA, nrow(allstatechanges[char.rows, ]))
			
			# As long as there are single edge replacements:
			if(length(edgeswithsinglematches) > 0) {
				
				# Find single match replacements and update them:
				for(j in edgeswithsinglematches) edgereplacements[which(allstatechanges[char.rows, "Edge"] == j)] <- singlematchreplacements[match(j, edgeswithsinglematches)]
				
			}
			
			# As long as there are multiple replacement edges:
			if(length(edgeswithmultiplematches) > 0) {
				
				# For each edge with multiple replacements:
				for(j in edgeswithmultiplematches) {
					
					# Get matching edges as edge numbers:
					matching.edges <- edgematches$matching.edges[[as.character(j)]]
					
					# Get matching edges as to-from node numbers:
					matching.edges.matrix <- tree$edge[matching.edges, ]
					
					# Get node ages for matching edges:
					matching.edges.ages <- cbind(tree.nodeages[matching.edges.matrix[, 1]], tree.nodeages[matching.edges.matrix[, 2]])
					
					# Get rows from allstatechanges which will need to be changed:
					rowstochange <- which(allstatechanges[char.rows, "Edge"] == j)
					
					# For each matching edge:
					for(k in 1:nrow(matching.edges.matrix)) {
						
						# Get maximum edge age:
						max.edge.age <- matching.edges.ages[k, 1]
						
						# Get minimum edge age:
						min.edge.age <- matching.edges.ages[k, 2]
						
						# Get changes on current edge:
						changes.on.current.edge <- intersect(which(allstatechanges[char.rows[rowstochange], "Age"] > min.edge.age), which(allstatechanges[char.rows[rowstochange], "Age"] <= max.edge.age))
						
						# Update edge replacements:
						edgereplacements[which(allstatechanges[char.rows, "Edge"] == j)[changes.on.current.edge]] <- matching.edges[k]
						
					}
					
				}
				
			}
			
			# Update state changes with edges from complete tree:
			allstatechanges[char.rows, "Edge"] <- edgereplacements

			# Anything with a row sum of zero in this list can be ignored:
			edgelinks <- matrix(FindLinkedEdges(tree)[edgematches$removed.edges, -edgematches$removed.edges], nrow=length(edgematches$removed.edges), dimnames=list(edgematches$removed.edges, colnames(FindLinkedEdges(tree))[-edgematches$removed.edges]))
			
			# Get list of missing edges that will have a node on a sampled edge:
			missingedges <- as.numeric(rownames(edgelinks)[which(apply(edgelinks, 1, sum) > 0)])
			
			# Case if missing edge terminates at root of pruned tree:
			if(length(sort(match(tree$edge[missingedges, 2], setdiff(tree$edge[sort(unlist(edgematches$matching.edges)), 1], tree$edge[sort(unlist(edgematches$matching.edges)), 2])))) > 0) {
				
				# Find node number of root in original tree:
				pruned.root.node <- setdiff(tree$edge[sort(unlist(edgematches$matching.edges)), 1], tree$edge[sort(unlist(edgematches$matching.edges)), 2])

				# Find root terminating edge:
				root.terminating.edge <- missingedges[match(pruned.root.node, tree$edge[missingedges, 2])]
				
				# Add root edge to allstatechanges:
				allstatechanges <- rbind(allstatechanges, cbind(sort(rep(MissingStringCharacters[[i]], Nsim)), rep(1:Nsim, ncol(root.states)), rep(root.terminating.edge, length(root.states)), rep(pruned.tree$root.time, length(root.states)), rep(NA, length(root.states)), as.vector(root.states)))
				
				# Update missing edges (remove edge prior to root or pruned tree):
				missingedges <- setdiff(missingedges, root.terminating.edge)
				
			# If root is sampled in pruned tree:
			} else {
				
				# Ensure root state (on edge "0") is recorded even if no other changes occur:
				allstatechanges <- rbind(allstatechanges, cbind(sort(rep(MissingStringCharacters[[i]], Nsim)), rep(1:Nsim, ncol(root.states)), rep(0, length(root.states)), rep(tree$root.time, length(root.states)), rep(NA, length(root.states)), as.vector(root.states)))
				
			}
			
			# List nodes from which missing edges emerge:
			missingedgestartnodes <- tree$edge[missingedges, 1]
			
			# Find edges at which missing start node is terminal to sample:
			edgestosample <- match(missingedgestartnodes, tree$edge[, 2])
			
			# For each edge to sample:
			for(j in edgestosample) {
				
				# Find pruned edge at which node emerges:
				pruned.edge <- names(which(unlist(lapply(lapply(edgematches$matching.edges, '==', j), sum)) == 1))
				
				# Get edges that make up pruned edge sampled:
				edge.string <- edgematches$matching.edges[pruned.edge][[1]]
				
				# For each character:
				for(k in MissingStringCharacters[[i]]) {
					
					# For each simulation:
					for(l in 1:Nsim) {
						
						# Get character map for edge:
						edge.charmap <- all.charmaps[[as.character(k)]][[l]][as.numeric(pruned.edge)][[1]]
						
						# Case if only a single state sampled on the edge (simple):
						if(length(edge.charmap) == 1) {
							
							# Get start state (only state on edge):
							start.state <- as.numeric(names(edge.charmap))
							
						}
						
						# Case if state changes along edge:
						if(length(edge.charmap) > 1) {
							
							# Get sample time for target node (end of target edge):
							sample.time <- cumsum(tree$edge.length[edge.string])[match(j, edge.string)]
							
							# Get state times:
							state.times <- cumsum(edge.charmap)
							
							# Get start state:
							start.state <- as.numeric(names(state.times[length(which(sample.time > state.times)) + 1]))
							
						}

						# Add shift to NA to allstatechanges:
						allstatechanges <- rbind(allstatechanges, c(k, l, missingedges[match(j, edgestosample)], tree.nodeages[missingedgestartnodes[match(j, edgestosample)]], start.state, NA))
						
					}
					
				}
				
			}
			
			# Get edge length per bin for pruned tree:
			edge.length.per.bin <- EdgeLengthsInBins(tree, time.bins, pruned.tree)
			
			# Store edge length per bin for each character:
			for(j in MissingStringCharacters[[i]]) {
			
				# Store edge length per bin:
				edge.lengths.per.bin[j, ] <- edge.length.per.bin$edge.length.in.bin
				
				# Store terminal edge length per bin:
				terminal.edge.lengths.per.bin[j, ] <- edge.length.per.bin$terminal.edge.length.in.bin
				
				# Store internal edge length per bin:
				internal.edge.lengths.per.bin[j, ] <- edge.length.per.bin$internal.edge.length.in.bin
			
			}
			
		}
		
	}
	
	# Add column names to edge lengths per bin:
	colnames(terminal.edge.lengths.per.bin) <- colnames(internal.edge.lengths.per.bin) <- colnames(edge.lengths.per.bin) <- names(edge.length.per.bin$edge.length.in.bin)
	
	# Get rid of pesky state rownames:
	rownames(allstatechanges) <- NULL
	
	# Sort all state changes by edge number:
	allstatechanges <- allstatechanges[order(allstatechanges[, "Edge"]), ]
	
	# Sort all state changes by character number:
	allstatechanges <- allstatechanges[order(allstatechanges[, "Character"]), ]

	# Sort all state changes by simulation number:
	allstatechanges <- allstatechanges[order(allstatechanges[, "Nsim"]), ]
	
	# Compile output as list:
	output <- list(allstatechanges, chartime, edge.lengths.per.bin, terminal.edge.lengths.per.bin, internal.edge.lengths.per.bin)
	
	# Add names to output:
	names(output) <- c("all.state.changes", "character.times", "edge.length.per.bin", "terminal.edge.length.per.bin", "internal.edge.length.per.bin")
	
	# Return output:
	return(output)
	
}
