# function that creates a transition layer, corrects diagnols, and creates a cost distance matrix
distanceMatrix <- function(x, transitionFunction, directions, symm, fromCoords, toCoords, dist.cutoff, cluster)
{
  ### modified transition function
  # calculate conductances among neighbors
  tr <- new("TransitionLayer",
            nrows=as.integer(nrow(x)),
            ncols=as.integer(ncol(x)),
            extent=extent(x),
            crs=projection(x, asText=FALSE),
            transitionMatrix = Matrix(0,ncell(x),ncell(x)),
            transitionCells = 1:ncell(x))
  transitionMatr <- transitionMatrix(tr)
  Cells <- which(!is.na(getValues(x)))
  adj <- adjacent(x, cells=Cells, pairs=TRUE, target=Cells, directions=directions)
  if(symm){adj <- adj[adj[,1] < adj[,2],]}
  dataVals <- cbind(getValues(x)[adj[,1]],getValues(x)[adj[,2]])
  # run in parallel using parapply
  transition.values <- parApply(cluster,dataVals,1,function(x) 1/mean(x))
  if(!all(transition.values>=0)){warning("transition function gives negative values")}
  transitionMatr[adj] <- as.vector(transition.values)
  if(symm)
  {
    transitionMatr <- forceSymmetric(transitionMatr)
  }
  transitionMatrix(tr) <- transitionMatr
  matrixValues(tr) <- "conductance"
  ### geocorrect the transition matrix
  # adjust diag. conductances
  trCorr <- geoCorrection(tr, type="c", multpl = FALSE, scl = FALSE)
  ### modified costdistance function
  ## create a list of sites only within dispersal distance (Euclidian) from each site
  toCoords.loc <- vector("list", nrow(toCoords)) # location of each retained coordinate in the distance matrix
  toCoords.list <- vector("list", nrow(toCoords)) # list of coordinates within dispersal distance of each individual point (list in order of points)
  for(j in 1:nrow(toCoords)){
    toCoords.loc[[j]] <- which(sqrt((fromCoords[j,1]-fromCoords[,1])^2 + 
                                      (fromCoords[j,2]-fromCoords[,2])^2) < dist.cutoff)
    toCoords.list[[j]] <- toCoords[toCoords.loc[[j]],]
  }
  # indicates which cells the fromCoords are in
  # remove duplicates with unique function
  fromCells <- cellFromXY(trCorr, fromCoords)
  
  # have to set up toCells a bit different
  toCells <- toCells_gdist <- vector("list", nrow(fromCoords)) # list to hold the toCells for each individual fromCell
  # indicate the toCells, but keep them within the list to hold true to the dispersal cutoff
  for(i in 1:length(toCoords.list)){
    toCells[[i]] <- cellFromXY(trCorr, toCoords.list[[i]])
    toCells_gdist[[i]] <- unique(toCells[[i]])
  }
  
  # make a list object that can be run through parLapply
  cell_list <- vector("list", length(fromCells))
  for(i in 1:length(cell_list)){
    cell_list[[i]] <- list(fromCells[i],toCells_gdist[[i]])
  }
  
  # make the transition matrix an object that igraph can read
  y <- gdistance::transitionMatrix(trCorr)
  # create an adjacencyGraph to get edge weights
  if(isSymmetric(y)) {m <- "undirected"} else{m <- "directed"}
  adjacencyGraph <- graph.adjacency(y, mode=m, weighted=TRUE)
  # reclassify edge weights as cost
  E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
  
  # create wrapper for shortest.paths function to go into lApply
  shortestPaths <- function(x){ shortest.paths(adjacencyGraph,
                                               v = x[[1]],
                                               to = x[[2]],
                                               mode = "out",
                                               algorithm = "dijkstra") }
  
  # export function and variables to clusters
  clusterExport(cl=cluster, list("adjacencyGraph", "shortestPaths", 
                                 "shortest.paths"),
                envir=environment())
  # use shortestPaths function in parallel
  costDist <- parLapply(cl=cluster,cell_list, shortestPaths)
  # remove exported variables to make room for next run
  clusterEvalQ(cluster, rm(list=ls()))
  # turn distances into full distance matrix
  D_mat <- matrix(NA, nrow(fromCoords), nrow(fromCoords))
  for(i in 1:nrow(fromCoords)){
    # create a vector so that removed cells are giving same value as its matching cell
    D_mat[toCoords.loc[[i]],i] <- left_join(x=data.frame(to=toCells[[i]]),
                                            y=data.frame(to=toCells_gdist[[i]],
                                                         val=costDist[[i]][1,]), by="to")$val
  }
  rownames(D_mat) <- rownames(fromCoords)
  colnames(D_mat) <- rownames(fromCoords)
  
  return(D_mat)
}
