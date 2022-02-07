metacomm_model <- function(OCN, parameters){

nNodes <- OCN$AG$nNodes
for(i in 1:length(parameters)) {assign(names(parameters)[i], parameters[[i]])} # unwrap parameters

distanceMatrix <- OCN$AG$downstreamPathLength + OCN$AG$downstreamLengthUnconnected + 
  weightUpstream*(t(OCN$AG$downstreamPathLength) + t(OCN$AG$downstreamLengthUnconnected)) # weighted distance matrix

speciesDB <- matrix(0,sum(habitatIndex),2) # initialize species database
                                           # col 1: ID of reach to which the subhabitat (row) is attributed
                                           # col 2: ID of species occupying the subhabitat
k <- 1
for (i in 1:nNodes){
  speciesDB[k:(k+habitatIndex[i]-1),1] <- i
  k <- k + habitatIndex[i]
}
speciesDB[,2] <- sample(speciesCount, sum(habitatIndex), replace=TRUE)

# Define dispersal kernel (both defined in terms of the medians of their distribution)
dispersalExp <- exp(-distanceMatrix*log(2)/medianDistanceExp)
dispersalCauchy <-  medianDistanceCauchy^2/(medianDistanceCauchy^2 + distanceMatrix^2)
dispersalExp[is.nan(dispersalExp)] <- 0       # patch for when medianDistanceExp = 0
dispersalCauchy[is.nan(dispersalCauchy)] <- 0 # patch for when medianDistanceCauchy = 0

dispersal <- dispersalExp + dispersalCauchy

dispersalKernel <- habitatIndex*dispersal # dispsersalKernel depends on habitat size and distance from extinction node
dispersalKernel  <- t(t(dispersalKernel)/colSums(dispersalKernel)) # normalize dispersalKernel (= probability to pick a given node)
                                                                   # equivalent to dividing each column of the matrix by its sum

cumProbMatrix <- rbind(numeric(nNodes),apply(dispersalKernel,2,cumsum)) # each column corresponds to a node; 
                                                                        # values are cumulative probabilities to pick a given node as destination 

if (save_alphaDiv_time){
  alphaDiv_time <- matrix(0,nTimesteps/100,nNodes) # initialize alphaDiv_time
  gammaDiv_time <- numeric(nTimesteps/100) # initialize gammaDiv_time
  } 

t0 <- Sys.time()
for (iter in 1:nTimesteps){
  unitExtinct <- sample(sum(habitatIndex),1)  # random extinction
  if (runif(1) < nu){ # speciation event
    speciesCount <- speciesCount + 1
    speciesDB[unitExtinct,2] <- speciesCount
  } else { # replacement from existing species
    nodeEmigration <- which(cumProbMatrix[,speciesDB[unitExtinct,1]] > runif(1))[1] - 1      # pick random node from which the immigrant arrives
    speciesDB[unitExtinct,2] <-  speciesDB[sample(which(speciesDB[,1]==nodeEmigration),1),2] # pick random subhabitat from that node
  }
  if ((iter %% 100)==0){ # write progress (and possibly evaluate alphaDiv_time)
    cat(sprintf("%.1f%% done  - Elapsed time: %.2f s\r",iter/nTimesteps*100,difftime(Sys.time(),t0,units="secs")))
    if (save_alphaDiv_time){
      PA_mat <- matrix(0,nNodes,speciesCount) # create presence/absence matrix
      for (i in 1:sum(habitatIndex)){
        PA_mat[speciesDB[i,1],speciesDB[i,2]] <- 1
      }
      PA_mat <- PA_mat[,colSums(PA_mat)!=0] # remove extinct species
      alphaDiv_time[iter/100,] <- rowSums(PA_mat)
      gammaDiv_time[iter/100] <- length(unique(speciesDB[,2]))
    }
  }
}
ll <- list(speciesDB=speciesDB, speciesCount=speciesCount) # store species database and total no. species generated
if (save_alphaDiv_time){ll[["alphaDiv_time"]] <- alphaDiv_time}
if (save_alphaDiv_time){ll[["gammaDiv_time"]] <- gammaDiv_time}

return(ll)
}