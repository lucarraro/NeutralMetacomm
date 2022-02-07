eval_evenness <- function(speciesDB){
  
  nNodes <- length(unique(speciesDB[,1]))
  evenness <- numeric(nNodes)
  for (indNode in 1:nNodes){
    subset <- which(speciesDB[,1]==indNode)
    speciesArray <- speciesDB[subset,2] 
    values <- unique(sort(speciesArray))
    counts <- numeric(length(values))
    for (ind in 1:length(values)){
      counts[ind] <- sum(speciesArray==values[ind])
    }
    p <- counts/sum(counts)
    H <- - sum(p*log(p)) # Shannon entropy
    evenness[indNode] <- H/log(length(values))
  }
  return(evenness)
}