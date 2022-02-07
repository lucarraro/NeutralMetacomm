rm(list=ls())
library(OCNet) # Generate OCNs
library(betapart) # Calculate beta diversity  
library(usedist) # Treat "dist" objects

setwd("C:/Users/carrarlu/Documents/Git/NeutralMetacomm") # change to current directory
source("metacomm_model.R") # Function running the neutral metacommunity model
source("eval_evenness.R") # Function to compute Pielou's evenness

## ALTERNATIVE IF LOADING OCNET DOESN'T WORK - uncomment the following 3 lines
# library(spam)
# library(fields)
# source("draw_thematic_OCN.R")  

# LOAD OCN ####
if (!file.exists("OCN.rda")){
  set.seed(2); 
  OCN <- create_OCN(250, 250, cellsize=100, typeInitialState = "V") # this command would generate the OCN
  # this is a 250x250 OCN, with cellsize = 100 (m). The catchment has a total drainage area of 625 km2
  OCN <- landscape_OCN(OCN, zMin=0, slope0=5e-3) # evaluate elevations and slopes 
  # (this is necessary for the application of the following functions, although elevations are not used)
  OCN <- aggregate_OCN(OCN, maxReachLength = 2500, thrA=1e6) # aggregate the OCN such that each reach has max length = 2.5 km; Threshold area = 1 km2
  OCN <- paths_OCN(OCN, includeUnconnectedPaths = TRUE)
  save(OCN,file="OCN.rda")
} else {load("OCN.rda")}

# some quantities that can be useful to analyze the resulting patterns
drainageArea <- OCN$AG$A # Drainage area vector
nNodes <- OCN$AG$nNodes # number of nodes (patches)
distanceMatrix <- OCN$AG$downstreamPathLength + OCN$AG$downstreamLengthUnconnected + 
  (t(OCN$AG$downstreamPathLength) + t(OCN$AG$downstreamLengthUnconnected)) # matrix of pairwise distances between nodes
                                                                           # It contains both downstream and upstream (unweighted) distances
dist_mat <- as.dist(distanceMatrix) # transform into "dist" object
distanceToOutlet <- OCN$AG$downstreamPathLength[,OCN$AG$outlet] # vector of distances to the outlet
wClosenessCentrality <- rowSums(exp(-distanceMatrix/2000)) # weighted closeness centrality ("how many patches are available close to a given node?")
                                                           # as the scale length (2000) increases, this variable passes from local characteristics (degree) 
                                                           # to global (total number of connected nodes) 

# METACOMMUNITY MODEL ####
# Define parameters
nu <- 1e-3 # Diversification rate
habitatIndex <- 100+numeric(nNodes) # Distribution of subhabitats across reaches
                                    # FLAT DISTRIBUTION: 100+numeric(nNodes) // SCALING WITH DRAINAGE AREA: round(drainageArea/310400)
medianDistanceExp <- 5000 # Median distance travelled according to the exponential distribution [in meters]
medianDistanceCauchy <- 0 # Same for Cauchy distribution
                          # The two above parameters define the dispersal kernel. If one is set to 0, then only the other part of the distribution is used
weightUpstream <- 5 # Weight for distances travelled upstream
speciesCount <- 1000 # Initial number of species
nTimesteps <- 2e5 # Total number of timesteps
save_alphaDiv_time <- FALSE # If TRUE, store matrix alphaDiv_time and vector gammaDiv_time showing time evolution of alpha & gamma diversity 
                            # (useful to check for convergence, but slows down the algorithm) 

# Wrap parameters
parameters <- list(nu=nu, habitatIndex=habitatIndex, medianDistanceExp=medianDistanceExp,
                   medianDistanceCauchy=medianDistanceCauchy, weightUpstream=weightUpstream,
                   speciesCount=speciesCount, nTimesteps=nTimesteps, save_alphaDiv_time=save_alphaDiv_time)

set.seed(1) # set seed for reproducibility
out <- metacomm_model(OCN, parameters) # Run model
for(i in 1:length(out)) {assign(names(out)[i], out[[i]])} # Unwrap results

# CALCULATE DIVERSITY INDICES ####
# build matrix of presence/absence
PA_mat <- matrix(0,nNodes,speciesCount)
for (i in 1:sum(habitatIndex)){
  PA_mat[speciesDB[i,1],speciesDB[i,2]] <- 1
}
PA_mat <- PA_mat[,colSums(PA_mat)!=0] # remove extinct species

alphaDiv <- rowSums(PA_mat) # alpha diversity
betaDiv_mat <- beta.pair(PA_mat, index.family="jaccard") # pairwise beta diversity matrix (Jaccard distance)
evenness <- eval_evenness(speciesDB) # Pielou's evenness

# alpha diversity map
x11(); draw_thematic_OCN(alphaDiv,OCN); title("alpha diversity (species richness)")

# evenness map
x11(); draw_thematic_OCN(evenness,OCN); title("Pielou's evenness")

# compare beta diversity on headwaters vs other streams (based on median drainage area)
p1 <- hist(dist_subset(betaDiv_mat$beta.jac, OCN$AG$A<=median(OCN$AG$A)),breaks=seq(0,1,0.02),plot=FALSE) # pairwise beta div for upstream reaches
p2 <- hist(dist_subset(betaDiv_mat$beta.jac, OCN$AG$A>median(OCN$AG$A)),breaks=seq(0,1,0.02),plot=FALSE)  # pairwise beta div for downstream reaches

x11(); plot(p1, col=rgb(0,0,1,1/4), xlim=c(0.5,1), 
      xlab="Jaccard distance", main="Beta diversity upstream/downstream") 
plot(p2, col=rgb(1,0,0,1/4), add=T)
legend("topleft",legend=c("upstream","downstream"),fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
abline(v=mean(dist_subset(betaDiv_mat$beta.jac, OCN$AG$A<=median(OCN$AG$A))),col="blue",lwd=1.5)
abline(v=mean(dist_subset(betaDiv_mat$beta.jac, OCN$AG$A>median(OCN$AG$A))),col="red",lwd=1.5)


# variation of beta diversity with distance
x11(); par(mfrow=c(1,3)) 
plot(dist_mat, betaDiv_mat$beta.jac,ylab="Jaccard distance",xlab="Along-stream distance [m]") # Total Jaccard distance
plot(dist_mat, betaDiv_mat$beta.jne,ylab="Nestedness",xlab="Along-stream distance [m]"); title("Beta diversity vs. distance") # Nestedness component of Jaccard distance
plot(dist_mat, betaDiv_mat$beta.jtu,ylab="Turnover",xlab="Along-stream distance [m]") # Turnover component of Jaccard distance


# rank abundance curves
x11(); plot(sort(colSums(PA_mat),decreasing=TRUE)/nNodes,type="l",
     ylab="Relative abundance",xlab="Rank",log="y",ylim=c(1e-5,1)) # rank abundance (based on presence at a site)
title("Rank abundance curves")
hh <- hist(speciesDB[,2], breaks=0.5:(speciesCount+0.5),plot=F) # rank abunance (based on true abundance)
points(sort(hh$counts[hh$counts>0],decreasing=TRUE)/sum(habitatIndex),col="red",type="l")
legend("topright",legend=c("Based on presence at-a-site","Based on true abundance"),lty=c(1,1),col=c("black","red"))

