#' ---
#' title: "Influentials, networks and public opinion formation"
#' author: Watts and Dodds 2007
#' output: github_document
#' ---

library(tidyverse)
library(igraph)
library(gridExtra)

options(scipen=999)  # turn-off scientific notation like 1e+48
theme_set(theme_bw())
set.seed(20130810)

dev.new()

#' *Note: This script is my effort to replicate the results from the paper
#' "Influentials, networks and public opinion formation", Watts and Dodds 2007.
#' This is a self-didactic attempt.*
#'
#' In marketing, there is a hugely popular idea that a small group of
#' "influencers" are able to drive large scale diffusion of products/information
#' on a scale beyond what can be achieved by ordinary individuals. In this
#' script, we take a look at this hypothesis and explore the boundary conditions
#' in which this idea is true.
#'
#' Before we begin our exploration, it is instructive to look at the degree
#' distribution of the underlying network on which we simulate the diffusion
#' process. At the outset, note that we focus only on networks that have a
#' balanced distribution around a mean for node degrees in this script. This
#' means that most nodes have the same degree (equal to the mean), and there are
#' a few nodes which have degrees spread about the mean. It is an open question
#' whether such a distribution of degrees is reflective of the real-world and is
#' not a subject of exploration here.

numNodes <- 1e4
degreeDistributions <- bind_rows(data_frame(d = degree(sample_smallworld(dim = 1, 
                                                                         size = numNodes, 
                                                                         nei = 2, 
                                                                         p = 0.5)), 
                                            m = rep(paste("Average degree: ", mean(d)), numNodes)),
                                 data_frame(d = degree(sample_smallworld(dim = 1, 
                                                                         size = numNodes, 
                                                                         nei = 6, 
                                                                         p = 0.5)), 
                                            m = rep(paste("Average degree: ", mean(d)), numNodes)),
                                 data_frame(d = degree(sample_smallworld(dim = 1, 
                                                                         size = numNodes, 
                                                                         nei = 12, 
                                                                         p = 0.5)), 
                                            m = rep(paste("Average degree: ", mean(d)), numNodes)),
                                 data_frame(d = degree(sample_smallworld(dim = 1, 
                                                                         size = numNodes, 
                                                                         nei = 16, 
                                                                         p = 0.5)), 
                                            m = rep(paste("Average degree: ", mean(d)), numNodes)))

ggplot(degreeDistributions) +
  geom_histogram(aes(x = d)) +
  labs(x = "Node degree", y = "Count", 
       title = "Histograms of degre distribution of sample small world graphs \n(10,000 nodes and varying mean degrees)") + 
  facet_wrap(~m, ncol = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#' We model the impact of seeding to specific individuals by running diffusion
#' simulations in two scenarios. In the first scenario, we seed the simulation
#' by activating a randomly chosen node from the network. In the second
#' scenario, we seed the most connected node (i.e., highest degree) in the
#' network.
#'
#' We make the following assumptions on the decision making behavior of the
#' nodes in a network:
#'
#' - A node in a graph is faced with a binary decision - to engage or to not
#' engage (e.g., with new products or discussions with friends)
#'
#' - They make this decision based on a simple rule: they compute the fraction
#' of their neighbors that have engaged, compare it with their personal
#' threshold, and engage if the fraction of engaged neighbors exceeds the
#' threshold. In other words, the thought process is on the lines, "If at least
#' 18% of my friends bought the new iPhone, I would want to buy it too".
#'
#' - At each update, we compute the vulnerable nodes by comparing the fraction
#' of their engaged neighbors to the threshold, compare this to the threshold
#' and activate the node if the fraction is greater than the threshold.
#'
#' # Model
#'
#' Create a new small world graph at each realization. This represents a
#' distribution of node degrees and is supposed to be reflective of the
#' influence distribution among the nodes in the network. Two kinds of seeding
#' is conducted - seeding the most influential person in the network (e.g., one
#' with maximum degree) and a randomly chosen person in the network. The network
#' then updates till no new node activations are possible at this configuration
#' of the network and the threshold. The simulation is executed a large number
#' of times and the cumulative number of activations is counted for each case.

#' *Hyper Parameters of the model*
#'
#' 1. Number of nodes in the Barabasi-Albert graph (n)
#'
#' 2. Average degree (z)
#'
#' 3. Threshold (distribution or a specific value)
#'
#' 4. Number of realizations
#'
#' *Output* A matrix of number of activations in random and influential seeding
#' and time steps to stability for each case at the end of each realization of
#' the simulation

diffusionSimulation <- function(n, z, threshold, numRealizations) {
  output <- matrix(nrow = numRealizations, ncol = 4)
  colnames(output) <- c("RandomActivations", "RandomTimesteps", "InfluentialActivations", "InfluentialTimesteps")
  k <-  floor(z/2)
  
  for (r in 1:numRealizations) {
    g <-  sample_smallworld(dim = 1, size = numNodes, nei = k, p = 0.5)
    A <- get.adjacency(g)
    degrees <- degree(g)
    thresholdNum <- threshold * degrees
    
    # Random Seeding begins here
    nodeStatus <- rep(0, n)
    seed <- sample(V(g), size = 1)
    nodeStatus[seed] <- 1
    
    t <- 1
    oldNodeStatus <- nodeStatus
    numEngagedNeighbors <- as.vector(A%*%nodeStatus)
    vulnerableNodes <- numEngagedNeighbors > thresholdNum
    nodeStatus[vulnerableNodes] <- 1
    
    while (max(nodeStatus != oldNodeStatus) > 0) {
      oldNodeStatus <- nodeStatus
      numEngagedNeighbors <- as.vector(A%*%nodeStatus)
      vulnerableNodes <- numEngagedNeighbors > thresholdNum
      nodeStatus[vulnerableNodes] <- 1
      t <- t + 1
    } 
    
    randomActivations <- sum(nodeStatus)
    randomTimesteps <- t
    
    # Seeding the most connected member begins here
    nodeStatus <- rep(0, n)
    seed <- which(degree(g) == max(centralization.degree(g)$res)) # Who has the maximum degree?
    nodeStatus[seed] <- 1
    
    t <- 1
    oldNodeStatus <- nodeStatus
    numEngagedNeighbors <- as.vector(A%*%nodeStatus)
    vulnerableNodes <- numEngagedNeighbors > thresholdNum
    nodeStatus[vulnerableNodes] <- 1
    
    while (max(nodeStatus != oldNodeStatus) > 0) {
      oldNodeStatus <- nodeStatus
      numEngagedNeighbors <- as.vector(A%*%nodeStatus)
      vulnerableNodes <- numEngagedNeighbors > thresholdNum
      nodeStatus[vulnerableNodes] <- 1
      t <- t + 1
    } 
    
    
    influentialActivations <- sum(nodeStatus)
    influentialTimesteps <- t
    
    output[r, ] <- c(randomActivations, randomTimesteps, influentialActivations, influentialTimesteps)
  }
  
  return (output)
}

#' Minimal working examples and analysis
#'
#' We illustrate our script by running it for a set of average degrees and
#' thresholds. We assume that all the nodes in the network have the same
#' fractional threshold for simplicity. We run the simulation at four points
#' that mimic a factorial design. We vary (z, threshold) as (4, 0.01), (4, 0.2),
#' (32, 0.01), (32, 0.4).
#' 

numNodes <- 1e4
system.time(data1 <- diffusionSimulation(n = numNodes, z = 4, rep(0.01, numNodes), numRealizations = 100))
system.time(data2 <- diffusionSimulation(n = numNodes, z = 4, rep(0.2, numNodes), numRealizations = 100))
system.time(data3 <- diffusionSimulation(n = numNodes, z = 32, rep(0.01, numNodes), numRealizations = 100))
system.time(data4 <- diffusionSimulation(n = numNodes, z = 32, rep(0.2, numNodes), numRealizations = 100))

results <- bind_rows(as.data.frame(data1), 
                     as.data.frame(data2),
                     as.data.frame(data3),
                     as.data.frame(data4))

glimpse(results)
resultSummary <- results %>% mutate(ElevationInActivations = InfluentialActivations/RandomActivations,
                                    RandGlobalCascade = map_lgl(results$RandomActivations, function(x) {return(x>1000)}),
                                    InfluentialGlobalCascade = map_lgl(results$InfluentialActivations, function(x) {return(x>1000)})) %>% 
                             summarize(MeanElevation = mean(ElevationInActivations),
                                       MedianElevation = median(ElevationInActivations),
                                       NumRandGlobalCascades = sum(RandGlobalCascade),
                                       NumInfluentialGlobalCascades = sum(InfluentialGlobalCascade))
print(resultSummary)

#' *Final Analysis*
#'
#' If we look at the mean elevation across all the simulation runs on the sample
#' space, it looks like influentials provide a big jump in elevation compared to
#' average individuals. This is a big `+1` for the influential hypothesis. A
#' direct implication of this finding is that most marketing dollars should be
#' directed towards identifying and 'influencing' the influencers.

# No wonder 'influencer marketing' was the buzz word for several years!
#   
# This is not the full picture though. Averages are severely misleading. 
# A key idea of the Tipping Point is that these 'special' individuals are able 
# to drive large swathes of ordinary individuals into adopting the idea because 
# of their endorsement. Such large scale sequence of adoptions following a single 
# adoption are called cascades. To be really effective, influentials should be 
# able to drive cascades that are large not only on average, but also in scale. 
# They have to be really, really, really big. The kind of big that poor ordinary 
# individuals cannot generate, ever. 

# These kind of super large cascades are called global cascades.
# 
# A reasonable definition of a global cascade is one that occupies at least 10% 
# of the network (this is based on prior literaure, but of course it is up to you 
# to pick a cut-off). 
# For the 10,000 node network we consider in this notebook, global cascades are 
# those that reach 1000 nodes. 

# On this metric, we see that influentials drive roughly 25% more global cascades 
# compared to randomly chosen individuals. 

# Is this result worth the effort?