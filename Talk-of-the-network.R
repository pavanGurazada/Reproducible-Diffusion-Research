#' ---
#' title: "Talk of the network, the tidy way "
#' author: Pavan Gurazada
#' output: github_document
#' ---

#' last update: Sun Mar 04 10:04:24 2018

#' An attempt to replicate the results from the influential paper, "Talk of the
#' network", using the igraph-tidygraph combination. To begin with, we explore
#' the tidygraph interface and see if the logic of the Julia version can be
#' translated. TO begin with we make little use of `tidygraph` and build the
#' first version on `igraph` and base R. Next step would ne to translate parts
#' of the code to Rcpp, and finally a different version using tidygraph. This
#' will enable comparisons between the approaches.
#'
#' A small example, to begin with:

library(tidygraph)
library(tidyverse)
library(igraph)

set.seed(20130810)

play_erdos_renyi(10, 0.5) %>% 
  activate(nodes) %>% 
  mutate(degree = centrality_degree()) %>% 
  activate(edges) %>% 
  mutate(centrality = centrality_edge_betweenness()) %>% 
  arrange(centrality)

#' As we can see from these results, the edge list serves as the core graph
#' structure by itself. tidygraph is an interesting wrapper around igraph
#' 
#' Explore this further using the examples from the manual

gr <- create_complete(5) %>% 
        activate(nodes) %>% 
        mutate(class = sample(c('a', 'b'), 5, replace = TRUE)) %>% 
        activate(edges) %>% 
        arrange(from)

#' Two tidy data frames keep track of the network data and the meta data
#' associated with the network. The `nodes` data frame can keep track of the
#' meta data, i.e., degree, etc, while the `edges` data frame can keep track of
#' the connections and the edge related meta data.
#'
#' The built-in graph creators are divided into two types - first, where there
#' is some randomization involved in building the graphs (`play_...`) and
#' second, where the graph creation process is deterministic (`create_...`).
#' 

create_star(10, directed = TRUE, mutual = TRUE) %>% 
  activate(edges) %>% 
  sample_frac(0.7) %>% 
  mutate(single_edge = !edge_is_mutual())

create_notable('chvatal') %>% 
  activate(nodes) %>% 
  mutate(neighborhood = local_members(mindist = 1))

#' Now, moving to the paper 'Talk of the network'
#'
#' The central logic of the diffusion simulation from the paper is as follows:
#' "Each individual belongs to a single personal network. Each network consists
#' of individuals who are connected by strong ties. In each period, individuals
#' also conduct a finte number of weak tie interactions outside their personal
#' networks... We divide the entire market equally into personal networks, in
#' which each individual can belong to one network. In addition, in each period,
#' every individual conducts random meetings with individuals external to his
#' personal network."
#'
#' So, to begin with we need several complete graphs, and then the individuals
#' within these subgraphs interact randomly. Our final data structure is hence 
#' a vector of several complete networks that are built based on the number of 
#' strong ties for each individual. Note that each individual in the network has 
#' a fixed number of strong ties (`s`) and weak ties (`w`).

initialize_graph <- function(n, n_strong_ties) {
  n_subgraphs <- floor(n/n_strong_ties)
  
  G <- lapply(1:n_subgraphs, function(g) make_full_graph(n_strong_ties))
 
  return(G)
}

#' The node status is stored as a vector of logical vectors. At the beginning of 
#' each simulation run, we call the following function to set the status of all 
#' the nodes to `FALSE`.

reset_node_status <- function(G) {
  
  node_status <- lapply(1:length(G), function(g) rep(FALSE, vcount(G[[g]])))

  return(node_status)
}

#' At each time step, we execute two tasks. First, we allow the nodes to mingle 
#' randomly with their strong ties and with weak ties from other sub-networks. 
#' At this point, we count the number of active strong and weak ties for each 
#' node. Then, we use this information to update the status of all the nodes in 
#' the network.
#' 
#' The first function counts the number of active strong ties within the node's 
#' sub-network. The second function executes the "random meetings" with weak 
#' ties as discussed in the paper. For each node we generate a random sample 
#' (without replacement) of size `w` from sub-networks other than its own. We 
#' then count the number of active ties in its own sub-network and among the 
#' random sample taken from the rest of the network.

count_active_str_ties <- function(G, node_network_id, node, node_status) {
  
  nbrs <- neighbors(G[[node_network_id]], node)
  
  nbr_status <- node_status[[node_network_id]][nbrs]
  
  return(sum(nbr_status))
}

# The following function is the one that should be an ideal candidate to be
# written in C++

random_meetings <- function(G, node_network_id, node, node_status, n_weak_ties) {
  all_network_ids <- 1:length(G)
  
  # Figure out the ids of networks other than the nodes original network
  
  other_network_ids <- all_network_ids[all_network_ids != node_network_id]
  
  possible_weak_ties <- data.frame(Network = integer(), Weak_Tie = integer())
  nsamples <- 1
  
  while(nsamples < n_weak_ties) {
    rand_network_id <- sample(other_network_ids, size = 1)
    rand_nbr <- sample(as.vector(V(G[[rand_network_id]])), size = 1)
    possible_weak_tie <- data.frame(Network = rand_network_id, Weak_Tie = rand_nbr)
    
    if (nrow(anti_join(possible_weak_ties, possible_weak_tie, by = c("Network", "Weak_Tie"))) == 0) {
      possible_weak_ties <- rbind(possible_weak_ties, possible_weak_tie)
      nsamples <- nsamples + 1
    }
  }
  
  return(possible_weak_ties)
}

count_active_wk_ties <- function(possible_weak_ties, node_status) {
  n_active_wk_ties <- 0
  
  for (i in 1:nrow(possible_weak_ties)) {
    if (node_status[[possible_weak_ties[i, 1]]][possible_weak_ties[i, 2]])
      n_active_wk_ties <- n_active_wk_ties + 1
  }
  
  return(n_active_wk_ties)
}

#' One more ideal candidate for C++

evolve <- function(G, node_status, n_weak_ties, alpha, beta_w, beta_s) {
  for (node_network_id in sample(1:length(G))) {
    for (node in as.vector(sample(V(G[[node_network_id]])))) {
      if (!node_status[[node_network_id]][node]) {
        n_active_str_ties <- count_active_str_ties(G, node_network_id, node, node_status)
        possible_weak_ties <- random_meetings(G, node_network_id, node, node_status, n_weak_ties)
        n_active_wk_ties <- count_active_wk_ties(possible_weak_ties, node_status)
        
        activation_prob <- 1 - (1 - alpha) * (1 - beta_w)^n_active_wk_ties * (1 - beta_s)^n_active_str_ties
        
        if (runif(n = 1) < activation_prob) node_status[[node_network_id]][node] = TRUE
      }
      
    }
  }
  return(node_status)
}

#' One more candidate for C++ is the function below

simulate <- function(parameter_space, n_nodes) {
  
  output <- data.frame(s = integer(), w = integer(), alpha = double(), 
                       beta_w = double(), beta_s = double(), 
                       t = integer(), num_engaged = integer())
  
  cat("\n Beginning simulation at: ", date())
  
  for (i in 1:nrow(parameter_space)) {
    s_i <- parameter_space$s[i]
    w_i <- parameter_space$w[i]
    alpha_i <- parameter_space$alpha[i]
    beta_w_i <- parameter_space$beta_w[i]
    beta_s_i <- parameter_space$beta_w[i]
    
    G <- initialize_graph(n_nodes, s_i)
    node_status <- reset_node_status(G)
    
    n_engaged <- reduce(node_status, sum)
    t_i <- 1
    
    while (n_engaged < floor(0.95 * n_nodes)) {
      node_status <- evolve(G, node_status, w_i, alpha_i, beta_w_i, beta_s_i)
      n_engaged <- reduce(node_status, sum)
      output <- rbind(output, data.frame(s = s_i, w = w_i, alpha = alpha_i, 
                                         beta_w = beta_w_i, beta_s = beta_s_i,
                                         t = t_i, num_engaged = n_engaged))
      t_i <- t_i + 1
    }
    cat("\n Finished setting: ", i, " at : ", date())
  }
  
  return (output)
}

parameter_space <- expand.grid(s = seq(5, 29, length.out = 3), 
                               w = seq(5, 29, length.out = 3), 
                               alpha = seq(0.0005, 0.01, length.out = 3), 
                               beta_w = seq(0.005, 0.015, length.out = 3), 
                               beta_s = seq(0.01, 0.07, length.out = 3))

results <- simulate(parameter_space, 3000)

