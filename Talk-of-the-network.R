#' ---
#' title: "Talk of the network, the tidy way "
#' author: Pavan Gurazada
#' last update: Sat Mar 03 06:54:34 2018
#' output: github_document
#' ---

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
  
  G <- vector("list", n_subgraphs)
  
  for (i in 1:n_subgraphs) {
    G[[i]] <- make_full_graph(n_strong_ties)
  }
  
  return(G)
}

#' The node status is stored as a vector of logical vectors. At the beginning of 
#' each simulation run, we call the following function to set the status of all 
#' the nodes to `FALSE`.

reset_node_status <- function(G) {
  node_status <- vector("list", length(G))
  
  for (g in 1:length(G))
    node_status[[g]] <- rep(FALSE, vcount(G[[g]]))
  
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
  
}

