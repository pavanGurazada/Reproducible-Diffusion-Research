#' ---
#' title: "Talk of the network, the tidy way "
#' author: Pavan Gurazada
#' date: Sat Mar 03 06:54:34 2018
#' output: github_document
#' ---

#' An attempt to replicate the results from the influential paper, "Talk of the 
#' network", using the igraph-tidygraph combination. 
#' 
#' A small example, to begin with:

library(tidygraph)

play_erdos_renyi(10, 0.5) %>% 
  activate(nodes) %>% 
  mutate(degree = centrality_degree()) %>% 
  activate(edges) %>% 
  mutate(centrality = centrality_edge_betweenness()) %>% 
  arrange(centrality)

#' As we can see from these results, the edge list serves as the core graph
#' structure by itself. tidygraph is an interesting wrapper around igraph

