#include <Rcpp.h>

using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::IntegerVector;

using std::vector;

// [[Rcpp::plugins("cpp11")]]

//[[Rcpp::export]]
vector<List> initialize(const List& g, const int& n_subgraphs) {
  vector<List> G;
  
  for (int i = 0; i < n_subgraphs; ++i) {
    G.push_back(g);
  }
  
  return G;
}

//[[Rcpp::export]]
vector<LogicalVector> reset(const vector<List>& G) {
  vector<LogicalVector> node_status;
  
  for (auto& g : G) {
    LogicalVector ns;
    
    for (int i = 0; i < g.length(); ++i) {
      ns.push_back(false);
    }
    
    node_status.push_back(ns);
  }
  
  return node_status;
}

//[[Rcpp::export]]
int count_active_str_ties(const vector<List>& G, 
                          const int& node_network_id, 
                          const int& node, 
                          const vector<LogicalVector>& node_status) {
  int n_active_str_ties = 0;
  
  IntegerVector nbrs = G[node_network_id][node];
  
  for (auto& nbr : nbrs) {
    if (node_status[node_network_id][nbr] == true)
      n_active_str_ties++;
  }
  
  return n_active_str_ties;
}

int random_meetings(const vector<List>& G, 
                    const int& node_network_id, 
                    const int& node, 
                    const vector<LogicalVector>& node_status, 
                    const int& n_weak_ties) {
  int n_active_wk_ties;
  
  IntegerVector all_network_ids;
  
  auto it = G.begin();
  while (it != G.end()) {
    
    
    ++it;
  }
  
  return n_active_wk_ties;
}


/*** R
library(igraph)

g <- make_full_graph(5, loops = FALSE)

G <- initialize(neighborhood(g), 2)

ns <- reset(G)

ns
*/
