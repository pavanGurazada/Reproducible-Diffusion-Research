#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

using std::pow;
using std::vector;
using std::random_shuffle;
using std::normal_distribution;
using std::count;

vector<int> neighbours(const sp_mat& A, const int& node) {
  vector<int> nbrs; 
  int n_nodes = A.n_cols;
  
  for (int i = 0; i < n_nodes; ++i) {
    if (A(i, node-1) == 1) nbrs.push_back(i+1);
  }
  
  return nbrs;
}

vector<int> reset_nodes(const sp_mat& A) {
  vector<int> node_status;
  int n_nodes = A.n_cols;
  
  for (int i = 0; i < n_nodes; ++i) {
    node_status.push_back(0);
  }
  
  return node_status;
}

vector<int> shuffled_nodes(const sp_mat& A) {
  vector<int> sn;
  int n_nodes = A.n_cols;
  
  for (int i = 1; i <= n_nodes; ++i) {
    sn.push_back(i);
  }
  
  random_shuffle(sn.begin(), sn.end());
  
  return sn;
}

double adoption_probability(const sp_mat& A, const int& brand, const int& node, 
                            const vector<int>& node_status,
                            const double& delta, const vector<double>& q) {
  
  vector<int> nbrs = neighbours(A, node);
  int n_active_nbrs = 0;
  
  for (auto& nbr : nbrs) {
    if (node_status[nbr-1] == brand) n_active_nbrs ++;
  }
  
  return 1 - (1 - delta) * pow(1 - q[node-1], n_active_nbrs);
}

vector<int> random_seed(const sp_mat& A) {
  
}


/*** R

#' *This script attempts to replicate the results from "Decomposing the value of 
#' word-of-mouth seeding programs: Acceleration vs expansion", Libai, Muller and 
#' Peres (2010).*
#'
#' A new product is introduced into a competitive market via a *seeding program*
#' in which an initial group of influentials (the seed) receives the product
#' early on so that their word of mouth begins to drive diffusion. 
#' From a group of $g$ customers, a firm derives a direct value, which is the value
#' of their purchases and the social value, which is the effect of this customer
#' on the purchases of the rest of the customers.


*/
