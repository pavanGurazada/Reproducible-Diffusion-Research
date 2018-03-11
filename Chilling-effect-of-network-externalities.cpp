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

vector<int> neighbours(const sp_mat& A, const int& n_nodes, const int& node) {
  vector<int> nbrs; 
  
  for (int i = 0; i < n_nodes; ++i) {
    if (A(i, node-1) == 1) nbrs.push_back(i+1);
  }
  
  return nbrs;
}

vector<bool> reset_nodes(const sp_mat& A, const int& n_nodes) {
  vector<bool> node_status;

  for (int i = 0; i < n_nodes; ++i) {
    node_status.push_back(false);
  }
  
  return node_status;
}

vector<double> initialize_threshold(const sp_mat& A, const int& n_nodes, const double& mu, const double& sigma) {
  std::random_device rd;
  std::mt19937 gen(rd());
  
  normal_distribution<double> d(mu, sigma);
  
  vector<double> threshold;
  
  for (int i = 0; i < n_nodes; ++i) {
    threshold.push_back(d(gen));
  }
  
  return threshold;
}

vector<int> shuffled_nodes(const sp_mat& A, const int& n_nodes) {
  vector<int> sn;
  
  for (int i = 1; i <= n_nodes; ++i) {
    sn.push_back(i);
  }
  
  random_shuffle(sn.begin(), sn.end());
  
  return sn;
}

double network_externalities_effect(const sp_mat& A, const int& n_nodes, const int& node, 
                                    vector<bool>& node_status, 
                                    const vector<double>& threshold, 
                                    const double& a, const double& b) {
  vector<int> nbrs = neighbours(A, n_nodes, node);

  if (count(node_status.begin(), node_status.end(), true)/n_nodes > threshold[node-1]) {
    
    int n_active_nbrs = 0;
    
    for (auto& nbr : nbrs) {
      if (node_status[nbr-1] == true) n_active_nbrs++;
    }
    
    return 1 - (1 - a) * pow(1 - b, n_active_nbrs);
    
  } 
  
  else {
      
      return a;
  }
}


void evolve(const sp_mat& A, const int& n_nodes,
            vector<bool>& node_status, 
            const vector<double>& threshold,
            const double& a, 
            const double& b) {
  
  vector<int> s_nodes = shuffled_nodes(A, n_nodes);
  
  for (auto& node : s_nodes) {
    if (R::runif(0, 1) < network_externalities_effect(A, n_nodes, node, node_status, threshold, a, b))
      node_status[node-1] = true;
  }

}

void evolve(const sp_mat& A, const int& n_nodes,
            vector<bool>& node_status, 
            const vector<double>& threshold, 
            const double& a) {
  vector<int> s_nodes = shuffled_nodes(A, n_nodes);
  
  for (auto& node : s_nodes) {
    if (R::runif(0, 1) < a)
      node_status[node-1] = true;
  }
  
}


// [[Rcpp::export]]
DataFrame simulate(const sp_mat& A,
                   const double& a_i, 
                   const double& b_i,
                   const double& mu_i,
                   const double& sigma_i) {
 
  int T = 30, n_realizations = 10;
  const int n_nodes = A.n_cols;

  vector<bool> node_status_e, node_status_ne; 
  
  vector<int> realizations, time_steps, n_engaged_e, n_engaged_ne;
  vector<double> a, b, mu, sigma;
  
  for (int r = 1; r <= n_realizations; ++r) {

    node_status_e = reset_nodes(A, n_nodes);
    node_status_ne = reset_nodes(A, n_nodes);
    
    vector<double> threshold = initialize_threshold(A, n_nodes, mu_i, sigma_i);
    
    for (int t = 1; t <= T; ++t) {

      evolve(A, n_nodes, node_status_e, threshold, a_i, b_i);
      evolve(A, n_nodes, node_status_ne, threshold, a_i);

      realizations.push_back(r);
      time_steps.push_back(t);
      n_engaged_e.push_back(count(node_status_e.begin(), node_status_e.end(), true));
      n_engaged_ne.push_back(count(node_status_ne.begin(), node_status_ne.end(), true));
      a.push_back(a_i);
      b.push_back(b_i);
      mu.push_back(mu_i);
      sigma.push_back(sigma_i);
    }
  }
  

  DataFrame output = DataFrame::create(Named("realization") = realizations,
                                       Named("time_steps") = time_steps,
                                       Named("a") = a,
                                       Named("b") = b,
                                       Named("mu") = mu,
                                       Named("sigma") = sigma,
                                       Named("n_engaged_e") = n_engaged_e,
                                       Named("n_engaged_ne") = n_engaged_ne);

  return output;
}


 
/*** R

#' ---
#' title: "The chilling effect of network externalities"
#' author: Pavan Gurazada
#' output: github_document
#' ---
#' last update: Sat Mar 10 21:48:39 2018

library(igraph)

parameter_space <- expand.grid(a = seq(0.005, 0.05, length.out = 5),
                               b = seq(0.05, 0.25, length.out = 5),
                               mu = seq(0.01, 0.2, length.out = 5),
                               sigma = seq(0.005, 0.025, length.out = 5))

results <- data.frame()

g <- erdos.renyi.game(625, 8/625)
A <- get.adjacency(g)

cat("\n Starting simulation at ", date(), "\n")
for (row in 1:nrow(parameter_space)) {


  output <- simulate(A,
                     parameter_space[row, 1], parameter_space[row, 2],
                     parameter_space[row, 3], parameter_space[row, 4])

  results <- rbind(results, output)
}

cat("\n Ending simulation at ", date(), "\n")

*/
