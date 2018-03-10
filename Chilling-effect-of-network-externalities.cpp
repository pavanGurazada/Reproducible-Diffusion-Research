#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

using std::pow;
using std::accumulate;
using std::vector;
using std::random_shuffle;
using std::normal_distribution;

vector<int> neighbours(const sp_mat& A, const int& node) {
  vector<int> nbrs; 
  int n = A.n_rows;
  
  for (int i = 0; i < n; ++i) {
    if (A(i, node-1) == 1) nbrs.push_back(i+1);
  }
  
  return nbrs;
}

vector<bool> reset_nodes(sp_mat& A) {
  vector<bool> node_status;
  int n = A.n_rows;
  
  for (int i = 0; i < n; ++i) {
    node_status.push_back(false);
  }
  
  return node_status;
}

vector<double> initialize_threshold(sp_mat& A, double mu, double sigma) {
  std::random_device rd;
  std::mt19937 gen(rd());
  
  normal_distribution<double> d(mu, sigma);
  
  vector<double> threshold;
  int n = A.n_rows;
  
  for (int i = 0; i < n; ++i) {
    threshold.push_back(d(gen));
  }
  
  return threshold;
}

vector<int> shuffled_nodes(const sp_mat& A) {
  vector<int> sn;
  int n = A.n_rows;
  
  for (int i = 1; i <= n; ++i) {
    sn.push_back(i);
  }
  
  random_shuffle(sn.begin(), sn.end());
  
  return sn;
}

double network_externalities_effect(const sp_mat& A, const int& node, 
                                    vector<bool>& node_status, 
                                    const vector<double>& threshold, 
                                    double a, double b) {
  vector<int> nbrs = neighbours(A, node);
  int n = A.n_rows;

  if (accumulate(node_status.begin(), node_status.end(), 0)/n > threshold[node-1]) {
    
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

vector<bool> evolve_e(const sp_mat& A, 
                      vector<bool>& node_status, 
                      const vector<double>& threshold,
                      double a, 
                      double b) {
  
  vector<int> s_nodes = shuffled_nodes(A);
  
  for (auto& node : s_nodes) {
    if (R::runif(0, 1) < network_externalities_effect(A, node, node_status, threshold, a, b))
      node_status[node-1] = true;
  }
  
  return node_status;
}

vector<bool> evolve_ne(const sp_mat& A, 
                       vector<bool>& node_status, 
                       const vector<double>& threshold, 
                       double a) {
  vector<int> s_nodes = shuffled_nodes(A);
  
  for (auto& node : s_nodes) {
    if (R::runif(0, 1) < a)
      node_status[node-1] = true;
  }
  
  return node_status;
}


// [[Rcpp::export]]
DataFrame simulate(sp_mat& A,
                   const double& a_i, 
                   const double& b_i,
                   const double& mu_i,
                   const double& sigma_i) {
 
  int T = 30, n_realizations = 10;

  vector<bool> node_status_e, node_status_ne; 
  
  vector<int> realizations, time_steps, n_engaged_e, n_engaged_ne;
  vector<double> a, b, mu, sigma;
  
  for (int r = 1; r <= n_realizations; ++r) {

    node_status_e = reset_nodes(A);
    node_status_ne = reset_nodes(A);
    
    vector<double> threshold = initialize_threshold(A, mu_i, sigma_i);
    
    for (int t = 1; t <= T; ++t) {
      evolve_e(A, node_status_e, threshold, a_i, b_i);
      evolve_ne(A, node_status_ne, threshold, a_i);
    
      realizations.push_back(r);
      time_steps.push_back(t);
      n_engaged_e.push_back(accumulate(node_status_e.begin(), node_status_e.end(), 0));
      n_engaged_ne.push_back(accumulate(node_status_ne.begin(), node_status_ne.end(), 0));
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
#' last update: Sat Mar 10 14:19:21 2018

library(igraph)

parameter_space <- expand.grid(a = seq(0.005, 0.05, length.out = 5),
                               b = seq(0.05, 0.25, length.out = 5),
                               mu = seq(0.01, 0.2, length.out = 5),
                               sigma = seq(0.005, 0.025, length.out = 5))

results <- data.frame()

g <- erdos.renyi.game(625, 8/625)
A <- get.adjacency(g)

for (row in 1:nrow(parameter_space)) {
  cat("\n Starting simulation on setting :", row, " at ", date(), "\n")

  output <- simulate(A,
                     parameter_space[row, 1], parameter_space[row, 2],
                     parameter_space[row, 3], parameter_space[row, 4])

  results <- rbind(results, output)
}

*/
