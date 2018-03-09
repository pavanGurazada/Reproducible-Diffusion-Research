#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using std::pow;
using std::accumulate;

// [[Rcpp::export]]
IntegerVector neighbours(const sp_mat& A, const int& node) {
  IntegerVector nbrs; 
  
  sp_mat::const_row_iterator start = A.begin_row(node-1);
  sp_mat::const_row_iterator end = A.end_row(node-1);
  
  for (sp_mat::const_row_iterator& i = start; i != end; ++i) {
      nbrs.push_back(i.col() + 1);
  }
  
  return nbrs;
}

// [[Rcpp::export]]
int n_nodes(const sp_mat& A) {
  int n = A.n_rows;
  
  return n;
}

// [[Rcpp::export]]
LogicalVector reset_nodes(sp_mat& A) {
  LogicalVector node_status;
  
  for (int i = 0; i < n_nodes(A); ++i) {
    node_status.push_back(false);
  }
  
  return node_status;
}

// [[Rcpp::export]]
IntegerVector shuffled_vertices(const sp_mat& A) {
  IntegerVector sv;
  
  for (int i = 1; i <= n_nodes(A); ++i) {
    sv.push_back(i);
  }
  
  return sample(sv, sv.length(), false);
}

// [[Rcpp::export]]
double network_externalities_effect(const sp_mat& A, const int& node, 
                                    LogicalVector& node_status, 
                                    const NumericVector& threshold, 
                                    double a, double b) {
  IntegerVector nbrs = neighbours(A, node);

  if (sum(node_status)/n_nodes(A) > threshold[node-1]) {
    
    int n_active_nbrs = 0;
    
    for (auto& nbr : nbrs) {
      if (node_status[nbr-1] == true) n_active_nbrs++;
    }
    
    return 1 - (1 - a) * pow(1 - b, n_active_nbrs);
    
  } else {
    
    return a;
  }
}

// [[Rcpp::export]]
LogicalVector evolve_e(const sp_mat& A, 
                       LogicalVector& node_status, 
                       const NumericVector& threshold,
                       double a, 
                       double b) {
  
  for (auto& node : shuffled_vertices(A)){
    if (randu() < network_externalities_effect(A, node, node_status, threshold, a, b))
      node_status[node-1] = true;
  }
  
  return node_status;
}

// [[Rcpp::export]]
LogicalVector evolve_ne(const sp_mat& A, 
                        LogicalVector& node_status, 
                        const NumericVector& threshold, 
                        double a) {
  
  for (auto& node : shuffled_vertices(A)) {
    if (randu() < a)
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
 
  int T = 30, n_realizations = 5;
  int data_points = T * n_realizations;
  
  LogicalVector node_status_e; 
  LogicalVector node_status_ne; 
  
  IntegerVector realizations(data_points), time_steps(data_points), n_engaged_e(data_points), n_engaged_ne(data_points);
  NumericVector a(data_points), b(data_points), mu(data_points), sigma(data_points);
  
  int i = 0; // indexing the vectors
  
  for (int r = 1; r <= n_realizations; ++r) {
    
    Rcout << "Realization: " << r << std::endl;
    
    node_status_e = reset_nodes(A);
    node_status_ne = reset_nodes(A);
    
    NumericVector threshold = rnorm(n_nodes(A), mu_i, sigma_i);
    
    for (int t = 1; t <= T; ++t) {
      evolve_e(A, node_status_e, threshold, a_i, b_i);
      evolve_ne(A, node_status_ne, threshold, a_i);
      //NumericVector v = {r, t, a, b, mu, sigma, sum(node_status_e), sum(node_status_ne)};
      
      Rcout << "Time step: " << t << std::endl;
      
      realizations[i] = r;
      time_steps[i] = t;
      n_engaged_e[i] = sum(node_status_e);
      n_engaged_ne[i] = sum(node_status_ne);
      //n_engaged_e[i] = accumulate(node_status_e.begin(), node_status_e.end(), 0);
      //n_engaged_ne[i] = accumulate(node_status_ne.begin(), node_status_ne.end(), 0);
      a[i] = a_i;
      b[i] = b_i;
      mu[i] = mu_i;
      sigma[i] = sigma_i;
      i++;
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
#' last update: Thu Mar 08 12:17:12 2018


# library(igraph)
# 
# parameter_space <- expand.grid(a = seq(0.005, 0.05, length.out = 5),
#                                b = seq(0.05, 0.25, length.out = 5),
#                                mu = seq(0.01, 0.2, length.out = 5),
#                                sigma = seq(0.005, 0.025, length.out = 5))
# 
# n_realizations <- 10
# 
# results <- data.frame()
# 
# g <- erdos.renyi.game(625, 8/625)
# A <- get.adjacency(g)
# 
# for (row in 1:nrow(parameter_space)) {
#   cat("\n Starting simulation on setting :", row, " at ", date(), |)
# 
#   output <- simulate(A,
#                      parameter_space[row, 1], parameter_space[row, 2],
#                      parameter_space[row, 3], parameter_space[row, 4])
# 
#   results <- rbind(results, output)
# }

*/