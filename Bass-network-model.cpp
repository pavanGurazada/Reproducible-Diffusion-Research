#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using Rcpp::DataFrame;
using Rcpp::Environment;
using Rcpp::Function;
using Rcpp::Named;
using Rcpp::S4;
using Rcpp::as;
using Rcpp::Rcout;

using arma::sp_mat;
using arma::vec;
using arma::randu;

using std::vector;

struct BassModel {
  sp_mat A;
  double p;
  double q;
};

vector<int> neighbours(const sp_mat& A, const int& node) {
  vector<int> nbrs; 
  
  for (std::size_t i = 0; i < A.n_cols; ++i) {
    if (A(i, node-1) == 1) nbrs.push_back(i+1);
  }
  
  return nbrs;
}

vector<int> shuffled_nodes(const sp_mat& A) {
  vector<int> sn;
  
  for (std::size_t i = 1; i <= A.n_cols; ++i) {
    sn.push_back(i);
  }
  
  random_shuffle(sn.begin(), sn.end());
  
  return sn;
}

vector<int> reset_nodes(BassModel& M) {
  vector<int> node_status;
  
  for (std::size_t i = 0; i < M.A.n_cols; ++i) {
    node_status.push_back(0);
  }
  
  return node_status;
}

double adoption_prob(vector<int>& node_status, BassModel& M, const int& node) {
  int n_adopted_nbrs = 0;
  
  for (auto& nbr : neighbours(M.A, node)) {
    
    if (node_status[nbr-1] == 1) n_adopted_nbrs++;
    
  }
  
  return M.q * n_adopted_nbrs/M.A.n_cols;
    
}

void evolve(vector<int>& node_status, BassModel& M){
  
  for (auto& node: shuffled_nodes(M.A)) {
    
    if (node_status[node-1] == 0) {
      if ((randu() < M.p) || (randu() < adoption_prob(node_status, M, node)))
        node_status[node-1] = 1;
    }
  }
}

sp_mat get_adjacency_er(int n, double p) {
  Environment igraph("package:igraph");
  Function game_er = igraph["erdos.renyi.game"];
  Function adj_mat = igraph["get.adjacency"];
  
  SEXP g = game_er(Named("n", n), Named("p", p));
  S4 Am = adj_mat(Named("g", g));
  
  sp_mat A = as<sp_mat>(Am);
  
  return A;
}

sp_mat get_adjacency_ba(int n, int m) {
  Environment igraph("package:igraph");
  Function game_ba = igraph["sample_pa"];
  Function adj_mat = igraph["get.adjacency"];
  
  SEXP g = game_ba(Named("n", n), Named("m", m));
  S4 Am = adj_mat(Named("g", g));
  
  sp_mat A = as<sp_mat>(Am);
  
  return A;
}

// [[Rcpp::export]]
DataFrame simulate_er(int n_i, int k_i, double p_i, double q_i, 
                      const int n_realizations = 10, const int T = 30) {
  
  vector<int> n, k, realizations, time_steps, n_engaged;
  vector<double> p, q;
  
  for (int r = 1; r <= n_realizations; ++r) {
    Rcout << "Starting the realization :" << r << std::endl;
    
    sp_mat A = get_adjacency_er(n_i, k_i/n_i);
    BassModel M = {A, p_i, q_i};
    
    vector<int> node_status = reset_nodes(M);
    
    for (int t = 1; t <= T; ++t) {
      evolve(node_status, M);
      
      n.push_back(n_i);
      k.push_back(k_i);
      p.push_back(p_i);
      q.push_back(q_i);
      realizations.push_back(r);
      time_steps.push_back(t);
      n_engaged.push_back(accumulate(node_status.begin(), node_status.end(), 0));
    }
  }
  
  DataFrame output = DataFrame::create(Named("n") = n,
                                       Named("k") = k,
                                       Named("p") = p,
                                       Named("q") = q,
                                       Named("realization") = realizations,
                                       Named("time_steps") = time_steps,
                                       Named("n_engaged") = n_engaged);
  
  return output;
}

// [[Rcpp::export]]
DataFrame simulate_ba(int n_i, int m_i, double p_i, double q_i, 
                      const int n_realizations = 10, const int T = 30) {
  
  vector<int> n, m, realizations, time_steps, n_engaged;
  vector<double> p, q;
  
  for (int r = 1; r <= n_realizations; ++r) {
    //Rcout << "Starting the realization :" << r << std::endl;
    
    sp_mat A = get_adjacency_ba(n_i, m_i);
    BassModel M = {A, p_i, q_i};
    
    vector<int> node_status = reset_nodes(M);
    
    for (int t = 1; t <= T; ++t) {
      evolve(node_status, M);
      
      n.push_back(n_i);
      m.push_back(m_i);
      p.push_back(p_i);
      q.push_back(q_i);
      realizations.push_back(r);
      time_steps.push_back(t);
      n_engaged.push_back(accumulate(node_status.begin(), node_status.end(), 0));
    }
  }
  
  DataFrame output = DataFrame::create(Named("n") = n,
                                       Named("m") = m,
                                       Named("p") = p,
                                       Named("q") = q,
                                       Named("realization") = realizations,
                                       Named("time_steps") = time_steps,
                                       Named("n_engaged") = n_engaged);
  
  return output;
}

/*** R
library(igraph)
results <- data.frame(n = integer(), m = integer(), p = numeric(), q = numeric(),
                      realization = integer(), time_steps = integer(), n_engaged = integer())

for (m in 1:10) {
  cat("[Refrigerators] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(400, m, 0.0026167, 0.215666)) # Refrigerators
  cat("[Freezers] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(220, m, 0.0181190, 0.17110)) # Freezers
  cat("[TV] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(970, m, 0.0278770, 0.25105)) # TV
  cat("[Softener] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(600, m, 0.0177030, 0.29695)) # Softener
  cat("[AC] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(170, m, 0.0103990, 0.41861)) # AC
  cat("[Dryer] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(151, m, 0.0172, 0.35688)) # Dryer
  cat("[Lawnmowers] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(450, m, 0.0091837, 0.33790)) # Lawnmower
  cat("[Bed] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(765, m, 0.0058760, 0.24387)) # Bed
  cat("[Coffee] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(589, m, 0.0171350, 0.30145)) # Coffee
  cat("[Iron] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(557, m, 0.0286320, 0.32791)) # Iron
  cat("[Player] Starting simulations for min degree : ", m, " at ", format(Sys.time(), "%H:%M"), "\n")
  results <- rbind(results, simulate_ba(220, m, 0.0247960, 0.65410)) # Player
}

*/
