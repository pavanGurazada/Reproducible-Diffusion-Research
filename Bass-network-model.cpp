#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

/**
 * @title The Bass model on a network
 * @author Pavan Gurazada
 * @licence GPL (>=2)
 * @summary Demonstrates the execution of a complex model of network diffusion.
 *   Replicates the results from "Agend-based models in marketing: Guidelines
 *   for rigor", Rand and Rust (2011).
 * 
 */

/**
 * # Introduction
 * 
 * Rand and Rust (2011) establish a blue-print for executing diffusion
 * simulation studies and formalize several steps that are required to be
 * followed in the execution of such models. The explicit time dependence of
 * diffusion simulations makes the execution times in R slow. While execution
 * times can be sped up using vectorized code, there are fundamental assumptions
 * (like that of synchronous updates) that violate the basic tenets of diffusion
 * processes.
 *
 * This script presents a framework to set up and execute diffusion simulations
 * using `Rcpp`. The `igraph` `R` library is used to generate built-in network
 * structures and the entire simulation is executed in `C++` for speed.
 *
 * Our conceptualization of diffusion simulations is composed of a set of verbs
 * - `initialize`, `reset`, `seed`, `evolve` and `simulate`(hat-tip to `dplyr`
 * and `particles`). The process of simulation (executed by the exported
 * function `simlate`) is at its core a repeated execution of the following
 * sequence:
 *
 * `reset(...) %>% seed(...) %>% evolve`
 *
 * Along the way, we also present several helper functions that build bridges between
 * the functions within R packages and C++.
 */

using namespace Rcpp;

/**
 * The model we wish to replicate is an extension of the traditional Bass model
 * where the diffusion of products is dependent on three parameters:
 *
 * 1. The social network of the users in the market (represented by an adjacency
 * matrix `A`) 
 * 2. The influence of external communication (e.g., advertising) on
 * adoption (represented by `p`) 
 * 3. The influence of word-of-mouth from friends
 * (represented by `q`)
 * 
 * Consequently, the model is conceptualized and stored as a struct with these
 * parameters that are initialized at each realization of the simulation.
 */

struct BassModel {
  arma::sp_mat A;
  double p;
  double q;
};

/**
 * # Overview of key functions
 * 
 * The first helper function `neighbors(...)` returns a vector of neighbors for a
 * given node. Parallel functionality exists in `igraph`, but we re-create
 * this C++ function to avoid the message passing overhead between R and C++.
 *
 * This function indexes into the corresponding column of the adjacency matrix
 * and returns the indices of all non-zero entries in the column
 */

std::vector<int> neighbours(const arma::sp_mat& A, const int& node) {
  std::vector<int> nbrs; 
  
  for (std::size_t i = 0; i < A.n_cols; ++i) {
    if (A(i, node-1) == 1) nbrs.push_back(i+1);
  }
  
  return nbrs;
}

/**
 * The second helper function `shuffled_nodes(...)` returns a sample vector
 * of nodes drawn at random from the original set of vertices in the network.
 * This is used to randomize the order of updation of the node status at each 
 * step of the diffusion process. 
 */

std::vector<int> shuffled_nodes(const arma::sp_mat& A) {
  std::vector<int> sn;
  
  for (std::size_t i = 1; i <= A.n_cols; ++i) {
    sn.push_back(i);
  }
  
  std::random_shuffle(sn.begin(), sn.end());
  
  return sn;
}

/**
 * The third helper function `reset_nodes(...)` returns a vector of 0's
 * indicating the false adoption status of all nodes at the beginning of each
 * diffusion process
 * 
 */

std::vector<int> reset_nodes(BassModel& M) {
  std::vector<int> node_status;
  
  for (std::size_t i = 0; i < M.A.n_cols; ++i) {
    node_status.push_back(0);
  }
  
  return node_status;
}

double adoption_prob(std::vector<int>& node_status, BassModel& M, const int& node) {
  int n_adopted_nbrs = 0;
  
  for (auto& nbr : neighbours(M.A, node)) {
    
    if (node_status[nbr-1] == 1) n_adopted_nbrs++;
    
  }
  
  return M.q * n_adopted_nbrs/M.A.n_cols;
    
}

void evolve(std::vector<int>& node_status, BassModel& M){
  
  for (auto& node: shuffled_nodes(M.A)) {
    
    if (node_status[node-1] == 0) {
      if ((arma::randu() < M.p) || (arma::randu() < adoption_prob(node_status, M, node)))
        node_status[node-1] = 1;
    }
  }
}

arma::sp_mat get_adjacency_er(int n, double p) {
  Environment igraph("package:igraph");
  Function game_er = igraph["erdos.renyi.game"];
  Function adj_mat = igraph["get.adjacency"];
  
  SEXP g = game_er(Named("n", n), Named("p", p));
  S4 Am = adj_mat(Named("g", g));
  
  arma::sp_mat A = as<arma::sp_mat>(Am);
  
  return A;
}

arma::sp_mat get_adjacency_ba(int n, int m) {
  Environment igraph("package:igraph");
  Function game_ba = igraph["sample_pa"];
  Function adj_mat = igraph["get.adjacency"];
  
  SEXP g = game_ba(Named("n", n), Named("m", m));
  S4 Am = adj_mat(Named("g", g));
  
  arma::sp_mat A = as<arma::sp_mat>(Am);
  
  return A;
}

// [[Rcpp::export]]
DataFrame simulate_er(int n_i, int k_i, double p_i, double q_i, 
                      const int n_realizations = 10, const int T = 30) {
  
  std::vector<int> n, k, realizations, time_steps, n_engaged;
  std::vector<double> p, q;
  
  for (int r = 1; r <= n_realizations; ++r) {
    Rcout << "Starting the realization :" << r << std::endl;
    
    arma::sp_mat A = get_adjacency_er(n_i, k_i/n_i);
    BassModel M = {A, p_i, q_i};
    
    std::vector<int> node_status = reset_nodes(M);
    
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
  
  std::vector<int> n, m, realizations, time_steps, n_engaged;
  std::vector<double> p, q;
  
  for (int r = 1; r <= n_realizations; ++r) {
    //Rcout << "Starting the realization :" << r << std::endl;
    
    arma::sp_mat A = get_adjacency_ba(n_i, m_i);
    BassModel M = {A, p_i, q_i};
    
    std::vector<int> node_status = reset_nodes(M);
    
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

result_summary <- results %>% group_by(n, m, p, q, time_steps) %>% 
                             summarize(mean_engaged = mean(n_engaged)) %>% 
                             group_by(n, m, p, q) %>% 
                             mutate(incremental_engaged = mean_engaged - lag(mean_engaged))

*/
