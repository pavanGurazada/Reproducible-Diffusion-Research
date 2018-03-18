#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using arma::sp_mat;
using arma::vec;

using std::vector;

struct BassModel {
  sp_mat A;
  float p;
  float q;
};

vector<int> neighbours(const sp_mat& A, const int& node) {
  vector<int> nbrs; 
  
  for (int i = 0; i < A.n_cols; ++i) {
    if (A(i, node-1) == 1) nbrs.push_back(i+1);
  }
  
  return nbrs;
}

vector<int> reset_nodes(BassModel& M) {
  vector<int> node_status;
  
  for (int i = 0; i < M.A.n_cols; ++i) {
    node_status.push_back(0);
  }
  
  return node_status;
}

double adoption_prob(vec& node_status, BassModel& M, int node) {
  int n_adopted_nbrs = 0;
  
  for (auto& nbr : neighbours(M.A, node)) {
    
    if (node_status[node] == 1) n_adopted_nbrs++;
    
  }
  
  return M.q * n_adopted_nbrs/M.A.n_cols;
    
}


void evolve(vec& node_status, BassModel& M){
  
  
  
}


/*** R
timesTwo(42)
*/
