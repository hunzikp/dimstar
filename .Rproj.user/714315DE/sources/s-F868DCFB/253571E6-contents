#define ARMA_64BIT_WORD 1

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

class Kernel {
public:

  mat ystar_sample;
  int M;
  vector<sp_mat> TM_vec;
  vector<sp_mat> W_vec;
  vector<sp_mat> O_vec;
  vector<sp_mat> W_t_vec;
  vector<sp_mat> O_t_vec;
  mat X;
  int G;
  int n_pairs;
  int NGT;
  int NG;

  Kernel(List W_ls, List TM_ls, List O_ls, mat X, int NG) : X(X), NG(NG)  {

    G = W_ls.size();
    for (int k = 0; k < G; ++k) {
      sp_mat W_k = W_ls.at(k);
      sp_mat TM_k = TM_ls.at(k);

      TM_vec.push_back(TM_k);
      W_vec.push_back(W_k);

      W_t_vec.push_back(W_k.submat(0, 0, NG-1, NG-1));
    }
    NGT = W_vec.at(0).n_cols;

    n_pairs = O_ls.size();
    for (int k = 0; k < n_pairs; ++k) {
      sp_mat O_k = O_ls.at(k);
      O_vec.push_back(O_k);

      O_t_vec.push_back(O_k.submat(0, 0, NG-1, NG-1));
    }
  };

  sp_mat make_A(NumericVector rho_vec, NumericVector gamma_vec, NumericVector lambda_vec) {
    sp_mat A = speye(NGT, NGT);
    for (int k = 0; k < G; ++k) {
      A -= rho_vec.at(k)*W_vec.at(k) + gamma_vec.at(k)*TM_vec.at(k);
    }
    for (int k = 0; k < n_pairs; ++k) {
      A -= lambda_vec.at(k)*O_vec.at(k);
    }
    return A;
  }

  sp_mat make_A_t(NumericVector rho_vec, NumericVector lambda_vec) {
    sp_mat A = speye(NG, NG);
    for (int k = 0; k < G; ++k) {
      A -= rho_vec.at(k)*W_t_vec.at(k);
    }
    for (int k = 0; k < n_pairs; ++k) {
      A -= lambda_vec.at(k)*O_t_vec.at(k);
    }
    return A;
  }

  void update_ystar_sample(mat _ystar_sample) {
    ystar_sample = _ystar_sample;
    M = ystar_sample.n_rows;
  }

  double get_E_kernel(NumericVector rho_vec, NumericVector gamma_vec, NumericVector lambda_vec,
                      colvec isigma, colvec beta) {

    double E_kernel = 0;
    mat Xbeta = X*beta;
    sp_mat A = make_A(rho_vec, gamma_vec, lambda_vec);

    for (int i = 0; i < M; ++i) {
      colvec ystar_i = ystar_sample.row(i).t();
      colvec eps = Xbeta - A*ystar_i;
      E_kernel += accu(isigma % pow(eps, 2));
    }
    E_kernel = E_kernel/double(M);

    return E_kernel;
  }

  double get_kernel_i(NumericVector rho_vec, NumericVector gamma_vec, NumericVector lambda_vec,
                      colvec isigma, colvec beta, int i) {

    mat Xbeta = X*beta;
    sp_mat A = make_A(rho_vec, gamma_vec, lambda_vec);

    colvec ystar_i = ystar_sample.row(i).t();
    colvec eps = Xbeta - A*ystar_i;
    double kern = accu(isigma % pow(eps, 2));

    return kern;
  }
};

RCPP_MODULE(kernel) {
  class_<Kernel>("Kernel")
  .constructor<List, List, List, arma::mat, int>()
  .method("update_ystar_sample", &Kernel::update_ystar_sample)
  .method("get_E_kernel", &Kernel::get_E_kernel)
  .method("make_A_t", &Kernel::make_A_t)
  .method("make_A", &Kernel::make_A)
  .method("get_kernel_i", &Kernel::get_kernel_i)
  ;
}


