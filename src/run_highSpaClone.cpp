#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

double calculateObject(arma::sp_mat &Y, arma::mat &V, arma::sp_mat &A, arma::mat &B, double lambda) {
  arma::sp_mat colsum_A = arma::sum(A, 1);
  arma::sp_mat D((int)V.n_rows, (int)V.n_rows);
  D.diag() = colsum_A;
  arma::sp_mat L = D - A;

  double Obj_NMF = norm(Y - B * V.t(), "fro");
  double traceTerm = trace(V.t() * L * V);
  double obj = Obj_NMF * Obj_NMF + lambda * traceTerm;

  return obj;
}

arma::mat calc_VD(arma::sp_mat &colsum_A, arma::mat &V, arma::uvec &indexStr){
  arma::vec colsum_A_dense = vec(colsum_A.col(0));
  arma::mat res = V.rows(indexStr);
  res.each_col() %= colsum_A_dense.elem(indexStr);
  return res;
}

// [[Rcpp::export]]
SEXP run_iter(arma::sp_mat &Y, SEXP BIn, arma::sp_mat &A, arma::mat V, arma::uvec vecStr, double lambda) {
  try {

    arma::mat B = as<arma::mat>(BIn);
    arma::sp_mat colsum_A = arma::sum(A, 1);
    arma::sp_mat D((int)V.n_rows, (int)V.n_rows);
    D.diag() = colsum_A;

    arma::uvec UniqueStr = arma::unique(vecStr);

    for (int istr : UniqueStr) {
      arma::uvec indexStr = find(vecStr == istr);
      arma::mat VD = calc_VD(colsum_A, V, indexStr);
      arma::mat VA = A.cols(indexStr).t() * V;
      arma::mat nomUpdateVStr = Y.cols(indexStr).t() * B + lambda * VA;
      arma::mat denomUpdateVStr = V.rows(indexStr) * B.t() * B + lambda * VD;

      arma::mat updateVStr = nomUpdateVStr / denomUpdateVStr;
      V.rows(indexStr) = V.rows(indexStr) % updateVStr;

    }

    // Objective value
    double obj = calculateObject(Y, V, A, B, lambda);

    return Rcpp::List::create(
      Rcpp::Named("V") = V,
      Rcpp::Named("Obj") = obj
    );
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
}

