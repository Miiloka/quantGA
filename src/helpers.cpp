#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double checkFun(const arma::vec& residuals, double tau) {
    if (tau <= 0.0 || tau >= 1.0) {
        Rcpp::stop("tau must be between 0 and 1");
    }
    
    double result = 0.0;
    int n = residuals.n_elem;
    
    for (int i = 0; i < n; i++) {
        if (residuals(i) >= 0.0) {
            result += residuals(i) * tau;
        } else {
            result += residuals(i) * (tau - 1.0);
        }
    }
    
    return result;
}

// [[Rcpp::export]]
double objFun(const arma::vec& beta, const arma::mat& X, const arma::vec& y, double tau) {
    if (beta.n_elem != X.n_cols) {
        return std::numeric_limits<double>::max();
    }
    
    try {
        arma::vec fitted_values = X * beta;
        
        if (!fitted_values.is_finite()) {
            return std::numeric_limits<double>::max();
        }
        
        arma::vec residuals = y - fitted_values;
        
        return checkFun(residuals, tau);
        
    } catch (...) {
        return std::numeric_limits<double>::max();
    }
}

// [[Rcpp::export]]
arma::vec popFitness(const arma::mat& pop, const arma::mat& X, const arma::vec& y, double tau) {
    int pop_size = pop.n_rows;
    arma::mat resid_mat = -pop * X.t();
    resid_mat.each_row() += y.t();

    arma::mat weights = arma::conv_to<arma::mat>::from(resid_mat >= 0.0);
    weights = weights * tau + (1.0 - weights) * (tau - 1.0);

    arma::vec fitness = arma::sum(resid_mat % weights, 1);

    return fitness;
}
