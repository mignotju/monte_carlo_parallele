/* 
 * File:   OptionAsiatique.cpp
 * Author: gononl
 * 
 * Created on September 17, 2016, 4:27 PM
 */
#include <iostream>

#include "OptionAsiatique.hpp"

using namespace std;

OptionAsiatique::OptionAsiatique() {
    T_ = 10;
    nbTimeSteps_ = 10;
    size_ = 2;
    lambda = pnl_vect_create_from_scalar(size_, 1/nbTimeSteps_);
    K = 100;
    sommeStd = pnl_vect_create(size_);
}

OptionAsiatique::OptionAsiatique(double maturity, int nbTimeSteps, int size, double strike, PnlVect* lambdaOption) {
    T_ = maturity;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    K = strike;
    lambda = lambdaOption;
    sommeStd = pnl_vect_create(size_);
}

double OptionAsiatique::payoff(const PnlMat* path) {
    double payoff = 0;
    pnl_mat_sum_vect(sommeStd, path, 'r');
    pnl_vect_div_scalar(sommeStd, (double) path->m);
    payoff = pnl_vect_scalar_prod(lambda, sommeStd);    
    payoff -= K;
    if (payoff < 0) {
        return 0;
    } else {
        return payoff;
    }
}