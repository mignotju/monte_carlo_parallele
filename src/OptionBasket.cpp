/* 
 * File:   OptionBasket.cpp
 * Author: gononl
 * 
 * Created on September 16, 2016, 10:03 AM
 */
#include <iostream>
#include "OptionBasket.hpp"

using namespace std;

OptionBasket::OptionBasket() {
    T_ = 10;
    nbTimeSteps_ = 10;
    size_ = 3;
    lambda = pnl_vect_create(size_);
    pnl_vect_set(lambda,0,0.2);
    pnl_vect_set(lambda,1,0.7);
    pnl_vect_set(lambda,2,0.1);
    K = 100;   
    ST = pnl_vect_create(size_);
}

OptionBasket::OptionBasket(double maturity, int nbTimeSteps, int size, double strike, PnlVect* lambdaOption) {
    T_ = maturity;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    K = strike;
    lambda = lambdaOption;
    ST = pnl_vect_create(size_);
}

double OptionBasket::payoff(const PnlMat* path) {
    double payoff = 0;
    pnl_mat_get_row(ST, path, path->m - 1);
    payoff = pnl_vect_scalar_prod(lambda, ST);

    payoff -= K;
    if (payoff < 0) {
        return 0;
    } else {
        return payoff;
    }   
}
 
