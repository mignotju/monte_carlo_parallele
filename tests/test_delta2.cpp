/* 
 * File:   test_delta2.cpp
 * Author: gononl
 *
 * Created on September 21, 2016, 2:41 PM
 */

#include <cstdlib>
#include <iostream>

#include "../src/MonteCarlo.hpp"
#include "../src/BlackScholesModel.hpp"
#include "../src/OptionBasket.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    // Black and Scholes
    BlackScholesModel *bs = new BlackScholesModel();
    bs->size_ = 1;
    bs->r_ = 0.1;
    bs->rho_ = 0;
    bs->sigma_ = pnl_vect_create_from_scalar(bs->size_, 0.2);
    bs->spot_ = pnl_vect_create_from_scalar(bs->size_, 100);
    
    // Option
    double maturity = 2;
    int nbTimeSteps = 2;
    int size = 1;
    double strike = 10;
    PnlVect *lambda = pnl_vect_create_from_scalar(size, 1);

    OptionBasket opt(maturity, nbTimeSteps, size, strike, lambda);
    
    PnlRng *rng = pnl_rng_new();
    PnlMat *past = pnl_mat_new();
    
    // Creation du vecteur deterministe 
    PnlVect *vect = pnl_vect_create_from_scalar(size, 0.5);
    
    // Monte Carlo
    MonteCarlo *mc = new MonteCarlo();
    mc->fdStep_ = 0.005;
    mc->mod_ = bs;
    mc->nbSamples_ = 1;
    mc->opt_ = &opt;
    mc->rng_ = rng;
    
    past = pnl_mat_create_from_scalar(1, 1, 100);

    PnlVect *vectDelta = pnl_vect_create(bs->size_);
    mc->delta(past, 0, vectDelta, vect);
    cout << "Delta final : " << endl;
    pnl_vect_print(vectDelta);
    return 0;
}

