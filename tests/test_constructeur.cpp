/* 
 * File:   test_constructeur.cpp
 * Author: gononl
 *
 * Created on September 23, 2016, 5:50 PM
 */
#include "../src/MonteCarlo.hpp"
#include "../src/BlackScholesModel.hpp"

#include <cstdlib>
#include <iostream>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    MonteCarlo *mc = new MonteCarlo();
    PnlMat *path = pnl_mat_create(mc->opt_->nbTimeSteps_ + 1, mc->mod_->size_);
    
    mc->mod_->asset(path, mc->opt_->T_, mc->opt_->nbTimeSteps_, mc->rng_);
    cout << "Chemin de simulation : " << endl;
    pnl_mat_print(path);
    
    PnlMat *past = pnl_mat_create(1, mc->mod_->size_);
    PnlVect *delta = pnl_vect_create(mc->mod_->size_);
    mc->delta(past, 0, delta);
    
    cout << "Delta : " << endl;
    pnl_vect_print(delta);
    return 0;
}

