/* 
 * File:   OptionPerformance.cpp
 * Author: gononl
 * 
 * Created on September 17, 2016, 5:06 PM
 */
#include <iostream>
#include <stdexcept>
#include "OptionPerformance.hpp"

using namespace std;

OptionPerformance::OptionPerformance() {
    T_ = 10;
    nbTimeSteps_ = 10;
    size_ = 2;
    lambda = pnl_vect_create_from_scalar(size_, 1/nbTimeSteps_);
    cours = pnl_vect_create(size_);
}

OptionPerformance::OptionPerformance(double maturity, int nbTimeSteps, int size, PnlVect* lambdaOption) {
    T_ = maturity;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    lambda = lambdaOption;
    cours = pnl_vect_create(size_);
}

double OptionPerformance::payoff(const PnlMat* path) {
    double numerateur = 0;
    double denominateur = 0;
    double sommePartielle = 0;

    for (int i = 1; i < path->m; i++) {
        pnl_mat_get_row(cours, path, i);
        numerateur = pnl_vect_scalar_prod(lambda, cours);
        pnl_mat_get_row(cours, path, i-1);
        denominateur = pnl_vect_scalar_prod(lambda, cours);
        numerateur /= denominateur;
        numerateur -= 1;
        if (numerateur > 0) {
            sommePartielle += numerateur;
        }
    }
    return (1 + sommePartielle);    
}