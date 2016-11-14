/* 
 * File:   test_MonteCarlo.cpp
 * Author: girerda
 *
 * Created on September 19, 2016, 11:20 AM
 */

#include "../src/MonteCarlo.hpp"
#include "../src/BlackScholesModel.hpp"
#include "../src/OptionBasket.hpp"
#include "../src/OptionPerformance.hpp"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>

/*
 * Simple C++ Test Suite
 */

void valid_price_at_0_test() {
    std::cout << "test_MonteCarlo test prix en 0 - Basket" << std::endl;
    //init option
    double maturity = 2;
    int nbTimeSteps = 2;
    int size = 1;
    double strike = 10;
    PnlVect *lambda = pnl_vect_create_from_scalar(size, 1);
    OptionBasket opt(maturity, nbTimeSteps, size, strike, lambda);
    
    //init modèle
    BlackScholesModel* bs = new BlackScholesModel();
    bs->size_ = 1;
    bs->r_ = 0.1;
    bs->rho_ = 0.3;
    bs->sigma_ = pnl_vect_create_from_scalar(bs->size_, 0.2);
    bs->spot_ = pnl_vect_create_from_scalar(bs->size_, 100);
    bs->mat_cholesky = pnl_mat_create_from_scalar(bs->size_, bs->size_, bs->rho_);
    pnl_mat_set_diag(bs->mat_cholesky, 1, 0);
    pnl_mat_chol(bs->mat_cholesky);
    
    //Calcul des valeurs attendues
    PnlVect *G = pnl_vect_create_from_scalar(bs->size_,0.5);
    PnlMat *path = pnl_mat_create(opt.nbTimeSteps_+1, bs->size_);
    
    MonteCarlo *mc = new MonteCarlo();
    mc->fdStep_ = 1;
    mc->mod_ = bs;
    mc->nbSamples_ = 1000;
    mc->opt_ = &opt;
    
    double expected = 0;
        
    for (int i=0; i<mc->nbSamples_; i++) {
        bs->asset(path, opt.T_, opt.nbTimeSteps_, G);
        
        expected += opt.payoff(path);
    }
    
    expected *= exp(-bs->r_*(opt.T_ - 0)) / mc->nbSamples_;
    
    double obtained=0;
    double ic=0;
    mc->price(obtained,ic,G);
    
    if (expected == obtained) {
        std::cout << "    TEST MonteCarlo 1 OK : "<< expected << " == " << obtained << std::endl;
    } else {
        std::cout << "    FAILED 0" << std::endl;
        std::cout << "    " << expected << " != " << obtained << std::endl;
    }
}

void valid_price_at_t_test() {
    std::cout << "test_MonteCarlo test prix en t=2 - Basket" << std::endl;
    //init option
    double maturity = 6;
    int nbTimeSteps = 4;
    int size = 1;
    double strike = 10;
    PnlVect *lambda = pnl_vect_create_from_scalar(size, 1);
    OptionBasket opt(maturity, nbTimeSteps, size, strike, lambda);
    
    //Matrice jusqu'a t connu : t=2
    PnlMat *past = pnl_mat_create(3,1);
    past->array[0,0] = 100;
    past->array[0,1] = 92.12;
    past->array[0,2] = 87.93;

    //Calcul des valeurs attendues

    //init modèle
    BlackScholesModel* bs = new BlackScholesModel();
    bs->size_ = 1;
    bs->r_ = 0.1;
    bs->rho_ = 0.3;
    bs->sigma_ = pnl_vect_create_from_scalar(bs->size_, 0.2);
    bs->spot_ = pnl_vect_create_from_scalar(bs->size_, 100);
    bs->mat_cholesky = pnl_mat_create_from_scalar(bs->size_, bs->size_, bs->rho_);
    pnl_mat_set_diag(bs->mat_cholesky, 1, 0);
    pnl_mat_chol(bs->mat_cholesky);
    
    MonteCarlo *mc = new MonteCarlo();
    mc->fdStep_ = 1;
    mc->mod_ = bs;
    mc->nbSamples_ = 10;
    mc->opt_ = &opt;
    
    PnlVect * G = pnl_vect_create_from_scalar(bs->size_,0.25);
    
    PnlMat *path = pnl_mat_create(opt.nbTimeSteps_+1,bs->size_);
    double expected=0;
    
    for (int i=0; i<mc->nbSamples_; i++) {
        
        bs->asset(path, 2, opt.T_, opt.nbTimeSteps_, G, past);
        
        expected += opt.payoff(path);
    }
    expected *= exp(-bs->r_*(opt.T_ - 2)) / mc->nbSamples_;
    
    double obtained=0;
    double ic=0;
    mc->price(past, 2, obtained, ic, G);
    
    if (expected == obtained) {
        std::cout << "    TEST MonteCarlo 2 OK : "<< expected << " == " << obtained << std::endl;
    } else {
        std::cout << "    FAILED 1" << std::endl;
        std::cout << "    " << expected << " != " << obtained << std::endl;
    }
}

void valid_price_at_t_test2() {
    std::cout << "test_MonteCarlo test prix avec 2 titres - Performance" << std::endl;
    //init option
    double maturity = 10;
    int nbTimeSteps = 5;
    int size = 2;
    PnlVect *lambda = pnl_vect_create_from_scalar(size, 2);
    OptionPerformance opt(maturity, nbTimeSteps, size, lambda);
    
    //Matrice jusqu'a t connu (t=4) 
    PnlMat *past = pnl_mat_create(3,2);
    pnl_mat_set(past,0,0,100);
    pnl_mat_set(past,1,0,92.12);
    pnl_mat_set(past,2,0,87.93);
    pnl_mat_set(past,0,1,100);
    pnl_mat_set(past,1,1,94.12);
    pnl_mat_set(past,2,1,85.93);
    
    //init modèle
    BlackScholesModel* bs = new BlackScholesModel();
    bs->size_ = 2;
    bs->r_ = 0.1;
    bs->rho_ = 0.4;
    bs->sigma_ = pnl_vect_create_from_scalar(bs->size_, 0.2);
    bs->spot_ = pnl_vect_create_from_scalar(bs->size_, 100);
    bs->mat_cholesky = pnl_mat_create_from_scalar(bs->size_, bs->size_, bs->rho_);
    pnl_mat_set_diag(bs->mat_cholesky, 1, 0);
    pnl_mat_chol(bs->mat_cholesky);
    
    MonteCarlo *mc = new MonteCarlo();
    mc->fdStep_ = 1;
    mc->mod_ = bs;
    mc->nbSamples_ = 10;
    mc->opt_ = &opt;
    
    PnlVect * G = pnl_vect_create_from_scalar(bs->size_,0.25);
    
    PnlMat *path = pnl_mat_create(opt.nbTimeSteps_+1,bs->size_);
    double expected=0;
    double t=4;
    
    for (int i=0; i<mc->nbSamples_; i++) {      
        bs->asset(path, t, opt.T_, opt.nbTimeSteps_, G, past);
        
        expected += opt.payoff(path);
    }
    expected *= exp(-bs->r_*(opt.T_ - t)) / mc->nbSamples_;
    
    double obtained=0;
    double ic=0;
    mc->price(past, t, obtained, ic, G);
    
    if (expected == obtained) {
        std::cout << "    TEST price avec 2 titres OK : "<< expected << " == " << obtained << std::endl;
    } else {
        std::cout << "    FAILED 2" << std::endl;
        std::cout << "    " << expected << " != " << obtained << std::endl;
    }
}

void valid_price_at_t_test3() {
    std::cout << "test_MonteCarlo test prix avec 2 titres - Performance" << std::endl;
    //init option
    double maturity = 5;
    int nbTimeSteps = 5;
    int size = 2;
    PnlVect *lambda = pnl_vect_create_from_scalar(size, 2);
    OptionPerformance opt(maturity, nbTimeSteps, size, lambda);
    
    //Matrice jusqu'a t connu (t=1.5) 
    PnlMat *past = pnl_mat_create(3,2);
    pnl_mat_set(past,0,0,100);
    pnl_mat_set(past,1,0,92.12);
    pnl_mat_set(past,2,0,87.93);
    pnl_mat_set(past,0,1,100);
    pnl_mat_set(past,1,1,94.12);
    pnl_mat_set(past,2,1,85.93);
    
    //init modèle
    BlackScholesModel* bs = new BlackScholesModel();
    bs->size_ = 2;
    bs->r_ = 0.1;
    bs->rho_ = 0.5;
    bs->sigma_ = pnl_vect_create_from_scalar(bs->size_, 0.3);
    bs->spot_ = pnl_vect_create_from_scalar(bs->size_, 100);
    bs->mat_cholesky = pnl_mat_create_from_scalar(bs->size_, bs->size_, bs->rho_);
    pnl_mat_set_diag(bs->mat_cholesky, 1, 0);
    pnl_mat_chol(bs->mat_cholesky);
    
    MonteCarlo *mc = new MonteCarlo();
    mc->fdStep_ = 1;
    mc->mod_ = bs;
    mc->nbSamples_ = 10;
    mc->opt_ = &opt;
    
    PnlVect * G = pnl_vect_create_from_scalar(bs->size_,0.35);
    
    PnlMat *path = pnl_mat_create(opt.nbTimeSteps_+1,bs->size_);
    double expected=0;
    double t=1.5;
    
    for (int i=0; i<mc->nbSamples_; i++) {      
        bs->asset(path, t, opt.T_, opt.nbTimeSteps_, G, past);
        
        expected += opt.payoff(path);
    }
    expected *= exp(-bs->r_*(opt.T_ - t)) / mc->nbSamples_;
    
    double obtained=0;
    double ic=0;
    mc->price(past, t, obtained, ic, G);
    
    if (expected == obtained) {
        std::cout << "    TEST price avec 2 titres, t non multiple de T/N OK : "<< expected << " == " << obtained << std::endl;
    } else {
        std::cout << "    FAILED 3" << std::endl;
        std::cout << "    " << expected << " != " << obtained << std::endl;
    }
}


int main(int argc, char** argv) {
    std::cout << "%TEST_STARTED% test1 (test_MonteCarlo)" << std::endl;
    valid_price_at_0_test();
    
    valid_price_at_t_test();
    
    valid_price_at_t_test2();
    
    valid_price_at_t_test3();
    std::cout << "%TEST_FINISHED% time=0 test1 (test_MonteCarlo)" << std::endl;
    return (EXIT_SUCCESS);
}

