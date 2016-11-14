/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main_t.cpp
 * Author: paviotch
 *
 * Created on September 20, 2016, 8:51 AM
 */


#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <cstdlib>
#include "parser.hpp"
#include "BlackScholesModel.hpp"
#include "Option.hpp"
#include "OptionBasket.hpp"
#include "MonteCarlo.hpp"
#include <iostream>
#include <string>
#include <ctime>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    MonteCarlo *monte_carlo = new MonteCarlo();
    monte_carlo->mod_ = new BlackScholesModel();
    //monte_carlo->opt_ = new OptionBasket();
    double maturity;
    int nbTimeSteps;
    double strike;

    char *infile = argv[1];
    Param *P = new Parser(infile);
    P->extract("option size", monte_carlo->mod_->size_);
    P->extract("spot", monte_carlo->mod_->spot_, monte_carlo->mod_->size_);
    P->extract("volatility", monte_carlo->mod_->sigma_, monte_carlo->mod_->size_);
    P->extract("interest rate", monte_carlo->mod_->r_);
    P->extract("correlation", monte_carlo->mod_->rho_);
        
    P->extract("maturity", maturity);
    P->extract("TimeStep Number", nbTimeSteps);
    
    P->extract("strike", strike);
    cout << "strike : " << strike << endl;
    
    P->extract("Sample Number", monte_carlo->nbSamples_);  
    
    PnlVect* lambda = pnl_vect_create(monte_carlo->mod_->size_);
    
    P->extract("payoff coefficients", lambda, monte_carlo->mod_->size_);
    monte_carlo->opt_ = new OptionBasket(maturity, nbTimeSteps, monte_carlo->mod_->size_, strike, lambda);
    
  
    double price=0;
    double ic=5;
    
    double step = monte_carlo->opt_->T_/nbTimeSteps;
   
    
     
    PnlMat *path = pnl_mat_new();
    double t = 1.5;
    monte_carlo->rng_= pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(monte_carlo->rng_, time(NULL));
    
    PnlMat *past = pnl_mat_create(floor(t/step)+2,monte_carlo->mod_->size_);
    
    /*std::srand(std::time(NULL)); 
    
    for(int i = 0; i< past->m; i++){
        for (int j = 0; j<past->n;j++){
            pnl_mat_set(past,i,j, (std::rand()%300));
        }
    }*/
    
    monte_carlo->mod_->asset(past,monte_carlo->opt_->T_,nbTimeSteps,monte_carlo->rng_);
    
    cout << "past : " << endl;
    pnl_mat_print(past);
    
    
    /*monte_carlo->mod_->asset(path,t,maturity,nbTimeSteps,rng,past);
    
    cout << "path : " << endl;
    pnl_mat_print(path);*/
    
    cout << "PRICE : " << endl;
    monte_carlo->price(past,t,price,ic);
    cout << "price : " << price << endl;
    cout << "ic : " << ic << endl;
    
    return 0;
}

