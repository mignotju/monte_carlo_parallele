/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main0.cpp
 * Author: mignotju
 *
 * Created on September 15, 2016, 11:50 AM
 */

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <cstdlib>
#include "parser.hpp"
#include "BlackScholesModel.hpp"
#include "Option.hpp"
#include "OptionBasket.hpp"
#include "MonteCarlo.hpp"
#include "OptionAsiatique.hpp"
#include "OptionPerformance.hpp"
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
    string type;

    char *infile = argv[1];
    Param *P = new Parser(infile);
    P->extract("option size", monte_carlo->mod_->size_);
    P->extract("spot", monte_carlo->mod_->spot_, monte_carlo->mod_->size_);
    P->extract("volatility", monte_carlo->mod_->sigma_, monte_carlo->mod_->size_);
    P->extract("interest rate", monte_carlo->mod_->r_);
     cout << "monte_carlo->mod_->r_ : " << monte_carlo->mod_->r_<< endl;
    P->extract("correlation", monte_carlo->mod_->rho_);
    cout << "monte_carlo->mod_->rho_ : " << monte_carlo->mod_->rho_<< endl;
    
    P->extract("maturity", maturity);
    P->extract("TimeStep Number", nbTimeSteps);
    
    P->extract("strike", strike);
    P->extract("option type",type);
    cout << "strike : " << strike << endl;
    
    P->extract("Sample Number", monte_carlo->nbSamples_);  
    cout << "NB SAMPLES : " << monte_carlo->nbSamples_ << endl;
    
    PnlVect* lambda = pnl_vect_create(monte_carlo->mod_->size_);
    
    P->extract("payoff coefficients", lambda, monte_carlo->mod_->size_);
    
    if (type.compare("asian") == 0){
            
        monte_carlo->opt_ = new OptionAsiatique(maturity, nbTimeSteps, monte_carlo->mod_->size_, strike, lambda);

    } else if (type.compare("basket") == 0){
        monte_carlo->opt_ = new OptionBasket(maturity, nbTimeSteps, monte_carlo->mod_->size_, strike, lambda);
    } else if (type.compare("performance") == 0){
        monte_carlo->opt_ = new OptionPerformance(maturity, nbTimeSteps, monte_carlo->mod_->size_, lambda);
    }
   // monte_carlo->opt_ = new OptionBasket(maturity, nbTimeSteps, monte_carlo->mod_->size_, strike, lambda);
    
  
    double price=0;
    double ic=5;
    monte_carlo->price(price,ic);
    
    cout << "price : " << price << endl;
    cout << "ic : " << ic << endl;
    

    return 0;
}

