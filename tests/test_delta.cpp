#include "../src/MonteCarlo.hpp"
#include "../src/BlackScholesModel.hpp"
#include "../src/OptionBasket.hpp"
#include "../src/OptionAsiatique.hpp"
#include "../src/OptionPerformance.hpp"
#include "../src/parser.hpp"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <time.h>

/*
 * Simple C++ Test Suite
 */
using namespace std;

void test_delta() {
    
    MonteCarlo *mc = new MonteCarlo();
    mc->nbSamples_ = 1;
    
    mc->mod_->size_ = 2;
    mc->mod_->r_ = 0.04879;
    mc->mod_->rho_ = 0.3;
    mc->mod_->sigma_ = pnl_vect_create_from_scalar(mc->mod_->size_, 0.2);
    mc->mod_->spot_ = pnl_vect_create_from_scalar(mc->mod_->size_, 100);
    mc->mod_->mat_cholesky = pnl_mat_create_from_scalar(mc->mod_->size_, mc->mod_->size_, mc->mod_->rho_);
    pnl_mat_set_diag(mc->mod_->mat_cholesky, 1, 0);
    pnl_mat_chol(mc->mod_->mat_cholesky);

    int nbTimeSteps = 3;
    PnlMat *past = pnl_mat_create(nbTimeSteps,mc->mod_->size_);
    double T = 6;
    double strike = 30;
    PnlVect *lambdaOption = pnl_vect_create_from_scalar(mc->mod_->size_, 0.5);
    mc->opt_ = new OptionBasket(T, nbTimeSteps, mc->mod_->size_, strike, lambdaOption);
    PnlVect *G = pnl_vect_create_from_scalar(mc->mod_->size_, 1.5);

    mc->mod_->asset(past, T, nbTimeSteps, G);

    PnlMat *true_past = pnl_mat_create(3, mc->mod_->size_);
    pnl_mat_set_subblock(true_past, past, 0, 0);


    mc->fdStep_ = 2;
    double t = 2.4;
    PnlMat *path = pnl_mat_create(mc->opt_->nbTimeSteps_+1,mc->mod_->size_);

    mc->mod_->asset(path, t, T, nbTimeSteps, G, true_past);

    PnlMat *path_plus = pnl_mat_new();
    PnlMat *vect_plus = pnl_mat_create(2,mc->mod_->size_);

    PnlMat *path_moins = pnl_mat_new();
    PnlMat *vect_moins = pnl_mat_create(2,mc->mod_->size_);
    
    PnlVect *delta_test = pnl_vect_create(mc->mod_->size_);
    
    double payoff_plus = 0;
    double payoff_moins = 0;   
    
    double constante = exp(-mc->mod_->r_*(T-t))/(2*mc->nbSamples_*mc->fdStep_);

    for (int i = 0; i < mc->mod_->size_; i++) {

        pnl_mat_clone(path_plus, path);

        pnl_mat_extract_subblock(vect_plus,path,floor(t),path->m - floor(t), i, 1);
        
        
        pnl_mat_mult_scalar(vect_plus, 1 + mc->fdStep_);

        pnl_mat_set(path_plus, floor(t), i, pnl_mat_get(vect_plus, 0, 0));
        pnl_mat_set(path_plus, floor(t) +1, i, pnl_mat_get(vect_plus, 1,0));

        payoff_plus = mc->opt_->payoff(path_plus);
   
        pnl_mat_clone(path_moins, path);
        pnl_mat_extract_subblock(vect_moins,path,floor(t),path->m - floor(t), i, 1);
        pnl_mat_mult_scalar(vect_moins, 1 - mc->fdStep_);
        pnl_mat_set(path_moins, floor(t), i, pnl_mat_get(vect_moins, 0, 0));
        pnl_mat_set(path_moins, floor(t) + 1, i, pnl_mat_get(vect_moins, 1,0));
        
        payoff_moins = mc->opt_->payoff(path_moins);

        pnl_vect_set(delta_test,i,constante * (1/pnl_mat_get(true_past,true_past->m-1,i))*(payoff_plus - payoff_moins));
        
        
    }


    PnlVect *delta = pnl_vect_create(2);
    mc->delta(true_past, t, delta,G);
    
     
     cout << " delta : " << endl;
     pnl_vect_print(delta);
     cout << " delta test : " << endl;
     pnl_vect_print(delta_test);

}

void test2(char** argv) {
    std::cout << "test_delta test 2" << std::endl;
    char *infile = argv[1];
    Param *P = new Parser(infile);
    MonteCarlo *monte_carlo = new MonteCarlo(P);
 
    PnlVect *delta = pnl_vect_create_from_zero(monte_carlo->mod_->size_);

    double t = 0;
    PnlMat *past = pnl_mat_create(1,monte_carlo->mod_->size_);
  
    pnl_mat_set_row(past,monte_carlo->mod_->spot_,0);
    cout << "past : " << endl;
    pnl_mat_print(past);
    
//    cout << "appel delta :" << endl; 
    monte_carlo->delta(past,t,delta);
//    cout << "fin delta :" << endl; 
        
    cout << "DELTA FINAL " << endl;
    pnl_vect_print(delta);
    
    pnl_vect_free(&delta);
    pnl_mat_free(&past);
}

void test3(){

    MonteCarlo *monte_carlo = new MonteCarlo();
    double t = 4.3;     
   
    PnlMat *past = pnl_mat_create(floor(t)+2,monte_carlo->mod_->size_);
    PnlVect *G = pnl_vect_create_from_scalar(monte_carlo->mod_->size_,0.7);

    monte_carlo->mod_->simul_market(past,floor(t)+1,floor(t)+1,G);

    PnlVect *delta = pnl_vect_create(monte_carlo->mod_->size_);  
 
    monte_carlo->delta(past,t,delta,G);
    cout << "DELTA : " << endl;
    pnl_vect_print(delta);
}


int main(int argc, char** argv) {
    std::cout << "%TEST_STARTED% test1 (test_delta)" << std::endl;
    test_delta();
    std::cout << "%TEST_FINISHED% time=0 test1 (test_delta)" << std::endl;
    std::cout << "%TEST_STARTED% test2 (test_delta)" << std::endl;
    test2(argv);
    std::cout << "%TEST_FINISHED% time=0 test1 (test_delta)" << std::endl;
    //test3();
    return (EXIT_SUCCESS);
}

