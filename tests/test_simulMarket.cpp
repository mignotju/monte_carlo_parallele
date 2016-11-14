/* 
 * File:   test_simulMarket.cpp
 * Author: girerda
 *
 * Created on September 22, 2016, 4:03 PM
 */

#include <stdlib.h>
#include <iostream>
#include "../src/BlackScholesModel.hpp"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"

/*
 * Simple C++ Test Suite
 */
using namespace std;

void test_simul_market_multidim_risqueneutre() {
    BlackScholesModel *bl = new BlackScholesModel();
    // On fixe les paramètres du modèle
    bl->size_ = 2;
    bl->r_ = 0.1;
    bl->rho_ = 0.3;
    bl->sigma_ = pnl_vect_create_from_scalar(bl->size_, 0.2);
    pnl_vect_set(bl->sigma_, 1, 0.3);
    bl->spot_ = pnl_vect_create_from_scalar(bl->size_, 100);
    pnl_vect_set(bl->spot_, 1, 50);
    bl->trend = pnl_vect_create_from_scalar(bl->size_,bl->r_);
    bl->mat_cholesky = pnl_mat_create_from_scalar(bl->size_, bl->size_, bl->rho_);
    pnl_mat_set_diag(bl->mat_cholesky, 1, 0);
    pnl_mat_chol(bl->mat_cholesky);
    
    // Maturité
    double maturity = 1;
    // Nombre de dates
    int nbTimeSteps = 1;
    
    PnlMat *path = pnl_mat_create(nbTimeSteps, bl->size_);
    PnlMat *path2 = pnl_mat_create(nbTimeSteps, bl->size_);
    
    // On remplace provisoirement le RNG par un vecteur pré-déterminé afin de vérifier "manuellement" les calculs
    PnlVect *vect = pnl_vect_create_from_scalar(bl->size_, 0.5);
    pnl_vect_set(vect, 1, 0.7);
    bl->asset(path, maturity, nbTimeSteps, vect);
    
    bl->simul_market(path2,maturity,nbTimeSteps,vect);

    if (pnl_mat_eq(path,path2)) {
        cout << "    TEST OK" << endl;
    } else {
        cout << "    Failed " << endl;
        cout << "path " << endl;
        pnl_mat_print(path);
        cout << "path2 " << endl;
        pnl_mat_print(path2);
    }
    
}

void test_simul_market_historique() {
    std::cout << "test_simul_market_historique" << std::endl;
    BlackScholesModel *bl = new BlackScholesModel();
    // On fixe les paramètres du modèle
    bl->size_ = 2;
    bl->r_ = 0.1;
    bl->rho_ = 0.3;
    bl->sigma_ = pnl_vect_create_from_scalar(bl->size_, 0.2);
    pnl_vect_set(bl->sigma_, 1, 0.3);
    bl->spot_ = pnl_vect_create_from_scalar(bl->size_, 100);
    pnl_vect_set(bl->spot_, 1, 50);
    bl->trend = pnl_vect_create_from_scalar(bl->size_,bl->r_);
    pnl_vect_set(bl->trend,1,0.085);
    
    // Maturité
    double maturity = 10;
    // Nombre de dates
    int nbTimeSteps = 11;
    
    PnlMat *path = pnl_mat_create(nbTimeSteps, bl->size_);
    PnlMat *path2 = pnl_mat_create(nbTimeSteps, bl->size_);
    
    // On remplace provisoirement le RNG par un vecteur pré-déterminé afin de vérifier "manuellement" les calculs
    PnlVect *vect = pnl_vect_create_from_scalar(bl->size_, 0.5);
    pnl_vect_set(vect, 1, 0.7);

    bl->simul_market(path2,maturity,nbTimeSteps,vect);
    
    PnlVect *expected = pnl_vect_new();
    PnlVect *obtained = pnl_vect_new();
    for (int i=0; i<bl->size_; i++) {
        bl->r_ = pnl_vect_get(bl->trend,i);
        bl->asset(path, maturity, nbTimeSteps, vect);
        pnl_mat_get_col(expected,path,i);
        pnl_mat_get_col(obtained,path2,i);
        
        if (pnl_vect_eq(expected,obtained)) {
            cout << "   TRAJECTOIRE " << i << " OK" << endl;
        } else {
        cout << "    FAILED" << endl;
            return;
        }
    }
        cout << "    TEST OK" << endl;
    
    
    
}

int main(int argc, char** argv) {
    
    std::cout << "%TEST_STARTED% test1 (test_simulMarket)" << std::endl;
    test_simul_market_multidim_risqueneutre();
    std::cout << "%TEST_FINISHED% time=0 test1 (test_simulMarket)" << std::endl;

    test_simul_market_historique();
    
    return (EXIT_SUCCESS);
}

