/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_asset.cpp
 * Author: paviotch
 *
 * Created on September 15, 2016, 11:28 AM
 */

#include <cstdlib>
#include <math.h>
#include <iostream>
#include "pnl/pnl_vector.h"

#include "../src/BlackScholesModel.hpp"

using namespace std;

/*
 * Test de la simulation réalisée par la fonction asset de BlackScholesModel
 */
int main(int argc, char** argv) {
    ///////////////////////////////////////////////////////////
    //////////////////// Tests à partir de 0 //////////////////
    ///////////////////////////////////////////////////////////
    ////////////////////////////////////
    /////// Test Uni-dimensionnel //////
    ////////////////////////////////////
    BlackScholesModel* bl = new BlackScholesModel();
    // On fixe les paramètres du modèle 
    bl->size_ = 1;
    bl->r_ = 0.1;
    bl->rho_ = 0.3;
    bl->sigma_ = pnl_vect_create_from_scalar(bl->size_, 0.2);
    bl->spot_ = pnl_vect_create_from_scalar(bl->size_, 100);
    bl->mat_cholesky = pnl_mat_create_from_scalar(bl->size_, bl->size_, bl->rho_);
    pnl_mat_set_diag(bl->mat_cholesky, 1, 0);
    pnl_mat_chol(bl->mat_cholesky);
    // Maturité
    double maturity = 1;
    // Nombre de dates
    int nbTimeSteps = 1;
    
    PnlMat *path = pnl_mat_create(nbTimeSteps+1,bl->size_);
   
    // On remplace provisoirement le RNG par un vecteur pré-déterminé afin de vérifier les calculs
    PnlVect *vect = pnl_vect_create_from_scalar(bl->size_, 0.5);
    bl->asset(path, maturity, nbTimeSteps, vect);
    
    // Affichage de la matrice des cours calculée
//    cout << "Matrice Path après calcul : " << endl;
//    pnl_mat_print(path);
    
    if (pnl_mat_get(path, 1, 0) != 100*exp(0.18)) {
        cout << "Simulation 1 incorrecte" << endl;
    } else {
        cout << "Simulation 1 correcte" << endl;
    }
    
    
    ////////////////////////////////////
    /////// Test Multi-dimensionnel ////
    ////////////////////////////////////
    bl = new BlackScholesModel();
    // On fixe les paramètres du modèle
    bl->size_ = 2;
    bl->r_ = 0.1;
    bl->rho_ = 0.3;
    bl->sigma_ = pnl_vect_create_from_scalar(bl->size_, 0.2);
    pnl_vect_set(bl->sigma_, 1, 0.3);
    bl->spot_ = pnl_vect_create_from_scalar(bl->size_, 100);
    pnl_vect_set(bl->spot_, 1, 50);
    bl->mat_cholesky = pnl_mat_create_from_scalar(bl->size_, bl->size_, bl->rho_);
    pnl_mat_set_diag(bl->mat_cholesky, 1, 0);
    pnl_mat_chol(bl->mat_cholesky);
    
    // Maturité
    maturity = 1;
    // Nombre de dates
    nbTimeSteps = 1;
    
    path = pnl_mat_create(nbTimeSteps+1,bl->size_);
    
    // On remplace provisoirement le RNG par un vecteur pré-déterminé afin de vérifier "manuellement" les calculs
    vect = pnl_vect_create_from_scalar(bl->size_, 0.5);
    pnl_vect_set(vect, 1, 0.7);
    bl->asset(path, maturity, nbTimeSteps, vect);
    
    
    double resultatAttendu = 50*exp(0.055+(0.3*((0.3*0.5)+(0.953939*0.7))));
    if (pnl_mat_get(path, 1, 0) != 100*exp(0.18) ||
        !((pnl_mat_get(path, 1, 1)/resultatAttendu < 1.01) && (pnl_mat_get(path, 1, 1)/resultatAttendu > 0.99) ) ) {
        cout << "Simulation 2 incorrecte" << endl;
    } else {
        cout << "Simulation 2 correcte" << endl;
    }
    
    ///////////////////////////////////////////////////////////
    //////////////////// Tests à partir de t //////////////////
    ///////////////////////////////////////////////////////////
    ////////////////////////////////////
    /////// Test Uni-dimensionnel //////
    ////////////////////////////////////
    // Test où le t tombe "pile" sur un pas de discrétisation
    bl = new BlackScholesModel();
    // On fixe les paramètres du modèle 
    bl->size_ = 1;
    bl->r_ = 0.1;
    bl->rho_ = 0;
    bl->sigma_ = pnl_vect_create_from_scalar(bl->size_, 0.3);
    bl->spot_ = pnl_vect_create_from_scalar(bl->size_, 100);
    bl->mat_cholesky = pnl_mat_create_from_scalar(bl->size_, bl->size_, bl->rho_);
    pnl_mat_set_diag(bl->mat_cholesky, 1, 0);
    pnl_mat_chol(bl->mat_cholesky);
    
    // Paramètres de la fonction asset à partir de t
    
    double t = 1;
    maturity = 2;
    nbTimeSteps = 2;
    path = pnl_mat_create(nbTimeSteps+1, bl->size_);
    PnlMat *past = pnl_mat_create_from_scalar(2, 1, 100);
    vect = pnl_vect_create_from_scalar(bl->size_, 0.5);
    
    bl->asset(path, t, maturity, nbTimeSteps, vect, past);
    
    PnlMat *verif = pnl_mat_create_from_scalar(3, 1, 100);
    pnl_mat_set(verif, 2, 0, 100*exp(0.205));

    
    if (pnl_mat_eq(path, verif)) {
        cout << "Simulation 3 correcte" << endl;
        
    } else {
        cout << "Simulation 3 incorrecte" << endl;
        cout << "matrice 1" << endl;
        pnl_mat_print(path);
        cout << "matrice 2" << endl;
        pnl_mat_print(verif);
        cout << "FIN" << endl;
    }
    
    ////////////////////////////////////////////////////////////
    // Test où le t ne coincide pas avec un pas de discrétisation
    bl = new BlackScholesModel();
    // On fixe les paramètres du modèle 
    bl->size_ = 1;
    bl->r_ = 0.1;
    bl->rho_ = 0;
    bl->sigma_ = pnl_vect_create_from_scalar(bl->size_, 0.2);
    bl->spot_ = pnl_vect_create_from_scalar(bl->size_, 100);
    bl->mat_cholesky = pnl_mat_create_from_scalar(bl->size_, bl->size_, bl->rho_);
    pnl_mat_set_diag(bl->mat_cholesky, 1, 0);
    pnl_mat_chol(bl->mat_cholesky);
    
    // Paramètres de la fonction asset à partir de t
    
    t = 6;
    maturity = 10;
    nbTimeSteps = 2;
    path = pnl_mat_create(nbTimeSteps+1, bl->size_);
    // Attention a la definition des tailles de ces arguments !
    past = pnl_mat_create_from_scalar(3, 1, 100);
    vect = pnl_vect_create_from_scalar(bl->size_, 0.5);
    
    bl->asset(path, t, maturity, nbTimeSteps, vect, past);
    
    verif = pnl_mat_create_from_scalar(3, 1, 100);
    pnl_mat_set(verif, 2, 0, 100*exp(0.52));
    
    if (pnl_mat_eq(path, verif)) {
        cout << "Simulation 4 correcte" << endl;
    } else {
        cout << "Simulation 4 incorrecte" << endl;
    }
    
    
    
    ////////////////////////////////////
    /////// Test Multi-dimensionnel ////
    ////////////////////////////////////
    // Test où t tombe "pile" sur un pas de discrétisation
    bl = new BlackScholesModel();
    // On fixe les paramètres du modèle
    bl->size_ = 2;
    bl->r_ = 0.1;
    bl->rho_ = 0;
    bl->sigma_ = pnl_vect_create_from_scalar(bl->size_, 0.2);
    pnl_vect_set(bl->sigma_, 1, 0.3);
    bl->spot_ = pnl_vect_create_from_scalar(bl->size_, 100);
    pnl_vect_set(bl->spot_, 1, 50);
    bl->mat_cholesky = pnl_mat_create_from_scalar(bl->size_, bl->size_, bl->rho_);
    pnl_mat_set_diag(bl->mat_cholesky, 1, 0);
    pnl_mat_chol(bl->mat_cholesky);
    
    // Maturité
    maturity = 2;
    t = 1;
    // Nombre de dates
    nbTimeSteps = 2;
    
    path = pnl_mat_create(nbTimeSteps+1, bl->size_);
    past = pnl_mat_create_from_scalar(2, 2, 100);
    pnl_mat_set_col(past, pnl_vect_create_from_scalar(nbTimeSteps + 1, 50), 1);

    // On remplace provisoirement le RNG par un vecteur pré-déterminé afin de vérifier "manuellement" les calculs
    vect = pnl_vect_create_from_scalar(bl->size_, 0.5);
    pnl_vect_set(vect, 1, 0.7);
    
    bl->asset(path, t, maturity, nbTimeSteps, vect, past);
    
    verif = pnl_mat_create_from_scalar(3, 2, 100);
    pnl_mat_set_col(verif, pnl_vect_create_from_scalar(path->m, 50), 1);
    pnl_mat_set(verif, 2, 0, 100*exp(0.18));
    pnl_mat_set(verif, 2, 1, 50*exp(0.265));
    
    if (pnl_mat_eq(path, verif)) {
        cout << "Simulation 5 correcte" << endl;
    } else {
        cout << "Simulation 5 incorrecte" << endl;
        cout << "matrice 1" << endl;
        pnl_mat_print(path);
        cout << "matrice 2" << endl;
        pnl_mat_print(verif);
        cout << "FIN" << endl;
    }
    
    return 0;
}

