/* 
 * File:   test_basket.cpp
 * Author: gononl
 *
 * Created on September 16, 2016, 10:04 AM
 */

#include <iostream>
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"

#include "../src/OptionBasket.hpp"

using namespace std;

/*
 * Test de la fonction de calcul du payoff
 */
int main(int argc, char** argv) {    
    //////////////////////////
    ////// Option 1 //////////
    //////////////////////////
    double maturity = 1;
    int nbTimeSteps = 1;
    int size = 1;
    double strike = 10;
    PnlVect *lambda = pnl_vect_create_from_scalar(size, 1);
    // Option Basket avec 1 actif sous-jacent avec un cours constant à 20 sur 2 dates
    OptionBasket basket1(maturity, nbTimeSteps, size, strike, lambda);
    // Création du path représentant le cours du sous-jacent
    PnlMat *path = pnl_mat_create_from_scalar(nbTimeSteps + 1, size, 20);
    
    //Calcul du payoff de l'option
    double payoff = basket1.payoff(path);
    if (payoff != 10) {
        cout << "Payoff Basket 1 incorrect" << endl;
    } else {
        cout << "Payoff Basket 1 correct" << endl;
    }
    
    //////////////////////////
    ////// Option 2 //////////
    //////////////////////////
    maturity = 1;
    nbTimeSteps = 1;
    size = 2;
    strike = 10;
    lambda = pnl_vect_create_from_scalar(size, 0.7);
    pnl_vect_set(lambda, 1, 0.3);
    // Option Basket avec 2 actifs sous-jacent sur 2 dates
    // Un cours constant a 15 et un autre a 25
    OptionBasket basket2(maturity, nbTimeSteps, size, strike, lambda);
    // Création du path représentant le cours du sous-jacent
    path = pnl_mat_create_from_scalar(nbTimeSteps + 1, size, 15);
    pnl_mat_set_col(path, pnl_vect_create_from_scalar(nbTimeSteps + 1, 25), size - 1);

    // Calcul du payoff de l'option
    payoff = basket2.payoff(path);
    if (payoff != 8) {
        cout << "Payoff Basket 2 incorrect" << endl;
    } else {
        cout << "Payoff Basket 2 correct" << endl;
    }
    
    //////////////////////////
    ////// Option 3 //////////
    //////////////////////////
    maturity = 1;
    nbTimeSteps = 1;
    size = 2;
    strike = 20;
    lambda = pnl_vect_create_from_scalar(size, 0.7);
    pnl_vect_set(lambda, 1, 0.3);
    // Option Basket avec 2 actifs sous-jacent sur 2 dates avec des cours constants a 10
    OptionBasket basket3(maturity, nbTimeSteps, size, strike, lambda);
    // Création du path représentant le cours du sous-jacent
    path = pnl_mat_create_from_scalar(nbTimeSteps + 1, size, 10);

    // Calcul du payoff de l'option
    payoff = basket3.payoff(path);
    if (payoff != 0) {
        cout << "Payoff Basket 3 incorrect" << endl;
    } else {
        cout << "Payoff Basket 3 correct" << endl;
    }
    
    return 0;
}

