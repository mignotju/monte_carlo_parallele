/* 
 * File:   test_performance.cpp
 * Author: lucas
 *
 * Created on September 17, 2016, 5:20 PM
 */

#include <cstdlib>
#include <iostream>

#include "../src/OptionPerformance.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    //////////////////////////
    ////// Option 1 //////////
    //////////////////////////
    double maturity = 1;
    int nbTimeSteps = 1;
    int size = 1;
    PnlVect *lambda = pnl_vect_create_from_scalar(size, 1);
    // Option Basket avec 1 actif sous-jacent avec un cours constant à 20 sur 2 dates
    OptionPerformance perf1(maturity, nbTimeSteps, size, lambda);
    // Création du path représentant le cours du sous-jacent
    PnlMat *path = pnl_mat_create_from_scalar(nbTimeSteps + 1, size, 20);
    
    //Calcul du payoff de l'option
    double payoff = perf1.payoff(path);
    if (payoff != 1) {
        cout << "Payoff Performance 1 incorrect" << endl;
    } else {
        cout << "Payoff Performance 1 correct" << endl;
    }
    
    //////////////////////////
    ////// Option 2 //////////
    //////////////////////////
    maturity = 1;
    nbTimeSteps = 1;
    size = 2;
    lambda = pnl_vect_create_from_scalar(size, 0.7);
    pnl_vect_set(lambda, 1, 0.3);
    // Option Basket avec 2 actifs sous-jacent sur 2 dates
    // Un cours constant a 15 et un autre a 25
    OptionPerformance perf2(maturity, nbTimeSteps, size, lambda);
    // Création du path représentant le cours du sous-jacent
    path = pnl_mat_create_from_scalar(nbTimeSteps + 1, size, 15);
    pnl_mat_set_col(path, pnl_vect_create_from_scalar(nbTimeSteps + 1, 25), size - 1);

    // Calcul du payoff de l'option
    payoff = perf2.payoff(path);
    if (payoff != 1) {
        cout << "Payoff Performance 2 incorrect" << endl;
    } else {
        cout << "Payoff Performance 2 correct" << endl;
    }
    
    //////////////////////////
    ////// Option 3 //////////
    //////////////////////////
    maturity = 1;
    nbTimeSteps = 1;
    size = 2;
    lambda = pnl_vect_create_from_scalar(size, 0.7);
    pnl_vect_set(lambda, 1, 0.3);

    OptionPerformance perf3(maturity, nbTimeSteps, size, lambda);
    // Création du path représentant le cours du sous-jacent
    path = pnl_mat_create_from_scalar(nbTimeSteps + 1, size, 10);
    pnl_mat_set_row(path, pnl_vect_create_from_scalar(size, 40), nbTimeSteps);
    // Calcul du payoff de l'option
    payoff = perf3.payoff(path);
    if (payoff != 4) {
        cout << "Payoff Performance 3 incorrect" << endl;
    } else {
        cout << "Payoff Performance 3 correct" << endl;
    }
    
    return 0;
}

