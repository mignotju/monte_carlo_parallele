#include "../src/MonteCarlo.hpp"
#include "../src/BlackScholesModel.hpp"
#include "../src/OptionBasket.hpp"
#include "../src/OptionPerformance.hpp"
#include "../src/OptionAsiatique.hpp"
#include "../src/Simulation.hpp"
#include "../src/parser.hpp"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

/*
 * Simple C++ Test Suite
 */
using namespace std;


void test1() {
    std::cout << "test_couverture test 1" << std::endl;
    Simulation *sim = new Simulation();
    sim->monte_carlo = new MonteCarlo();

    PnlVect *val_pf = pnl_vect_create(sim->nbTimeStepH + 1);
    PnlVect *price = pnl_vect_create(sim->nbTimeStepH + 1);
    double err = 0;

    sim->simu_couverture(val_pf, err, price);
       
    cout << "prix pf : " << endl;
    pnl_vect_print(val_pf);
  
    cout << "erreur : " << err << endl;
}


void test2(char** argv) {

    std::cout << "test_couverture test 2" << std::endl;

    char *infile = argv[1];
    Param *P = new Parser(infile);
    
    Simulation *sim = new Simulation(P);

    PnlVect *val_pf = pnl_vect_create_from_zero(sim->nbTimeStepH + 1);
    PnlVect *price = pnl_vect_create_from_zero(sim->nbTimeStepH + 1);
    double err = 0;

    sim->simu_couverture(val_pf, err, price);

    cout << "prix pf : " << endl;
    pnl_vect_print(val_pf);

    cout << "erreur : " << err << endl;
    pnl_vect_free(&val_pf);
    pnl_vect_free(&price);

}

int main(int argc, char** argv) {
    std::cout << "%TEST_STARTED% test1 (test_couverture)" << std::endl;
    test1();
    test2(argv);
    std::cout << "%TEST_FINISHED% time=0 test1 (test_couverture)" << std::endl;
    return (EXIT_SUCCESS);
}

