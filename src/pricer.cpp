/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   pricer.cpp
 * Author: paviotch
 *
 * Created on September 26, 2016, 11:05 AM
 */

#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#include "Simulation.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) 
{
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	clock_t t;
	t = clock();
	const char* extension = "-c";

	if ((argc == 3) && (strcmp(argv[1], extension) == 0))
	{
		if (0 == rank)
		{
		char *infile = argv[2];
		Param *P = new Parser(infile);

		Simulation *sim = new Simulation(P);
		PnlVect * val_pf = pnl_vect_create(sim->nbTimeStepH + 1);
		PnlVect * price = pnl_vect_create(sim->nbTimeStepH + 1);
		double err = 0;

		sim->simu_couverture(val_pf, err, price);

		cout << "erreur P&L : " << err << endl;
		}
	}
	else if (argc == 2)
	{

		char *infile = argv[1];
		Param *P = new Parser(infile);

		Simulation *sim = new Simulation(P);

		if (0 == rank)
		{
			double prix = 0;
			double ic = 0;

			sim->monte_carlo->price_master(prix, ic);

			cout << "prix en 0 : " << prix << endl;
			cout << "largeur de l'intervalle de confiance en 0 pour le prix : " << ic << endl;

		}
		else
		{
			sim->monte_carlo->price_slave();
		}
		MPI_Barrier(MPI_COMM_WORLD);

		if (0 == rank)
		{
			PnlMat * past = pnl_mat_create(1, sim->monte_carlo->mod_->size_);
			pnl_mat_set_row(past, sim->monte_carlo->mod_->spot_, 0);
			PnlVect *delta = pnl_vect_create(sim->monte_carlo->mod_->size_);

			sim->monte_carlo->delta(past, 0, delta);
			cout << "delta en 0 : " << endl;
			pnl_vect_print(delta);

			t = clock() - t;
			cout << "Temps d'exécution du programme : " << 
				((float)t)/CLOCKS_PER_SEC << " secondes." << endl;
		}
	}
	else
	{
		cout << "Veuillez tapez une des commandes suivantes"
			" : ./pricer fichier ou ./pricer -c fichier" << endl;
	}

	MPI_Finalize();

	return 0;
}

