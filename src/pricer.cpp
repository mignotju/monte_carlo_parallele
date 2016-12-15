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

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	clock_t t;
	t = clock();
	const char* extension = "-c";

	// cas avec l'extension -c
	if ((argc == 3) && (strcmp(argv[1], extension) == 0))
	{
		if (0 == rank)
		{
			char *infile = argv[2];
			Param *P = new Parser(infile);

			Simulation *sim = new Simulation(P, false);
			PnlVect * val_pf = pnl_vect_create(sim->nbTimeStepH + 1);
			PnlVect * price = pnl_vect_create(sim->nbTimeStepH + 1);
			double err = 0;

			sim->simu_couverture(val_pf, err, price);

			cout << "erreur P&L : " << err << endl;
		}
	}

	// cas sans precision
	else if (argc == 2)
	{
		char *infile = argv[1];
		Param *P = new Parser(infile);

		Simulation *sim = new Simulation(P, true);

		double prix = 0;
		double ic = 0;

		sim->monte_carlo->price(prix, ic);
		MPI_Barrier(MPI_COMM_WORLD);

		if (0 == rank)
		{
			cout << "prix en 0 : " << prix << endl;
			cout << "largeur de l'intervalle de confiance en 0 pour le prix : "
				<< ic << endl;

			t = clock() - t;
			cout << "Temps d'exécution du programme : " <<
				((float)t)/CLOCKS_PER_SEC << " secondes." << endl;
		}
	}

	// cas avec precision
	else if (argc == 3)
	{
		char *infile = argv[1];
		Param *P = new Parser(infile);

		double precision = atof(argv[2]);

		Simulation *sim = new Simulation(P, true);

		double prix = 0;
		double ic = 0;
		bool precision_reached = false;

		sim->monte_carlo->price(prix, ic, precision);
		MPI_Barrier(MPI_COMM_WORLD);

		if (0 == rank)
		{
			cout << "prix en 0 : " << prix << endl;
			cout << "largeur de l'intervalle de confiance en 0 pour le prix : "
				<< ic << endl;
			cout << "nombre de tirages requis : " << sim->monte_carlo->nbSamples_ << endl;

			t = clock() - t;
			cout << "Temps d'exécution du programme : " <<
				((float)t)/CLOCKS_PER_SEC << " secondes." << endl;
		}
	}
	else
	{
		if (0 == rank)
		{
			cout << "Veuillez tapez une des commandes suivantes"
				<< " : ./pricer fichier ou ./pricer -c fichier" << endl;
		}
	}

	MPI_Finalize();

	return 0;
}
