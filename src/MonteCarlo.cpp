/*
 * File:   MonteCarlo.cpp
 * Author: cpaviot
 *
 * Created on 17 septembre 2016, 13:49
 */

#include "OptionBasket.hpp"
#include "OptionAsiatique.hpp"
#include "OptionPerformance.hpp"
#include "MonteCarlo.hpp"
#include "parser.hpp"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <mpi.h>
#include <omp.h>

using namespace std;

void MonteCarlo::packingSizePnlVect(int &bufsize, PnlVect* V) {
	int count;
	MPI_Pack_size(1,MPI_INT, MPI_COMM_WORLD,&count);
	bufsize += count;
	MPI_Pack_size(V->size,MPI_DOUBLE,MPI_COMM_WORLD,&count);
	bufsize += count;
}

void MonteCarlo::packingPnlVect(char* buf, int bufsize, int &pos, PnlVect* V) {
	MPI_Pack(&(V->size), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
	MPI_Pack(V->array, V->size, MPI_DOUBLE, buf, bufsize, &pos, MPI_COMM_WORLD);
}

void MonteCarlo::unpackingPnlVect(char* buf, int bufsize, int &pos, PnlVect* &V)  {
	int n;
	MPI_Unpack(buf,bufsize, &pos,&n,1,MPI_INT,MPI_COMM_WORLD);
	V = pnl_vect_create_from_zero(n);
	MPI_Unpack(buf,bufsize, &pos,V->array,n,MPI_DOUBLE,MPI_COMM_WORLD);
}


// MonteCarlo::MonteCarlo(Param *P, bool parallel)
// {
// 	if (parallel) {
// 		int rank,size;
// 		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 		MPI_Comm_size(MPI_COMM_WORLD, &size);
//
// 		if (0 == rank) {
// 			mod_ = new BlackScholesModel(P, true);
//
// 			P->extract("fd step", fdStep_);
//
// 			//Option
// 			double maturity = 0;
// 			int nbTimeSteps = 0;
// 			double strike = 0;
// 			string type = "";
// 			P->extract("maturity", maturity);
// 			P->extract("TimeStep Number", nbTimeSteps);
// 			P->extract("strike", strike);
// 			P->extract("option type", type);
// 			PnlVect* lambda = pnl_vect_create(mod_->size_);
// 			P->extract("payoff coefficients", lambda, mod_->size_);
//
// 			P->extract("Sample Number", nbSamples_);
//
// 			/*Getting pack size*/
// 			int bufsize=0, pos=0, count;
// 			MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count);
// 			bufsize += count;
// 			MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 			bufsize += count;
// 			MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count);
// 			bufsize += count;
// 			MPI_Pack_size(sizeof(type), MPI_BYTE, MPI_COMM_WORLD, &count);
// 			bufsize += count;
// 			packingSizePnlVect(bufsize, lambda);
// 			MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 			bufsize += count;
//
// 			/*Creating pack*/
// 			char* buf = new char[bufsize];
//
// 			MPI_Pack(&(maturity), 1, MPI_DOUBLE, buf, bufsize,
// 					&pos, MPI_COMM_WORLD);
// 			MPI_Pack(&(nbTimeSteps), 1, MPI_INT, buf, bufsize,
// 					&pos, MPI_COMM_WORLD);
// 			MPI_Pack(&(strike), 1, MPI_DOUBLE, buf, bufsize,
// 					&pos, MPI_COMM_WORLD);
// 			MPI_Pack(&(type), sizeof(type), MPI_BYTE, buf, bufsize,
// 					&pos, MPI_COMM_WORLD);
// 			packingPnlVect(buf, bufsize, pos, lambda);
// 			MPI_Pack(&(nbSamples_), 1, MPI_INT, buf, bufsize,
// 					&pos, MPI_COMM_WORLD);
//
//
// 			for (int i = 1; i < size; i++) {
// 				MPI_Send(buf, bufsize, MPI_PACKED, i, 0, MPI_COMM_WORLD);
// 			}
// 			delete(buf);
//
//
// 			if (type.compare("asian") == 0)
// 			{
// 				opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
// 			}
// 			else if (type.compare("basket") == 0)
// 			{
// 				opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
// 			}
// 			else if (type.compare("performance") == 0)
// 			{
// 				opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
// 			}
//
// 		} else {
// 			mod_ = new BlackScholesModel(true);
//
// 			int bufsize;
// 	  	MPI_Status status;
//
// 	  	MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
// 	  	MPI_Get_count(&status, MPI_PACKED, &bufsize);
//
// 	  	char* buf = new char[bufsize];
//
// 	  	MPI_Recv(buf, bufsize, MPI_PACKED, 0, 0, MPI_COMM_WORLD, NULL);
//
// 	  	int pos = 0;
// 			double maturity = 0;
// 			int nbTimeSteps = 0;
// 			double strike = 0;
// 			string type = "";
// 			PnlVect* lambda = pnl_vect_create(mod_->size_);
//
// 	  	MPI_Unpack(buf, bufsize, &pos, &(maturity), 1, MPI_DOUBLE,
// 	  	 		MPI_COMM_WORLD);
// 	  	MPI_Unpack(buf, bufsize, &pos, &(nbTimeSteps), 1, MPI_INT,
// 	  			MPI_COMM_WORLD);
// 	  	MPI_Unpack(buf, bufsize, &pos, &(strike), 1, MPI_DOUBLE,
// 	  			MPI_COMM_WORLD);
// 	  	unpackingPnlVect(buf, bufsize, pos, lambda);
//
// 			delete(buf);
//
// 			if (type.compare("asian") == 0)
// 			{
// 				opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
// 			}
// 			else if (type.compare("basket") == 0)
// 			{
// 				opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
// 			}
// 			else if (type.compare("performance") == 0)
// 			{
// 				opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
// 			}
// 		}
//
// 	} else {
// 		mod_ = new BlackScholesModel(P, false);
// 	}
//
//
// 	if (parallel)
// 	{
// 		int rank,size;
// 		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 		MPI_Comm_size(MPI_COMM_WORLD, &size);
// 		int seed = time(NULL);
//
// 		if (0 == rank)
// 		{
// 			int count;
// 			PnlRng ** array_rng = pnl_rng_dcmt_create_array_id(1, size, seed, &count);
// 			while (count != size)
// 			{
// 				cout << "Nombre de générateurs créés incorrects !" << endl;
// 				array_rng = pnl_rng_dcmt_create_array_id(1, size, seed, &count);
// 			}
//
// 			for (int i = 1; i < count; i++)
// 			{
// 				int count, bufsize=0, pos=0;
// 				PnlRng* my_rng = array_rng[i];
// 				pnl_rng_sseed(my_rng, seed);
//
// 				/*Getting pack size*/
// 				MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 				bufsize += count;
// 				MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 				bufsize += count;
// 				MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 				bufsize += count;
// 				MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 				bufsize += count;
// 				MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 				bufsize += count;
// 				MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count);
// 				bufsize += count;
// 				MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
// 				bufsize += count;
// 				MPI_Pack_size(my_rng->size_state, MPI_BYTE, MPI_COMM_WORLD, &count);
// 				bufsize += count;
//
// 				/*Creating pack*/
// 				char* buf = new char[bufsize];
//
// 				MPI_Pack(&(my_rng->type), 1, MPI_INT, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
// 				MPI_Pack(&(my_rng->rand_or_quasi), 1, MPI_INT, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
// 				MPI_Pack(&(my_rng->dimension), 1, MPI_INT, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
// 				MPI_Pack(&(my_rng->counter), 1, MPI_INT, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
// 				MPI_Pack(&(my_rng->has_gauss), 1, MPI_INT, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
// 				MPI_Pack(&(my_rng->gauss), 1, MPI_DOUBLE, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
// 				MPI_Pack(&(my_rng->size_state), 1, MPI_INT, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
// 				MPI_Pack(my_rng->state, my_rng->size_state, MPI_BYTE, buf, bufsize,
// 						&pos, MPI_COMM_WORLD);
//
// 				MPI_Send(buf, bufsize, MPI_PACKED, i, 0, MPI_COMM_WORLD);
//
// 				delete(buf);
// 			}
//
// 			rng_ = array_rng[0];
// 			pnl_rng_sseed(rng_, seed);
// 		}
// 		else
// 		{
// 			int bufsize;
// 			MPI_Status status;
//
// 			MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
// 			MPI_Get_count(&status, MPI_PACKED, &bufsize);
//
// 			char* buf = new char[bufsize];
//
// 			MPI_Recv(buf, bufsize, MPI_PACKED, 0, 0, MPI_COMM_WORLD, NULL);
//
// 			int pos = 0;
// 			int type;
//
// 			MPI_Unpack(buf, bufsize, &pos, &(type), 1, MPI_INT, MPI_COMM_WORLD);
// 			rng_ = pnl_rng_create(type);
//
// 			MPI_Unpack(buf, bufsize, &pos, &(rng_->rand_or_quasi), 1, MPI_INT,
// 					MPI_COMM_WORLD);
// 			MPI_Unpack(buf, bufsize, &pos, &(rng_->dimension), 1, MPI_INT,
// 					MPI_COMM_WORLD);
// 			MPI_Unpack(buf, bufsize, &pos, &(rng_->counter), 1, MPI_INT,
// 					MPI_COMM_WORLD);
// 			MPI_Unpack(buf, bufsize, &pos, &(rng_->has_gauss), 1, MPI_INT,
// 					MPI_COMM_WORLD);
// 			MPI_Unpack(buf, bufsize, &pos, &(rng_->gauss), 1, MPI_DOUBLE,
// 					MPI_COMM_WORLD);
// 			MPI_Unpack(buf, bufsize, &pos, &(rng_->size_state), 1, MPI_INT,
// 					MPI_COMM_WORLD);
// 			MPI_Unpack(buf, bufsize, &pos, rng_->state, rng_->size_state, MPI_BYTE,
// 					MPI_COMM_WORLD);
//
// 			pnl_rng_sseed(rng_, seed);
// 			delete(buf);
// 		}
// 	}
// 	else
// 	{
// 		rng_ = pnl_rng_create(PNL_RNG_MERSENNE);
// 		pnl_rng_sseed(rng_, time(NULL));
// 	}
//
// 	shiftPlus_ = pnl_mat_new();
// 	shiftMoins_ = pnl_mat_new();
// 	path_ = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
// }
//

//--------------------------------------------------------------------------------------------------------------------------------------


/* Méthode sans pack/unpack, où chaque thread créé son générateur avec un id différent,
* donc les générateurs sont indépendants
* plus rapide a priori
*/
MonteCarlo::MonteCarlo(Param *P, bool parallel)
{
	if (parallel) {
		int rank,size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		if (0 == rank) {
			mod_ = new BlackScholesModel(P, true);

			P->extract("fd step", fdStep_);

			//Option
			double maturity = 0;
			int nbTimeSteps = 0;
			double strike = 0;
			string type = "";
			int type_length;
			P->extract("maturity", maturity);
			P->extract("TimeStep Number", nbTimeSteps);
			P->extract("strike", strike);
			P->extract("option type", type);
			type_length = type.size() + 1;
			PnlVect* lambda = pnl_vect_create(mod_->size_);
			P->extract("payoff coefficients", lambda, mod_->size_);

			P->extract("Sample Number", nbSamples_);

			/*Getting pack size*/
			int bufsize=0, pos=0, count;
			MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count);
			bufsize += count;
			MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
			bufsize += count;
			MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count);
			bufsize += count;
			MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
			bufsize += count;
			MPI_Pack_size(type_length, MPI_CHAR, MPI_COMM_WORLD, &count);
			bufsize += count;
			packingSizePnlVect(bufsize, lambda);
			MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
			bufsize += count;

			/*Creating pack*/
			char* buf = new char[bufsize];

			MPI_Pack(&(maturity), 1, MPI_DOUBLE, buf, bufsize,
					&pos, MPI_COMM_WORLD);
			MPI_Pack(&(nbTimeSteps), 1, MPI_INT, buf, bufsize,
					&pos, MPI_COMM_WORLD);
			MPI_Pack(&(strike), 1, MPI_DOUBLE, buf, bufsize,
					&pos, MPI_COMM_WORLD);
			MPI_Pack(&(type_length), 1, MPI_INT, buf, bufsize,
					&pos, MPI_COMM_WORLD);
			MPI_Pack(&(type), type_length, MPI_CHAR, buf, bufsize,
					&pos, MPI_COMM_WORLD);
			packingPnlVect(buf, bufsize, pos, lambda);
			MPI_Pack(&(nbSamples_), 1, MPI_INT, buf, bufsize,
					&pos, MPI_COMM_WORLD);


			for (int i = 1; i < size; i++) {
				MPI_Send(buf, bufsize, MPI_PACKED, i, 0, MPI_COMM_WORLD);
			}
			delete(buf);


			if (type.compare("asian") == 0)
			{
				opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
			}
			else if (type.compare("basket") == 0)
			{
				opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
			}
			else if (type.compare("performance") == 0)
			{
				opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
			}

		} else {
			mod_ = new BlackScholesModel(true);

			int bufsize;
	  	MPI_Status status;

	  	MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
	  	MPI_Get_count(&status, MPI_PACKED, &bufsize);

	  	char* buf = new char[bufsize];

	  	MPI_Recv(buf, bufsize, MPI_PACKED, 0, 0, MPI_COMM_WORLD, NULL);

	  	int pos = 0;
			double maturity = 0;
			int nbTimeSteps = 0;
			double strike = 0;
			int type_length = 0;
			PnlVect* lambda = pnl_vect_create(mod_->size_);

	  	MPI_Unpack(buf, bufsize, &pos, &(maturity), 1, MPI_DOUBLE,
	  	 		MPI_COMM_WORLD);
			cout << "maturity : " << maturity << endl;
	  	MPI_Unpack(buf, bufsize, &pos, &(nbTimeSteps), 1, MPI_INT,
	  			MPI_COMM_WORLD);
					cout << "nbTimesSteps : " << nbTimeSteps << endl;

	  	MPI_Unpack(buf, bufsize, &pos, &(strike), 1, MPI_DOUBLE,
	  			MPI_COMM_WORLD);
					cout << "strike : " << strike << endl;
			MPI_Unpack(buf, bufsize, &pos, &(type_length), 1, MPI_INT,
			  	MPI_COMM_WORLD);
			cout << "type length : " << type_length << endl;
			//char type[type_length];
			char* type;
			type = (char*) malloc(type_length*sizeof(char));
			MPI_Unpack(buf, bufsize, &pos, type, type_length, MPI_CHAR,
			  	MPI_COMM_WORLD);
			for (int i = 0; i < type_length; i++) {
				cout << " " << type[i] << " ";
			}
			cout << "type : " << type << endl;
			// string type_s = string s(type);
	  	unpackingPnlVect(buf, bufsize, pos, lambda);
			pnl_vect_print(lambda);

			delete(buf);

			// if (type_s.compare("asian") == 0)
			// {
			// 	opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
			// }
			// else if (type_s.compare("basket") == 0)
			// {
			// 	opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
			// }
			// else if (type_s.compare("performance") == 0)
			// {
			// 	opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
			// }
		}

	} else {
		mod_ = new BlackScholesModel(P, false);
	}


	if (parallel)
	{
		int rank,size;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		int seed = time(NULL);

		if (0 == rank)
		{
			rng_ = pnl_rng_create(PNL_RNG_MERSENNE);
			pnl_rng_sseed(rng_, time(NULL));
		}
		else
		{
			rng_ = pnl_rng_dcmt_create_id(rank, seed);
			pnl_rng_sseed(rng_, seed);
		}
	}
	else
	{
		rng_ = pnl_rng_create(PNL_RNG_MERSENNE);
		pnl_rng_sseed(rng_, time(NULL));
	}

	shiftPlus_ = pnl_mat_new();
	shiftMoins_ = pnl_mat_new();
	path_ = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
}


void MonteCarlo::price(double &prix, double &ic)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (0 == rank)
	{
		price_master(prix, ic);
	}
	else
	{
		price_slave();
	}
}



void MonteCarlo::price_master(double &prix, double &ic)
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int slaves = size-1;
	double sum = 0;
	double sumSq = 0;
	double res[2];

	for (int i = 0; i < slaves; i++)
	{

		MPI_Recv(res, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
		sum += res[0];
		sumSq += res[1];
	}

	double variance = getVariance(sum, sumSq, 0);
	prix = getPrice(sum, 0);
	ic = getIntervalleConfiance(variance);

}



void MonteCarlo::price_slave()
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//on calcule le nombre de samples par slave
	int slaves = size - 1;
	int samples = nbSamples_/slaves;
	if (slaves == rank)
	{
		samples += (nbSamples_%slaves);
	}

	double res[2] = {0, 0};
	double payoff;

	PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

	for (int i = 0; i < samples; i++)
	{
		pnl_mat_set_all(path, 0);
		mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
		payoff = opt_->payoff(path);
		res[0] += payoff;
		res[1] += payoff*payoff;
	}

	pnl_mat_free(&path);

	MPI_Send(res, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

void MonteCarlo::price(double &prix, double &ic, double precision)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	nbSamples_ = 0;

	if (0 == rank)
	{
		price_master_precision(prix, ic, precision);
	}
	else
	{
		price_slave_precision();
	}
}


void MonteCarlo::price_master_precision(double &prix, double &ic, double precision)
{
	int size;
  int precision_reached = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int slaves = size-1;
  double sum = 0;
  double sumSq = 0;
	double res[2];

	 while (precision_reached != 1) {
		 nbSamples_ += 1000;
		for (int i = 0; i < slaves; i++)
		{
			MPI_Recv(res, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
			sum += res[0];
			sumSq += res[1];
		}

		double variance = getVariance(sum, sumSq, 0);
		prix = getPrice(sum, 0);
		ic = getIntervalleConfiance(variance);

		if (ic < precision) {
			precision_reached = 1;
		}

		for (int i = 1; i <= slaves; i++) {
			MPI_Send(&precision_reached, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	 }

}

void MonteCarlo::price_slave_precision()
{
	int rank, size;
  int precision_reached = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//on calcule le nombre de samples par slave
	int slaves = size - 1;
	int samples = 1000/slaves;
	if (slaves == rank)
	{
		samples += (1000%slaves);
	}

	while (precision_reached != 1) {

		double res[2] = {0, 0};
		double payoff;

		PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

		for (int i = 0; i < samples; i++)
		{
			pnl_mat_set_all(path, 0);
			mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
			payoff = opt_->payoff(path);
			res[0] += payoff;
			res[1] += payoff*payoff;
		}

		pnl_mat_free(&path);

		MPI_Send(res, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

		MPI_Recv(&precision_reached, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
	 }

}


void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic)
{
	double sum = 0;
	double tmp = sum;
	double sum_square = 0;
	double variance = 0;
	double payoff = 0;

	PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
	for (int i = 0; i < nbSamples_; i++)
	{
		mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);
		payoff = opt_->payoff(path);
		sum += payoff;
		sum_square += pow(payoff, 2);
	}

	pnl_mat_free(&path);

	tmp = sum;
	variance = getVariance(tmp, sum_square, t);
	prix = getPrice(sum, t);
	ic = getIntervalleConfiance(variance);
}



void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta)
{
	double h = fdStep_;
	double payoff = 0;
	double prec = 0;

	PnlVect *sum_square = pnl_vect_create_from_zero(delta->size);
	PnlVect *ic = pnl_vect_create_from_zero(delta->size);

	pnl_vect_set_all(delta, 0);

	// Moyenne des payoffs
	for (int j = 0; j < nbSamples_; j++)
	{
		// Simulation du path
		if (t == 0)
		{
			this->mod_->asset(path_, this->opt_->T_, this->opt_->nbTimeSteps_, rng_);
		}
		else
		{
			this->mod_->asset(path_, t, this->opt_->T_, this->opt_->nbTimeSteps_, rng_,
					past);
		}

		// Shift_path
		for (int i = 0; i < this->mod_->size_; i++)
		{
			// Création des trajectoires shiftées
			this->mod_->shiftAsset(shiftPlus_, path_, i, h, t,
					this->opt_->T_ / this->opt_->nbTimeSteps_);
			this->mod_->shiftAsset(shiftMoins_, path_, i, -h, t,
					this->opt_->T_ / this->opt_->nbTimeSteps_);
			payoff = this->opt_->payoff(shiftPlus_) - this->opt_->payoff(shiftMoins_);
			prec = pnl_vect_get(delta, i);
			pnl_vect_set(delta, i, prec + payoff);

			//pour l'intervalle de confiance
			pnl_vect_set(sum_square, i, pnl_vect_get(sum_square, i) + pow(payoff, 2));
		}
	}

	//Pour l'intervalle de confiance en 0
	if (t == 0)
	{
		for (int i = 0; i < this->mod_->size_; i++)
		{
			pnl_vect_set(ic, i,
					sqrt((pnl_vect_get(sum_square, i) / nbSamples_) -
						(pow(pnl_vect_get(delta, i), 2) / pow(nbSamples_, 2))));
		}
	}

	double coeff = exp(-mod_->r_ * (opt_->T_ - t)) / (2 * nbSamples_ * h);

	pnl_vect_mult_scalar(ic, 2 * coeff * nbSamples_ * 1.96 / sqrt(nbSamples_));
	pnl_vect_mult_scalar(delta, coeff);
	PnlVect *copy = pnl_vect_create_from_zero(past->n);
	pnl_mat_get_row(copy, past, past->m - 1);

	pnl_vect_div_vect_term(delta, copy);

	pnl_vect_div_vect_term(ic, copy);
	if (t == 0)
	{
		cout << "largeur des intervalles de confiance pour DELTA en t=0 : " << endl;
		pnl_vect_print(ic);
	}

	pnl_vect_free(&copy);
	pnl_vect_free(&sum_square);
	pnl_vect_free(&ic);

}



MonteCarlo::~MonteCarlo()
{
	pnl_rng_free(&rng_);
	pnl_mat_free(&shiftMoins_);
	pnl_mat_free(&shiftPlus_);
	pnl_mat_free(&path_);
}

/*  ----------- fonctions auxiliaires de factorisation du code ----------  */

double MonteCarlo::getVariance(double sum, double sum_square, double t)
{
	// cout << "sum : " << sum << endl;
	// cout << "sumSq : " << sum_square << endl;
	sum /= nbSamples_;
	sum = pow(sum, 2);

	sum_square /= nbSamples_;
	return exp(-2 * mod_->r_ * (opt_->T_ - t))*(sum_square - sum);
}



double MonteCarlo::getIntervalleConfiance(double variance)
{
	return 2 * 1.96 * sqrt(variance / nbSamples_);
}



double MonteCarlo::getPrice(double sum, double t)
{
	sum *= exp(-mod_->r_ * (opt_->T_ - t)) / nbSamples_;
	return sum;
}
