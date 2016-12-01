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

MonteCarlo::MonteCarlo()
{
    fdStep_ = 0.01;
    mod_ = new BlackScholesModel();
    nbSamples_ = 500;
    opt_ = new OptionBasket();

    int seed = time(NULL);

    int count;

    PnlRng ** array_rng = pnl_rng_dcmt_create_array_id(0,omp_get_num_threads()-1, seed, &count);
    if (count != omp_get_num_threads()) {
      std::cout << "Nombre de générateurs créés incorrects !" << std::endl;
    }
    for (int i = 0; i < count; i++) {
      pnl_rng_sseed(array_rng[i], seed);
    }

    rng_ = array_rng[omp_get_thread_num()];
    
    shiftPlus_ = pnl_mat_new();
    shiftMoins_ = pnl_mat_new();
    path_ = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
}



MonteCarlo::MonteCarlo(Param *P)
{
    mod_ = new BlackScholesModel(P);
    P->extract("fd step", fdStep_);

    //Option
    double maturity = 0;
    int nbTimeSteps = 0;
    double strike = 0;
    string type = "";
    P->extract("maturity", maturity);
    P->extract("TimeStep Number", nbTimeSteps);
    P->extract("strike", strike);
    P->extract("option type", type);
    PnlVect* lambda = pnl_vect_create(mod_->size_);
    P->extract("payoff coefficients", lambda, mod_->size_);

    P->extract("Sample Number", nbSamples_);

    if (type.compare("asian") == 0) {
        opt_ = new OptionAsiatique(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (type.compare("basket") == 0) {
        opt_ = new OptionBasket(maturity, nbTimeSteps, mod_->size_, strike, lambda);
    } else if (type.compare("performance") == 0) {
        opt_ = new OptionPerformance(maturity, nbTimeSteps, mod_->size_, lambda);
    }

    int seed = time(NULL);

    int count;
    PnlRng ** array_rng = pnl_rng_dcmt_create_array_id(0,omp_get_num_threads()-1, seed, &count);
    if (count != omp_get_num_threads()) {
      std::cout << "Nombre de générateurs créés incorrects !" << std::endl;
    }
    for (int i = 0; i < count; i++) {
      pnl_rng_sseed(array_rng[i], seed);
    }

    rng_ = array_rng[omp_get_thread_num()];

    shiftPlus_ = pnl_mat_new();
    shiftMoins_ = pnl_mat_new();
    path_ = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
}



void MonteCarlo::price(double &prix, double &ic) {
    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff = 0;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);


    for (int i = 0; i < nbSamples_; i++) {
        pnl_mat_set_all(path, 0);
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);
    }

    pnl_mat_free(&path);

    tmp = sum;
    variance = getVariance(tmp, sum_square, 0);
    prix = getPrice(sum, 0);
    //    std::cout << "VARIANCE : " << variance << std::endl;
    ic = getIntervalleConfiance(variance);

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
	int samples = nbSamples_;
	if (slaves == rank)
	{
		samples += (nbSamples_%slaves);
	}

	double res[2] = {0, 0};
	double payoff;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

	for (int i = 0; i < samples; i++) {
        pnl_mat_set_all(path, 0);
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        payoff = opt_->payoff(path);
        res[0] += payoff;
        res[1] += pow(payoff, 2);
    }

    pnl_mat_free(&path);

	MPI_Send(res, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}



void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic)
{
    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff = 0;

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
    for (int i = 0; i < nbSamples_; i++) {
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
    for (int j = 0; j < nbSamples_; j++) {
        // Simulation du path
        if (t == 0) {
            this->mod_->asset(path_, this->opt_->T_, this->opt_->nbTimeSteps_, rng_);
        } else {
            this->mod_->asset(path_, t, this->opt_->T_, this->opt_->nbTimeSteps_, rng_, past);
        }

        // Shift_path
        for (int i = 0; i < this->mod_->size_; i++) {
            // Création des trajectoires shiftées
            this->mod_->shiftAsset(shiftPlus_, path_, i, h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            this->mod_->shiftAsset(shiftMoins_, path_, i, -h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            payoff = this->opt_->payoff(shiftPlus_) - this->opt_->payoff(shiftMoins_);
            prec = pnl_vect_get(delta, i);
            pnl_vect_set(delta, i, prec + payoff);

            //pour l'intervalle de confiance
            pnl_vect_set(sum_square, i, pnl_vect_get(sum_square, i) + pow(payoff, 2));
        }
    }
    //Pour l'intervalle de confiance en 0
    if (t == 0) {
        for (int i = 0; i < this->mod_->size_; i++) {
            pnl_vect_set(ic, i, sqrt((pnl_vect_get(sum_square, i) / nbSamples_) - (pow(pnl_vect_get(delta, i), 2) / pow(nbSamples_, 2))));
        }
    }

    double coeff = exp(-mod_->r_ * (opt_->T_ - t)) / (2 * nbSamples_ * h);

    pnl_vect_mult_scalar(ic, 2 * coeff * nbSamples_ * 1.96 / sqrt(nbSamples_));
    pnl_vect_mult_scalar(delta, coeff);
    PnlVect *copy = pnl_vect_create_from_zero(past->n);
    pnl_mat_get_row(copy, past, past->m - 1);

    pnl_vect_div_vect_term(delta, copy);

    pnl_vect_div_vect_term(ic, copy);
    if (t == 0) {
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

/*  ----------- fonctions déterministes pour les tests ----------  */

void MonteCarlo::price(double &prix, double &ic, PnlVect *G)
{

    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff = 0;
    prix = 0;
    ic = 0;

    PnlMat *path = pnl_mat_create_from_ptr(1, mod_->size_, mod_->spot_->array);

    for (int i = 0; i < nbSamples_; i++) {

        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, G);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);
    }
    pnl_mat_free(&path);

    tmp = sum;
    tmp /= nbSamples_;
    tmp = pow(tmp, 2);
    sum_square /= nbSamples_;

    variance = exp(-2 * mod_->r_ * opt_->T_)*(sum_square - tmp);
    //std::cout << "VARIANCE : " << variance << std::endl;
    sum *= exp(-mod_->r_ * opt_->T_) / nbSamples_;
    prix = sum;
    ic = 2 * 1.96 * sqrt(variance / nbSamples_);
}



void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic, PnlVect *G)
{
    double sum = 0;
    double tmp = sum;
    double sum_square = 0;
    double variance = 0;
    double payoff;
    prix = 0;
    ic = 0;


    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    for (int i = 0; i < nbSamples_; i++) {
        mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, G, past);
        payoff = opt_->payoff(path);
        sum += payoff;
        sum_square += pow(payoff, 2);

    }

    pnl_mat_free(&path);
    tmp = sum;
    variance = getVariance(tmp, sum_square, t);
    prix = getPrice(sum, t);
    //    std::cout << "VARIANCE : " << variance << std::endl;
    ic = getIntervalleConfiance(variance);
}



void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *vect)
{

    double h = fdStep_;
    double payoff = 0;
    double prec = 0;
    PnlMat *path = pnl_mat_create_from_zero(this->opt_->nbTimeSteps_ + 1, this->mod_->size_);
    PnlMat *shiftPlus = pnl_mat_new();
    PnlMat *shiftMoins = pnl_mat_new();


    pnl_vect_set_all(delta, 0);
    // Moyenne des payoffs
    for (int j = 0; j < nbSamples_; j++) {
        // Simulation du path
        if (t == 0) {
            this->mod_->asset(path, this->opt_->T_, this->opt_->nbTimeSteps_, vect);
        } else {
            this->mod_->asset(path, t, this->opt_->T_, this->opt_->nbTimeSteps_, vect, past);
        }

        // Shift_path
        for (int i = 0; i < this->mod_->size_; i++) {
            // Création des trajectoires shiftées

            this->mod_->shiftAsset(shiftPlus, path, i, h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            this->mod_->shiftAsset(shiftMoins, path, i, -h, t, this->opt_->T_ / this->opt_->nbTimeSteps_);
            payoff = this->opt_->payoff(shiftPlus) - this->opt_->payoff(shiftMoins);
            prec = pnl_vect_get(delta, i);
            pnl_vect_set(delta, i, prec + payoff);
        }
    }

    double coeff = exp(-mod_->r_ * (opt_->T_ - t)) / (2 * nbSamples_ * h);
    pnl_vect_mult_scalar(delta, coeff);
    PnlVect *copy = pnl_vect_create_from_zero(past->n);
    pnl_mat_get_row(copy, past, past->m - 1);
    pnl_vect_div_vect_term(delta, copy);

    pnl_vect_free(&copy);
    pnl_mat_free(&shiftPlus);
    pnl_mat_free(&shiftMoins);
    pnl_mat_free(&path);
}

