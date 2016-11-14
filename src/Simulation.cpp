
#include "MonteCarlo.hpp"
#include "BlackScholesModel.hpp"
#include "Simulation.hpp"
#include <iostream>
#include "parser.hpp"

using namespace std;

Simulation::Simulation() {
    rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng,time(NULL));
    nbTimeStepH = 60;
    monte_carlo = new MonteCarlo();
    
}

Simulation::Simulation(Param *P) {
    rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    P->extract("hedging dates number", nbTimeStepH);
    monte_carlo = new MonteCarlo(P);

}

void Simulation::simu_couverture(PnlVect *val_pf, double &erreur_couverture, PnlVect *price) {

    PnlVect * V = pnl_vect_create(nbTimeStepH + 1);

    double T = monte_carlo->opt_->T_;
    int N = monte_carlo->opt_->nbTimeSteps_;

    PnlMat *path = pnl_mat_create(nbTimeStepH + 1, monte_carlo->mod_->size_);
    int nb_asset = monte_carlo->mod_->size_;

    simul_market(nb_asset, path);


    //initialisation : calcul de V0
    double prix = 0;
    double ic = 0;
    PnlVect *delta = pnl_vect_create(monte_carlo->mod_->size_);
    PnlVect *delta_prec = pnl_vect_create(monte_carlo->mod_->size_);

    //Calcul de p0 dans prix
    monte_carlo->price(prix, ic);

    pnl_vect_set(price,0,prix);
    pnl_vect_set(val_pf,0,prix);
    
    //Calcul de delta0  
    PnlMat *past = pnl_mat_create(1, monte_carlo->mod_->size_);
    pnl_mat_set_row(past, monte_carlo->mod_->spot_, 0);

    monte_carlo->delta(past, 0, delta);

    pnl_vect_clone(delta_prec, delta);

    //Calcul du portefeuille de couverture en 0
    double V0 = prix - pnl_vect_scalar_prod(delta, monte_carlo->mod_->spot_);

    pnl_vect_set(V, 0, V0);

    PnlVect *cours_date = pnl_vect_new();
    PnlVect *copy_delta = pnl_vect_new();
    PnlVect *row_to_add = pnl_vect_new();

    for (int i = 1; i < nbTimeStepH + 1; i++) {
        if (monte_carlo->opt_->nbTimeSteps_ != 0) {

            pnl_mat_get_row(row_to_add, path, i);

            if (!((i * monte_carlo->opt_->nbTimeSteps_) % nbTimeStepH == 0)) {


                if ((((i - 1) * monte_carlo->opt_->nbTimeSteps_) % nbTimeStepH == 0)) {

                    pnl_mat_add_row(past, past->m, row_to_add);

                } else {
                    pnl_mat_set_row(past, row_to_add, past->m - 1);
                }

            } else {
                pnl_mat_set_row(past, row_to_add, past->m - 1);
            }

            //A décommenter si on veut remplir le vecteur de prix théoriques de l'option
            // monte_carlo->price(past, i*T/nbTimeStepH, prix, ic);
            // pnl_vect_set(price,i,prix);

            monte_carlo->delta(past, i * T / nbTimeStepH, delta);
            pnl_vect_clone(copy_delta, delta);
            pnl_vect_minus_vect(copy_delta, delta_prec);

            pnl_mat_get_row(cours_date, path, i);

            pnl_vect_set(V, i, pnl_vect_get(V, i - 1) * exp((monte_carlo->mod_->r_ * T) / nbTimeStepH) -
                    pnl_vect_scalar_prod(copy_delta, cours_date));

            pnl_vect_set(val_pf, i, pnl_vect_get(V, i) + pnl_vect_scalar_prod(delta, cours_date));

            pnl_vect_clone(delta_prec, delta);
        }
    }

    cout << "Delta Final " << endl;
    pnl_vect_print(delta);

    erreur_couverture = pnl_vect_get(V, nbTimeStepH) +
            pnl_vect_scalar_prod(delta, cours_date) - monte_carlo->opt_->payoff(path);
    
    pnl_vect_free(&V);
    pnl_mat_free(&path);
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_prec);
    pnl_mat_free(&past);
    pnl_vect_free(&cours_date);
    pnl_vect_free(&copy_delta);
    pnl_vect_free(&row_to_add);
}


void Simulation::simu_couverture(PnlVect *val_pf, double &erreur_couverture, PnlVect *price, PnlVect *G) {

    PnlVect * V = pnl_vect_create(nbTimeStepH + 1);

    double T = monte_carlo->opt_->T_;
    int N = monte_carlo->opt_->nbTimeSteps_;

    PnlMat *path = pnl_mat_create(nbTimeStepH + 1, monte_carlo->mod_->size_);
    int nb_asset = monte_carlo->mod_->size_;

    this->monte_carlo->mod_->simul_market(path,T,nbTimeStepH,G);


    //initialisation : calcul de V0
    double prix = 0;
    double ic = 0;
    PnlVect *delta = pnl_vect_create(monte_carlo->mod_->size_);
    PnlVect *delta_prec = pnl_vect_create(monte_carlo->mod_->size_);

    //Calcul de p0 dans prix
    monte_carlo->price(prix, ic);

    pnl_vect_set(price,0,prix);
    pnl_vect_set(val_pf,0,prix);
    
    //Calcul de delta0  
    PnlMat *past = pnl_mat_create(1, monte_carlo->mod_->size_);
    pnl_mat_set_row(past, monte_carlo->mod_->spot_, 0);

    monte_carlo->delta(past, 0, delta,G);

    pnl_vect_clone(delta_prec, delta);

    //Calcul du portefeuille de couverture en 0
    double V0 = prix - pnl_vect_scalar_prod(delta, monte_carlo->mod_->spot_);

    pnl_vect_set(V, 0, V0);

    PnlVect *cours_date = pnl_vect_new();
    PnlVect *copy_delta = pnl_vect_new();
    PnlVect *row_to_add = pnl_vect_new();

    for (int i = 1; i < nbTimeStepH + 1; i++) {
        if (monte_carlo->opt_->nbTimeSteps_ != 0) {

            pnl_mat_get_row(row_to_add, path, i);

            if (!((i * monte_carlo->opt_->nbTimeSteps_) % nbTimeStepH == 0)) {


                if ((((i - 1) * monte_carlo->opt_->nbTimeSteps_) % nbTimeStepH == 0)) {

                    pnl_mat_add_row(past, past->m, row_to_add);

                } else {
                    pnl_mat_set_row(past, row_to_add, past->m - 1);
                }

            } else {
                pnl_mat_set_row(past, row_to_add, past->m - 1);
            }

            //A décommenter si on veut remplir le vecteur de prix théoriques de l'option
            // monte_carlo->price(past, i*T/nbTimeStepH, prix, ic,G);
            // pnl_vect_set(price,i,prix);

            monte_carlo->delta(past, i * T / nbTimeStepH, delta,G);
            pnl_vect_clone(copy_delta, delta);
            pnl_vect_minus_vect(copy_delta, delta_prec);

            pnl_mat_get_row(cours_date, path, i);

            pnl_vect_set(V, i, pnl_vect_get(V, i - 1) * exp((monte_carlo->mod_->r_ * T) / nbTimeStepH) -
                    pnl_vect_scalar_prod(copy_delta, cours_date));

            pnl_vect_set(val_pf, i, pnl_vect_get(V, i) + pnl_vect_scalar_prod(delta, cours_date));

            pnl_vect_clone(delta_prec, delta);
        }
    }

    cout << "Delta Final " << endl;
    pnl_vect_print(delta);

    erreur_couverture = pnl_vect_get(V, nbTimeStepH) +
            pnl_vect_scalar_prod(delta, cours_date) - monte_carlo->opt_->payoff(path);
    
    pnl_vect_free(&V);
    pnl_mat_free(&path);
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_prec);
    pnl_mat_free(&past);
    pnl_vect_free(&cours_date);
    pnl_vect_free(&copy_delta);
    pnl_vect_free(&row_to_add);

}

void Simulation::simul_market(int nbAsset, PnlMat *market) {
    monte_carlo->mod_->simul_market(market, monte_carlo->opt_->T_, nbTimeStepH, rng);
}

//void Simulation::simul_market(int nbAssets, PnlMat *market)
//{
//    const char *marketFile = "market.dat";
//    PnlMat *market_from_file = pnl_mat_create_from_file(marketFile);
//    pnl_mat_extract_subblock(market, market_from_file, 0, market_from_file->m, 0, nbAssets);
//    pnl_mat_free(&market_from_file);
//}