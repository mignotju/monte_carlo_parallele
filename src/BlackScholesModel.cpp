#include "BlackScholesModel.hpp"
#include <math.h>
#include <iostream>
#include "parser.hpp"

using namespace std;

void BlackScholesModel::packingSizePnlVect(int &bufsize, PnlVect* V) {
	int count;
	MPI_Pack_size(1,MPI_INT, MPI_COMM_WORLD,&count);
	bufsize += count;
	MPI_Pack_size(V->size,MPI_DOUBLE,MPI_COMM_WORLD,&count);
	bufsize += count;
}

void BlackScholesModel::packingPnlVect(char* buf, int bufsize, int &pos, PnlVect* V) {
	MPI_Pack(&(V->size), 1, MPI_INT, buf, bufsize, &pos, MPI_COMM_WORLD);
	MPI_Pack(V->array, V->size, MPI_DOUBLE, buf, bufsize, &pos, MPI_COMM_WORLD);
}

void BlackScholesModel::unpackingPnlVect(char* buf, int bufsize, int &pos, PnlVect* &V)  {
	int n;
	MPI_Unpack(buf,bufsize, &pos,&n,1,MPI_INT,MPI_COMM_WORLD);
	V = pnl_vect_create_from_zero(n);
	MPI_Unpack(buf,bufsize, &pos,V->array,n,MPI_DOUBLE,MPI_COMM_WORLD);
}

// void BlackScholesModel::packingSizePnlMat(int &bufsize, int &pos, PnlMat* M) {
// 	int info, count;
// 	info=MPI_Pack_size(1,MPI_INT, MPI_COMM_WORLD,&count);
// 	// if (info) return (info);
// 	bufsize += count;
// 	info=MPI_Pack_size(1,MPI_INT, MPI_COMM_WORLD,&count);
// 	// if (info) return (info);
// 	bufsize += count;
// 	info=MPI_Pack_size(1,MPI_INT, MPI_COMM_WORLD,&count);
// 	// if (info) return (info);
// 	bufsize += count;
// 	info=MPI_Pack_size(M->mn,MPI_DOUBLE,MPI_COMM_WORLD,&count);
// 	// if (info) return (info);
// 	bufsize += count;
// }
//
// void BlackScholesModel::packingPnlMat(char* buf, int bufsize, int pos, PnlMat* M) {
// 	int info;
// 	info=MPI_Pack(&(M->n),1,MPI_INT,buf,bufsize,&pos,MPI_COMM_WORLD);
// 	// if (info) return info;
// 	info=MPI_Pack(&(M->m),1,MPI_INT,buf,bufsize,&pos,MPI_COMM_WORLD);
// 	// if (info) return info;
// 	info=MPI_Pack(&(M->mn),1,MPI_INT,buf,bufsize,&pos,MPI_COMM_WORLD);
// 	// if (info) return info;
// 	info=MPI_Pack(M->array,M->mn,MPI_DOUBLE,buf,bufsize,&pos,MPI_COMM_WORLD);
// }
//
// void BlackScholesModel::unpackingPnlMat(char * buf, int bufsize, int pos, PnlMat* M)  {
// 	int n,m, info;
// 	info=MPI_Unpack(buf,bufsize,&pos,&n,1,MPI_INT,MPI_COMM_WORLD);
// 	// if (info) return info;
// 	info=MPI_Unpack(buf,bufsize,&pos,&m,1,MPI_INT,MPI_COMM_WORLD);
// 	// if (info) return info;
// 	info=MPI_Unpack(buf,bufsize,&pos,&(M->mn),1,MPI_INT,MPI_COMM_WORLD);
// 	// if (info) return info;
// 	pnl_mat_resize(M,n,m);
// 	info=MPI_Unpack(buf,bufsize,&pos,M->array,M->mn,MPI_DOUBLE,MPI_COMM_WORLD);
// 	// if (info) return info;
// }

BlackScholesModel::BlackScholesModel() {
    size_ = 3;
    r_ = 0.04879;
    rho_ = 0.3;
    sigma_ = pnl_vect_create_from_scalar(size_, 0.2);
    spot_ = pnl_vect_create_from_scalar(size_, 100);
    trend = pnl_vect_create_from_scalar(size_, 0.04879);
    G = pnl_vect_new();

    //Calculation of the Cholesky matrix of the correlation matrix
    mat_cholesky = pnl_mat_create_from_scalar(size_, size_, rho_);
    pnl_mat_set_diag(mat_cholesky, 1, 0);
    pnl_mat_chol(mat_cholesky);

    clone_past_ = pnl_mat_new();
    subBlock_ = pnl_mat_new();

}

BlackScholesModel::BlackScholesModel(Param *P, bool parallel) {

    G = pnl_vect_new();
    P->extract("option size", size_);
    P->extract("spot", spot_, size_);
    P->extract("volatility", sigma_, size_);
    P->extract("interest rate", r_);
    P->extract("correlation", rho_);
    P->extract("trend", trend, size_);

    if (parallel) {

      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      /*Getting pack size*/
      int bufsize=0, pos=0, count;
      MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count);
  		bufsize += count;
  		MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count);
  		bufsize += count;
  		MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count);
  		bufsize += count;

      packingSizePnlVect(bufsize, sigma_);
  		packingSizePnlVect(bufsize, spot_);
  		packingSizePnlVect(bufsize, trend);

      /*Creating pack*/
  		char* buf = new char[bufsize];

  		MPI_Pack(&(size_), 1, MPI_INT, buf, bufsize,
  				&pos, MPI_COMM_WORLD);
  		MPI_Pack(&(r_), 1, MPI_DOUBLE, buf, bufsize,
  				&pos, MPI_COMM_WORLD);
  		MPI_Pack(&(rho_), 1, MPI_DOUBLE, buf, bufsize,
  				&pos, MPI_COMM_WORLD);

  		packingPnlVect(buf, bufsize, pos, sigma_);
  		packingPnlVect(buf, bufsize, pos, spot_);
  		packingPnlVect(buf, bufsize, pos, trend);

      for (int i = 1; i < size; i++) {
  	 		MPI_Send(buf, bufsize, MPI_PACKED, i, 0, MPI_COMM_WORLD);
  		}
      delete(buf);
    }

    //Calculation of the Cholesky matrix of the correlation matrix
    mat_cholesky = pnl_mat_create_from_scalar(size_, size_, rho_);
    pnl_mat_set_diag(mat_cholesky, 1, 0);
    pnl_mat_chol(mat_cholesky);

    clone_past_ = pnl_mat_new();
    subBlock_ = pnl_mat_new();
}

BlackScholesModel::BlackScholesModel(bool slave) {
  if (slave) {

    int bufsize;
  	MPI_Status status;

  	MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
  	MPI_Get_count(&status, MPI_PACKED, &bufsize);

  	char* buf = new char[bufsize];

  	MPI_Recv(buf, bufsize, MPI_PACKED, 0, 0, MPI_COMM_WORLD, NULL);

  	int pos = 0;

  	MPI_Unpack(buf, bufsize, &pos, &(size_), 1, MPI_INT,
  	 		MPI_COMM_WORLD);
  	MPI_Unpack(buf, bufsize, &pos, &(r_), 1, MPI_DOUBLE,
  			MPI_COMM_WORLD);
  	MPI_Unpack(buf, bufsize, &pos, &(rho_), 1, MPI_DOUBLE,
  			MPI_COMM_WORLD);
  	unpackingPnlVect(buf, bufsize, pos, sigma_);
		pnl_vect_print(sigma_);
  	unpackingPnlVect(buf, bufsize, pos, spot_);
  	unpackingPnlVect(buf, bufsize, pos, trend);


  	pnl_vect_print(spot_);
		pnl_vect_print(trend);

  	delete(buf);

    //Calculation of the Cholesky matrix of the correlation matrix
    mat_cholesky = pnl_mat_create_from_scalar(size_, size_, rho_);
    pnl_mat_set_diag(mat_cholesky, 1, 0);
    pnl_mat_chol(mat_cholesky);

    clone_past_ = pnl_mat_new();
    subBlock_ = pnl_mat_new();
  }
}

void BlackScholesModel::asset(PnlMat* path, double T, int nbTimeSteps, PnlRng* rng) {

    PnlVect *cours_date = pnl_vect_new();
    pnl_vect_clone(cours_date, spot_);

    PnlVect *mat_chol_row = pnl_vect_create(size_);

    for (int date = 0; date < nbTimeSteps + 1; date++) {

        if (date == 0) {
            pnl_mat_set_row(path, cours_date, date);

        } else {

            pnl_vect_rng_normal(G, size_, rng);

            for (int i = 0; i < size_; i++) {

                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);
                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((r_ - pow(pnl_vect_get(sigma_, i), 2) / 2)*(T / nbTimeSteps) +
                        pnl_vect_get(sigma_, i) * sqrt(T / nbTimeSteps) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));

            }

            pnl_mat_set_row(path, cours_date, date);

        }
    }

    pnl_vect_free(&cours_date);
    pnl_vect_free(&mat_chol_row);
}

void BlackScholesModel::asset(PnlMat *path, double t, double T, int nbTimeSteps,
        PnlRng *rng, const PnlMat *past) {

    double step = T / nbTimeSteps;


    //PnlMat * clone_past_ = pnl_mat_new();
    pnl_mat_clone(clone_past_, past);


    if (fmod(t, T / nbTimeSteps) != 0) {
        pnl_mat_del_row(clone_past_, clone_past_->m - 1);
    }

    pnl_mat_set_subblock(path, clone_past_, 0, 0);


    PnlVect *last_row = pnl_vect_create(past->n);
    pnl_mat_get_row(last_row, past, (past->m) - 1);


    PnlVect *cours_date = pnl_vect_new();
    pnl_vect_clone(cours_date, last_row);

    PnlVect *mat_chol_row = pnl_vect_new();

    for (int date = clone_past_->m; date < nbTimeSteps + 1; date++) {

        pnl_vect_rng_normal(G, size_, rng);

        if (date == clone_past_->m && fmod(t, step) != 0) {
            for (int i = 0; i < size_; i++) {

                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);

                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((r_ - pow(pnl_vect_get(sigma_, i), 2) / 2)*((past->m - 1) * step - t) +
                        pnl_vect_get(sigma_, i) * sqrt((past->m - 1) * step - t) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));

            }
        } else {
            for (int i = 0; i < size_; i++) {

                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);

                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((r_ - pow(pnl_vect_get(sigma_, i), 2) / 2)*(step) +
                        pnl_vect_get(sigma_, i) * sqrt(step) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));

            }
        }

        pnl_mat_set_row(path, cours_date, date);
    }

    //pnl_mat_free(&clone_past);
    pnl_vect_free(&last_row);
    pnl_vect_free(&cours_date);
    pnl_vect_free(&mat_chol_row);

}

void BlackScholesModel::shiftAsset(PnlMat* shift_path, const PnlMat* path, int d, double h, double t, double timestep) {

    int ind_t;
    if (fmod(t, timestep) != 0) {
        ind_t = floor(t / timestep) + 1;
    } else {
        ind_t = t / timestep;
    }

    //Extraction de la sous colonne à modifier

    pnl_mat_extract_subblock(subBlock_, path, ind_t, path->m - ind_t, d, 1);

    //shift : *(1+h)
    pnl_mat_mult_scalar(subBlock_, 1 + h);

    //Remplacement de la sous colonne modifiée
    pnl_mat_resize(shift_path,path->m,path->n);
    pnl_mat_clone(shift_path, path);
    pnl_mat_set_subblock(shift_path, subBlock_, ind_t, d);

}

void BlackScholesModel::simul_market(PnlMat *path, double T, int H, PnlRng *rng) {

    PnlVect *cours_date = pnl_vect_new();
    pnl_vect_clone(cours_date, spot_);

    PnlVect *mat_chol_row = pnl_vect_new();

    for (int date = 0; date < H + 1; date++) {

        if (date != 0) {
            pnl_vect_rng_normal(G, size_, rng);
            for (int i = 0; i < size_; i++) {
                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);

                //Formule de récurrence avec trend pour la proba historique
                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((pnl_vect_get(trend, i) - pow(pnl_vect_get(sigma_, i), 2) / 2)*(T / H) +
                        pnl_vect_get(sigma_, i) * sqrt(T / H) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));
            }
        }

        pnl_mat_set_row(path, cours_date, date);
    }

    //libération de la mémoire
    pnl_vect_free(&cours_date);
    pnl_vect_free(&mat_chol_row);
}

/*  -----------  Fonctions déterministes pour les tests -------------------   */

void BlackScholesModel::asset(PnlMat* path, double T, int nbTimeSteps, PnlVect* G) {

    PnlVect *cours_date = pnl_vect_new();
    pnl_vect_clone(cours_date, spot_);

    PnlVect *mat_chol_row = pnl_vect_new();

    for (int date = 0; date < nbTimeSteps + 1; date++) {
        if (date == 0) {
            pnl_mat_set_row(path, cours_date, date);
        } else {
            for (int i = 0; i < size_; i++) {
                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);
                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((r_ - pow(pnl_vect_get(sigma_, i), 2) / 2)*(T / nbTimeSteps) +
                        pnl_vect_get(sigma_, i) * sqrt(T / nbTimeSteps) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));

            }
            pnl_mat_set_row(path, cours_date, date);
        }

    }

}

void BlackScholesModel::asset(PnlMat *path, double t, double T, int nbTimeSteps,
        PnlVect *G, const PnlMat *past) {

    //longueur d'un pas
    cout << T << endl;
    cout << nbTimeSteps;
    double step = T / nbTimeSteps;
    pnl_mat_set_all(path, 0);

    PnlMat * clone_past = pnl_mat_new();
    pnl_mat_clone(clone_past, past);


    if (fmod(t, T / nbTimeSteps) != 0) {
        pnl_mat_del_row(clone_past, clone_past->m - 1);
    }

    pnl_mat_set_subblock(path, clone_past, 0, 0);


    PnlVect *last_row = pnl_vect_create(past->n);
    pnl_mat_get_row(last_row, past, (past->m) - 1);

    PnlVect *cours_date = pnl_vect_new();
    pnl_vect_clone(cours_date, last_row);


    PnlVect *mat_chol_row = pnl_vect_new();
    for (int date = clone_past->m; date < nbTimeSteps + 1; date++) {

        if (date == clone_past->m && fmod(t, step) != 0) {

            for (int i = 0; i < size_; i++) {

                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);

                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((r_ - pow(pnl_vect_get(sigma_, i), 2) / 2)*((past->m - 1) * step - t) +
                        pnl_vect_get(sigma_, i) * sqrt((past->m - 1) * step - t) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));

            }
        } else {

            for (int i = 0; i < size_; i++) {

                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);

                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((r_ - pow(pnl_vect_get(sigma_, i), 2) / 2)*(step) +
                        pnl_vect_get(sigma_, i) * sqrt(step) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));

            }
        }


        pnl_mat_set_row(path, cours_date, date);

    }


    pnl_mat_free(&clone_past);
    pnl_vect_free(&last_row);
    pnl_vect_free(&cours_date);

}

void BlackScholesModel::simul_market(PnlMat *path, double T, int H, PnlVect *G) {

    PnlVect *cours_date = pnl_vect_new();
    pnl_vect_clone(cours_date, spot_);

    PnlVect *mat_chol_row = pnl_vect_new();

    for (int date = 0; date < H + 1; date++) {

        if (date != 0) {

            for (int i = 0; i < size_; i++) {

                pnl_mat_get_row(mat_chol_row, mat_cholesky, i);

                //Formule de récurrence avec trend pour la proba historique
                pnl_vect_set(cours_date, i, pnl_vect_get(cours_date, i) *
                        exp((pnl_vect_get(trend, i) - pow(pnl_vect_get(sigma_, i), 2) / 2)*(T / H) +
                        pnl_vect_get(sigma_, i) * sqrt(T / H) *
                        pnl_vect_scalar_prod(mat_chol_row, G)));
            }
        }

        pnl_mat_set_row(path, cours_date, date);
    }

    //libération de la mémoire
    pnl_vect_free(&cours_date);
    pnl_vect_free(&mat_chol_row);

}

BlackScholesModel::~BlackScholesModel() {
    pnl_vect_free(&sigma_);
    pnl_vect_free(&spot_);
    pnl_vect_free(&trend);
    pnl_vect_free(&G);
    pnl_mat_free(&mat_cholesky);
}
