#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "parser.hpp"

/*!
 * \file BlackScholesModel.hpp
 * \brief Modèle de Black-Scholes
 */
class BlackScholesModel
{
public:
    int size_; /// nombre d'actifs du modèle
    double r_; /// taux d'intérêt
    double rho_; /// paramètre de corrélation
    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *spot_; /// valeurs initiales du sous-jacent
    PnlVect *trend; /// pour la probabilité historique
    PnlVect *G; /// vecteur gaussien
    PnlMat *mat_cholesky; //matrice de cholesky
    
    PnlMat * clone_past_; //Variable temporaire utilisée dans asset
    PnlMat *subBlock_;
    
    /**
     * \brief Constructeur par defaut
     */
    BlackScholesModel();
    
    /**
     * \brief Constructeur permettant d'extraire des données d'un fichier
     * @param P  parseur d'extraction des données 
     */
    BlackScholesModel(Param *P);
    
    /**
     * \brief Destructeur
     */
    ~BlackScholesModel();
    
    /**
     * A utiliser pour les tests
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (N+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     * @param[in] G vecteur simulé à la main pour les tests
     */
    void asset(PnlMat *path, double T, int nbTimeSteps, PnlVect *G);
    
    /**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (N+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     */
    void asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);

    /**
     * Calcule une trajectoire du sous-jacent connaissant le
     * passé jusqu' à la date t
     *
     * @param[out] path  contient une trajectoire du sous-jacent
     * donnée jusqu'à l'instant T par la matrice past
     * @param[in] t date jusqu'à laquelle on connait la trajectoire.
     * t n'est pas forcément une date de discrétisation
     * @param[in] nbTimeSteps nombre de pas de constatation
     * @param[in] T date jusqu'à laquelle on simule la trajectoire
     * @param[in] past trajectoire réalisée jusqu'a la date t
     */
    void asset(PnlMat *path, double t, double T, int nbTimeSteps,
               PnlRng *rng, const PnlMat *past);
    
    /**
     * Calcule une trajectoire du sous-jacent connaissant le
     * passé jusqu' à la date t en supprimant l'aléatoire par le
     * parametre G
     * 
     * @param[out] path  contient une trajectoire du sous-jacent
     * donnée jusqu'à l'instant T par la matrice past
     * @param[in] t  date jusqu'à laquelle on connait la trajectoire.
     * t n'est pas forcément une date de discrétisation
     * @param[in] T  date jusqu'à laquelle on simule la trajectoire
     * @param[in] nbTimeSteps  nbTimeSteps nombre de pas de constatation
     * @param[in] G  vecteur fixe permettant les calculs "manuels"
     * @param[in] past  trajectoire réalisée jusqu'a la date t
     */
    void asset(PnlMat *path, double t, double T, int nbTimeSteps,
        PnlVect *G, const PnlMat *past);

    /**
     * Shift d'une trajectoire du sous-jacent
     *
     * @param[in]  path contient en input la trajectoire
     * du sous-jacent
     * @param[out] shift_path contient la trajectoire path
     * dont la composante d a été shiftée par (1+h)
     * à partir de la date t.
     * @param[in] t date à partir de laquelle on shift
     * @param[in] h pas de différences finies
     * @param[in] d indice du sous-jacent à shifter
     * @param[in] timestep pas de constatation du sous-jacent
     */
    void shiftAsset(PnlMat *shift_path, const PnlMat *path,
                    int d, double h, double t, double timestep);
    
    /**
     * Crée une simulation représentant un marché
     * 
     * @param[out] path contient les trajectoires des actifs présents sur le marché
     * @param[in] T  maturité
     * @param[in] H  nombre de rebalancement du portefeuille
     * @param[in] rng  générateur de nombres aléatoires
     */
    void simul_market(PnlMat *path, double T, int H, PnlRng *rng);
    
    /**
     * Crée une simulation représentant un marché pour les tests
     * i.e. sans aléatoire
     * 
     * @param[out] path  path contient les trajectoires des actifs présents sur le marché
     * @param[in] T  maturité
     * @param[in] H  nombre de rebalancement du portefeuille
     * @param[in] G  vecteur fixe représentant l'aléatoire pour nos tests
     */
    void simul_market(PnlMat *path, double T, int H, PnlVect *G);    
};


