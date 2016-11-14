/* 
 * File:   Simulation.hpp
 * Author: mignotju
 *
 * Created on September 21, 2016, 4:24 PM
 */
#pragma once

#include "MonteCarlo.hpp"
#include "pnl/pnl_random.h"
#include "parser.hpp"

/*!
 * \file Simulation.hpp
 * \brief Simulation 
 */
class Simulation
{
public:
    MonteCarlo* monte_carlo;
    int nbTimeStepH; /// nombre de dates de rebalancement
    PnlRng *rng;
    
    /**
     * \brief Constructeur par defaut
     */
    Simulation();
    
    /**
     * \brief Constructeur à partir d'un fichier
     * @param P  parseur pour l'extraction de données dans un fichier
     */
    Simulation(Param *P);
    
    /**
     * \brief Destructeur
     */
    ~Simulation();
    
    /**
     * Fonction de simulation de la couverture
     * 
     * @param[out] V  valeur du portefeuille de couverture
     * @param[out] erreur_couverture  erreur de couverture
     * @param[out] price  prix théorique de l'option
     */
    void simu_couverture(PnlVect *V, double &erreur_couverture, PnlVect *price);
    
    /**
     * Fonction de simulation de la couverture pour nos tests
     * 
     * @param[out] V  valeur du portefeuille de couverture
     * @param[out] erreur_couverture  erreur de couverture
     * @param[out] price
     * @param[in] G  vecteur fixe représentant l'aléatoire pour nos tests
     */
    void simu_couverture(PnlVect *V, double &erreur_couverture, PnlVect *price, PnlVect *G);
    
    /**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] market contient une trajectoire du modèle.
     * C'est une matrice de taille (H+1) x d
     * @param[in] nbAsset nombre de titres
     */
    void simul_market(int nbAsset, PnlMat *path);

};


