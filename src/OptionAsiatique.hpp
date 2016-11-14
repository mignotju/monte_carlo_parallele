/* 
 * File:   OptionAsiatique.hpp
 * Author: gononl
 *
 * Created on September 17, 2016, 4:27 PM
 */

#ifndef OPTIONASIATIQUE_HPP
#define	OPTIONASIATIQUE_HPP

#include "Option.hpp"

#include "pnl/pnl_vector.h"

/*!
 * \file OptionAsiatique.hpp
 * \brief Option Asiatique discrète
 */
class OptionAsiatique : public Option {
public:
    double K; // Prix d'exercice (Strike)
    PnlVect *sommeStd;

    // Constructeurs
    /**
     * \brief Constructeur par defaut
     */
    OptionAsiatique();
    
    /**
     * \brief Constructeur complet
     * 
     * @param maturity  maturité
     * @param nbTimeSteps  nombre de pas entre 0 et la maturité
     * @param size  nombre de sous-jacents
     * @param strike  prix d'exercice de l'option
     * @param lambda  vecteur représentant la composition de l'option
     */
    OptionAsiatique(double maturity, int nbTimeSteps, int size, double strike, PnlVect* lambda);
    
    /**
    * Calcule la valeur du payoff sur la trajectoire
    *
    * @param[in] path est une matrice de taille (N+1) x d
    * contenant une trajectoire du modèle telle que créée
    * par la fonction asset.
    * @return phi(trajectoire)
    */
    virtual double payoff(const PnlMat *path);
};

#endif	/* OPTIONASIATIQUE_HPP */

