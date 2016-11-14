/* 
 * File:   OptionBasket.hpp
 * Author: gononl
 *
 * Created on September 16, 2016, 10:03 AM
 */

#ifndef OPTIONBASKET_HPP
#define	OPTIONBASKET_HPP

#include "Option.hpp"

#include "pnl/pnl_vector.h"

/*!
 * \file OptionBasket.hpp
 * \brief Option Basket
 */
class OptionBasket : public Option {
public:
    double K; // Prix d'exercice (Strike)
    PnlVect *ST;
    
    /**
     * \brief Constructeur par defaut
     */
    OptionBasket();
    
    /**
     * \brief Constructeur complet
     *  
     * @param maturity  maturité
     * @param nbTimeSteps  nombre de pas entre 0 et la maturité
     * @param size  nombre de sous-jacents 
     * @param strike  prix d'exercice de l'option
     * @param lambda  composition de l'option
     */
    OptionBasket(double maturity, int nbTimeSteps, int size, double strike, PnlVect* lambda);
    
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

#endif	/* OPTIONBASKET_HPP */

