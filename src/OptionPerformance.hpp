/* 
 * File:   OptionPerformance.hpp
 * Author: gononl
 *
 * Created on September 17, 2016, 5:06 PM
 */

#ifndef OPTIONPERFORMANCE_HPP
#define	OPTIONPERFORMANCE_HPP

#include "Option.hpp"

/*!
 * \file OptionPerformance.hpp
 * \brief Option Performance sur panier
 */
class OptionPerformance : public Option {
public:    
    PnlVect *cours;
    
    /**
     * \brief Constructeur par defaut
     */
    OptionPerformance();
      
    /**
     * \brief Constructeur complet
     * 
     * @param maturity  maturité
     * @param nbTimeSteps  nombre de pas entre 0 et T 
     * @param size  nombre de sous-jacents de l'option
     * @param lambda  vecteur représentant la composition de l'option
     */
    OptionPerformance(double maturity, int nbTimeSteps, int size, PnlVect* lambda);
    
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

#endif	/* OPTIONPERFORMANCE_HPP */

