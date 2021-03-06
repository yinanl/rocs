/**
 *  dcdc.hpp
 *
 *  Dynamics of a boost dcdc converter.
 *
 *  Created by Yinan Li on Aug. 21, 2019.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _dcdc_h
#define _dcdc_h

#include "src/interval_vector.h"


/* Parameters of the model */
const double tau = 0.5;
const double xc = 70.0;
const double xl = 3.0;
const double rc = 0.005;
const double rl = 0.05;
const double r0 = 1.0;
const double vs = 1.0;

arma::mat I = arma::eye<arma::mat>(2, 2);
arma::vec b = {vs/xl, 0};
arma::mat A1 = {{-rl/xl, 0}, {0, -1/(xc*(rc+r0))}};
arma::mat F1 = arma::expmat(A1 * tau);
arma::vec g1 = arma::inv(A1) * (F1 - I) * b;
arma::mat A2 = { {(-1/xl)*(rl+r0*rc/(r0+rc)),(-1/xl)*(r0/(r0+rc))}, {(1/xc)*(r0/(r0+rc)),(-1/xc)*(1/(r0+rc))} };
arma::mat F2 = arma::expmat(A2 * tau);
arma::vec g2 = arma::inv(A2) * (F2 - I) * b;


/* Discrete-time dynamics of DCDC converter */
struct dcde {
    static const int n = 2;  // state dimension
    static const int m = 2;  // number of modes
    
    /**
     * Constructors: 
     * real-valued (arma::vec) and interval-valued (rocs::ivec)
     * @param[out] y the next state after the sampling time.
     * @param[in] x the current state.
     * @param[in] m the mode.
     */
    dcde(arma::vec &y, const arma::vec &x, const int m) {
	switch (m) {
	case 1:
	    y = F1*x+g1;
	    break;
	case 2:
	    y = F2*x+g2;
	    break;
	default:
	    break;
	}
    }

    dcde(rocs::ivec &y, const rocs::ivec &x, const int m) {
	switch (m) {
	case 1:
	    y = linmap(F1, g1, x);
	    break;
	case 2:
	    y = linmap(F2, g2, x);
	    break;
	default:
	    break;
	}
    }
    
};

#endif
