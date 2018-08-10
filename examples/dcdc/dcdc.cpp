/**
 *  dcdc.cpp
 *
 *  A DCDC converter invariance control example.
 *
 *  Created by Yinan Li on May 10, 2018.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"


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


int main()
{
    /* set the state space */
    // double xlb[] = {-2, -1.5};
    // double xub[] = {2, 3};
    double xlb[] = {-2, 0.70};
    double xub[] = {2, 1.50};

    /* define the control system */
    rocs::DTSwSys<dcde> dcdcInv("dcdc", tau, dcde::n, dcde::m);
    dcdcInv.init_workspace(xlb, xub);
    
    /* set the specifications */
    double glb[] = {1.15, 1.09};
    double gub[] = {1.55, 1.17};

    /* solve the problem */
    rocs::CSolver solver(&dcdcInv);
    solver.init(rocs::GOAL, glb, gub);
    solver.invariance_control(&dcdcInv, 0.001, rocs::RELMAXG);
    solver.print_controller_info();


    /* save the problem data and the solution */
    rocs::matWriter wtr("data_dcdcInv.mat");
    wtr.open();
    wtr.write_problem_setting(dcdcInv, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    return 0;
}
