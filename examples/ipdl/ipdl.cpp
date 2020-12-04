/**
 *  ipdl.cpp
 *
 *  The example of stabilization of an inverted pendulum.
 *
 *  Created by Yinan Li on May 10, 2018.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"


const double m = 0.2;
const double b = 0.1;
const double J = 0.006;
const double g = 9.8;
const double l = 0.3;

/* user defined dynamics */
struct ipdlode {

    static const int n = 2;  // system dimension
    static const int nu = 1;  // control dimension
    
    double a1; // = (m*g*l) / (J + m*l*l);
    double a2; // = b / (J + m*l*l);
    double a3; // = l / (J + m*l*l);

    /* template constructor
     * @param[out] dx
     * @param[in] x
     * @param u
     */
    template<typename S>
    ipdlode(S *dx, const S *x, rocs::Rn u) :
	a1((m*g*l) / (J + m*l*l)), a2(b / (J + m*l*l)), a3(l / (J + m*l*l)) {
	dx[0] = x[1];
	dx[1] = a1*sin(x[0]) - a2*x[1] + a3*cos(x[0])*u[0];
    }
};


int main() {

    /* set the state space */
    double xlb[2] = {-2, -3.2};
    double xub[2] = {2, 3.2};
    
    /* set the control values */
    double ulb[1] = {-10};
    double uub[1] = {10};
    double mu[1] = {0.05};

    /* set the sampling time and disturbance */
    double t = 0.01;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    
    /* define the control system */
    rocs::CTCntlSys<ipdlode> ipdl("inverted pendulum", t, ipdlode::n, ipdlode::nu, delta, &controlparams);
    ipdl.init_workspace(xlb, xub);
    ipdl.init_inputset(mu, ulb, uub);
    ipdl.allocate_flows();

    /* set the specifications */
    double glb[] = {-0.05, -0.01};
    double gub[] = {0.05, 0.01};

    /* solve the problem */
    rocs::CSolver solver(&ipdl, 0, rocs::RELMAX, 100);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    double ei[]{0.001, 0.001};
    double er[]{0.04, 0.04};
    double ermin[]{0.002, 0.002};
    solver.cobuchi(&ipdl, ei, er, true, ermin);
    solver.print_controller_info();

    /* save the problem data and the solution */
    rocs::matWriter wtr("data_ipdlCobuchi.mat");
    wtr.open();
    wtr.write_problem_setting(ipdl, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();
    
    ipdl.release_flows();
    
    return 0;
}
