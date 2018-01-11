/**
 *  invertpdl.hpp
 *
 *  Inverted pendulum [Michigan control tutorial](http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling)
 *
 *  Created by yinan li on April 26, 2017.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include "../../src_itvl/csolver.h"

#include "invertpdl.hpp"




int main()
{
    /* state and input space */
    const int XD = 2;
    const int UD = 1;
    double xlb[XD] = {-2, -3.2};
    double xub[XD] = {2, 3.2};
    double ulb[UD] = {-10};
    double uub[UD] = {10};
    double mu[UD] = {0.05};
    double tau = 0.01;  // sampling time
    double dt = 0.002;  // ode integration step

    /* inverted pendulum dynamics */
    ipdl *ptrIpdl = new ipdl(UD, ulb, uub, mu, tau, dt);

    /* specification */
    double glb[] = {-0.05, -0.01};
    double gub[] = {0.05, 0.01};

    /* define control problem */
    CntlProb ipdlReachStay("invertpdl", XD, UD, xlb, xub, ptrIpdl);

    /* create a solver */
    CSolver *solver = new CSolver(&ipdlReachStay, 100);
    solver->init(GOAL, glb, gub);
    
    /* solve */
    // solver->reach_stay(0.001, RELMAXG, 0.04, RELMAXW, 0.002, true);
    solver->cobuchi(0.001, RELMAXG, 0.04, RELMAXW, 0.002, true);

    solver->print_controller_info();

    /* save the specification and controller */
    ipdlReachStay.write2mat_settings("data_ipdl_spec.mat");
    // solver->serialize_controller("data_ipdl_ctree_cobuchi.mat");
    solver->write2mat_controller("data_ipdl_cbox_cobuchi.mat");
    

    delete solver;
    delete ptrIpdl;

    return 1;
}
