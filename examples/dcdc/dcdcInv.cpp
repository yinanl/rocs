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
#include "src/hdf5io.h"
// #include "src/matlabio.h"

#include "dcdc.hpp"


int main()
{
    /* set the state space */
    double xlb[] = {-2, 0.70};
    double xub[] = {2, 1.50};

    /* define the control system */
    rocs::DTSwSys<dcde> dcdcInv("dcdc", tau, dcde::n, dcde::m);
    dcdcInv.init_workspace(xlb, xub);
    
    /* set the specifications */
    double glb[] = {1.15, 1.09};
    double gub[] = {1.55, 1.17};

    /* solve the problem */
    rocs::CSolver solver(&dcdcInv, 0, rocs::RELMAX, UINT_MAX);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    
    double e[]{0.001,0.0002};
    solver.invariance_control(&dcdcInv, e);
    solver.print_controller_info();


    /* save the problem data and the solution */
    // rocs::matWriter wtr("data_dcdcInv.mat");
    // wtr.open();
    // wtr.write_problem_setting(dcdcInv, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string datafile = "controller_dcdcInv.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTSwSys<dcde> >(dcdcInv);
    ctlrWtr.write_ivec_array(solver._goal, "G");
    ctlrWtr.write_ivec_array(solver._obs, "xobs");
    ctlrWtr.write_sptree_controller(solver);

    return 0;
}
