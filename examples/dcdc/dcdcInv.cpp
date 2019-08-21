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

#include "dcdc.hpp"


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
    solver.init_goal_area();
    
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
