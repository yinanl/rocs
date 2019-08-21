/**
 *  dcdc_reachstay.cpp
 *
 *  A DCDC converter reach-and-stay control example.
 *
 *  Created by Yinan Li on April, 2019.
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
    double xlb[] = {0.649, 0.9898};
    double xub[] = {1.65, 1.19};

    /* define the control system */
    rocs::DTSwSys<dcde> dcdcRS("dcdc", tau, dcde::n, dcde::m);
    dcdcRS.init_workspace(xlb, xub);
    
    /* set the specifications */
    double glb[] = {1.1, 1.08};
    double gub[] = {1.6, 1.18};

    /* solve the problem */
    rocs::CSolver solver(&dcdcRS);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();

    // solver.cobuchi(&dcdcRS, 0.001, rocs::RELMAXG, 0.04, rocs::RELMAXW, 0.005, true);
    solver.cobuchi(&dcdcRS, 0.001, rocs::RELMAXG, 0.001, rocs::RELMAXW);
    // solver.reachstay_control(&dcdcRS, 0.001, rocs::RELMAXG, 0.04, rocs::RELMAXW, 0.005, true);
    // solver.reachability_control(&dcdcRS, 0.003, rocs::ABSMAX);
    solver.print_controller_info();


    /* save the problem data and the solution */
    rocs::matWriter wtr("data_dcdcCoBuchi.mat");
    // rocs::matWriter wtr("data_dcdcReachStay.mat");
    // rocs::matWriter wtr("data_dcdcReach.mat");
    wtr.open();
    wtr.write_problem_setting(dcdcRS, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    return 0;
}
