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
    rocs::CSolver solver(&dcdcRS, rocs::RELMAX);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();

    double ei[]{0.001, 0.0002};
    double er[]{0.04, 0.008};
    double emin[]{0.005, 0.001};
    // solver.cobuchi(&dcdcRS, ei, er, true, emin);
    solver.cobuchi(&dcdcRS, ei, ei);
    // solver.reachstay_control(&dcdcRS, ei, er, true, emin);
    solver.print_controller_info();


    /* save the problem data and the solution */
    rocs::matWriter wtr("data_dcdcCoBuchi.mat");
    // rocs::matWriter wtr("data_ReachStay.mat");
    wtr.open();
    wtr.write_problem_setting(dcdcRS, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    return 0;
}
