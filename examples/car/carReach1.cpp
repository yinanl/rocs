/**
 *  carReach1.cpp (a single obstacle)
 *
 *  Reachability control of a vehicle.
 *
 *  Created by Yinan Li on May 12, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>

#include "src/system.hpp"
#include "src/csolver.h"
#include "src/hdf5io.h"
// #include "src/matlabio.h"

#include "car.hpp"


int main()
{
    /* set the state space */
    double xlb[] = {7.3, 0, -M_PI};
    double xub[] = {10, 2, M_PI};

    /* set the control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* define the control system */
    rocs::DTCntlSys<carde> carReach("reach goal", h, carde::n, carde::m);
    carReach.init_workspace(xlb, xub);
    carReach.init_inputset(mu, ulb, uub);

    /* set the specifications */
    double glb[] = {9, 0, -M_PI};  // goal area
    double gub[] = {9.5, 0.5, M_PI};

    double olb[] = {8, 0.3, -M_PI};  // obstacle
    double oub[] = {8.4, 1.2, M_PI};

    /* solve the problem */
    rocs::CSolver solver(&carReach);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    solver.init(rocs::AVOID, olb, oub);
    solver.init_avoid_area();
    double eta[] = {0.2, 0.2, 0.2};
    // double emin[] = {0.1, 0.1, 0.1};
    solver.reachability_control(&carReach, eta);
    solver.print_controller_info();

    /* save the problem data and the solution */
    // rocs::matWriter wtr("data_carReach1.mat");
    // wtr.open();
    // wtr.write_problem_setting(carReach, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string datafile = "controller_carReach1.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTCntlSys<carde> >(carReach);
    ctlrWtr.write_ivec_array(solver._goal, "G");
    ctlrWtr.write_ivec_array(solver._obs, "xobs");
    ctlrWtr.write_sptree_controller(solver);

    return 0;
}
