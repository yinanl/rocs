/**
 *  carReach2.cpp (multiple obstacles)
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
    double xub[] = {10, 10, M_PI};

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
    double OBS[14][3] = {
			{8.2, 0, -M_PI},
			{8.4, 8.5, M_PI}, // block 1

			{8.4, 8.3, -M_PI},
			{9.3, 8.5, M_PI},

			{9.3, 7.1, -M_PI},
			{10.0, 7.3, M_PI},

			{8.4, 5.9, -M_PI},
			{9.3, 6.1, M_PI},

			{9.3, 4.7, -M_PI},
			{10.0, 4.9, M_PI},

			{8.4, 3.5, -M_PI},
			{9.3, 3.7, M_PI},

			{9.3, 2.3, -M_PI},
			{10.0, 2.5, M_PI}
    }; // obstacles


    /* solve the problem */
    rocs::CSolver solver(&carReach);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    for (size_t i = 0; i < 7; ++i) {
	solver.init(rocs::AVOID, OBS[2*i], OBS[2*i+1]);
    }
    solver.init_avoid_area();
    double eta[] = {0.2, 0.2, 0.2};
    solver.reachability_control(&carReach, eta);
    solver.print_controller_info();

    /* save the problem data and the solution */
    // rocs::matWriter wtr("data_carReach2.mat");
    // wtr.open();
    // wtr.write_problem_setting(carReach, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string datafile = "controller_carReach2.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTCntlSys<carde> >(carReach);
    ctlrWtr.write_ivec_array(solver._goal, "G");
    ctlrWtr.write_ivec_array(solver._obs, "xobs");
    ctlrWtr.write_sptree_controller(solver);

    return 0;
}
