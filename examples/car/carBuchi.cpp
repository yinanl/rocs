/**
 *  carBuchi.cpp
 *
 *  Buchi objective of a vehicle.
 *  
 *  Created by Yinan Li on May 12, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"

#include "car.hpp"


int main()
{
    /* State space */
    double xlb[] = {7.3, 0, -M_PI};
    double xub[] = {10, 2, M_PI};

    /* Control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* Control system */
    rocs::DTCntlSys<carde> carBuchi("Buchi", h, carde::n, carde::m);
    carBuchi.init_workspace(xlb, xub);
    carBuchi.init_inputset(mu, ulb, uub);

    /* Specifications */
    double glb[] = {9, 0, -M_PI};  // goal area
    double gub[] = {9.5, 0.5, M_PI};
    
    double olb[] = {8, 0.3, -M_PI};  // obstacle
    double oub[] = {8.4, 1.2, M_PI};
    
    /* Solve the problem */
    rocs::CSolver solver(&carBuchi);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    solver.init(rocs::AVOID, olb, oub);
    solver.init_avoid_area();
    const double eta[]{0.2, 0.2, 0.2};
    solver.buchi(&carBuchi, eta);
    solver.print_controller_info();
    
    /* Save the specification and controller */
    rocs::matWriter wtr("data_carBuchi.mat");
    wtr.open();
    wtr.write_problem_setting(carBuchi, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();
    
    return 0;
}
