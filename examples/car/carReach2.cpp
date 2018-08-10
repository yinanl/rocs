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

#include "src/matlabio.h"

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
    
    double olb1[] = {8.2, 0, -M_PI};  // obstacles
    double oub1[] = {8.4, 8.5, M_PI};
    
    double olb2[] = {8.4, 8.3, -M_PI};
    double oub2[] = {9.3, 8.5, M_PI};
    
    double olb3[] = {9.3, 7.1, -M_PI};
    double oub3[] = {10.0, 7.3, M_PI};
    
    double olb4[] = {8.4, 5.9, -M_PI};
    double oub4[] = {9.3, 6.1, M_PI};

    double olb5[] = {9.3, 4.7, -M_PI};
    double oub5[] = {10.0, 4.9, M_PI};
    
    double olb6[] = {8.4, 3.5, -M_PI};
    double oub6[] = {9.3, 3.7, M_PI};

    double olb7[] = {9.3, 2.3, -M_PI};
    double oub7[] = {10.0, 2.5, M_PI};

    
    /* solve the problem */
    rocs::CSolver solver(&carReach);
    solver.init(rocs::GOAL, glb, gub);
    solver.init(rocs::AVOID, olb1, oub1);
    solver.init(rocs::AVOID, olb2, oub2);
    solver.init(rocs::AVOID, olb3, oub3);
    solver.init(rocs::AVOID, olb4, oub4);
    solver.init(rocs::AVOID, olb5, oub5);
    solver.init(rocs::AVOID, olb6, oub6);
    solver.init(rocs::AVOID, olb7, oub7);
    solver.reachability_control(&carReach, 0.2, rocs::ABSMAX);
    solver.print_controller_info();
    
    /* save the problem data and the solution */
    rocs::matWriter wtr("data_carReach2.mat");
    wtr.open();
    wtr.write_problem_setting(carReach, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();
    
    return 0;
}
