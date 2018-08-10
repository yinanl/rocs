/**
 *  carParking.cpp
 *
 *  Automatic parallel parking of a vehicle.
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




int main() {

    /* set the state space */
    double alb = -M_PI/2.5; double aub = M_PI/2.5;
    
    double xlb[] = {0, 0, alb};
    double xub[] = {8, 4, aub};

    /* set the control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};
    
    /* define the control system */
    rocs::DTCntlSys<carde> carParking("parallel parking", h, carde::n, carde::m);
    carParking.init_workspace(xlb, xub);
    carParking.init_inputset(mu, ulb, uub);

    /* load problem settings */
    double ceps = 0.02; // precision for collision area over-approximation
    double glb[] = {L+L/2.0+ceps, 0.5, -M_PI*3.0/180.0};
    double gub[] = {L+L/2.0+d-2*ceps, 0.6, M_PI*3.0/180.0};
    
    rocs::CSolver solver(&carParking);
    solver.init(rocs::GOAL, glb, gub);
    solver.init(rocs::AVOID, &cst_v5<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v6<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v7<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v8<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v1rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v1front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v1curb<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v2rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v2front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v2curb<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v3rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v3front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v3curb<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v4rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v4front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v4curb<rocs::ivec>, ceps);
    solver._ctlr.retract();

    // /* setting 1: wide space */
    // double glb[] = {3.25, 0.5, -M_PI*3.0/180.0};  // goal area
    // double gub[] = {3.75, 0.6, M_PI*3.0/180.0};
    
    // double olb1[] = {0.0, 0.0, alb};  // rear car: [0,2]x[0,1]
    // double oub1[] = {3.0, 1.7, aub};
    // double olb2[] = {5.0, 0.0, alb};  // front car: [6,8]x[0,1]
    // double oub2[] = {8.0, 1.7, aub};

    // /* setting 2: narrow space */
    // double glb[] = {3.02, 0.5, -M_PI*3.0/180.0};  // goal area
    // double gub[] = {3.22, 0.6, M_PI*3.0/180.0};
    
    // double olb1[] = {0.0, 0.0, alb};  // rear car: [0,2]x[0,1]
    // double oub1[] = {3.0, 1.7, aub};
    // double olb2[] = {3.3, 0.0, alb};  // front car: [4.3,6.3]x[0,1]
    // double oub2[] = {6.3, 1.7, aub};

    // rocs::CSolver solver(&carParking);
    // solver.init(rocs::GOAL, glb, gub);
    // solver.init(rocs::AVOID, olb1, oub1);
    // Solver.init(rocs::AVOID, olb2, oub2);

    // solver.cobuchi(&carParking, 0.1, rocs::RELMAXG, 0.2, rocs::RELMAXW, 0.04, true);
    // solver.print_controller_info();


    /* solve the problem */
    solver.cobuchi(&carParking, 0.04, rocs::RELMAXG, 0.2, rocs::RELMAXW, 0.06, true);
    solver.print_controller_info();

    /* save data to file */
    rocs::matWriter wtr("data_carParking.mat");
    wtr.open();
    wtr.write_real_number(H, "H");
    wtr.write_real_number(L, "L");
    wtr.write_real_number(D, "D");
    wtr.write_real_number(d, "d");
    wtr.write_problem_setting(carParking, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();
    
    
    return 0;
}
