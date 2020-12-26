/**
 *  carFullmap.cpp
 *
 *  Reachability control of a vehicle, the example taken from
 *  
 *
 *  Created by Yinan Li on March 29, 2019.
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
    const double theta = 3.5;
    double xlb[] = {0, 0, -theta};
    double xub[] = {10, 10, theta};

    /* set the control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* define the control system */
    rocs::DTCntlSys<carde> carReach("reach goal", h, carde::n, carde::m);
    carReach.init_workspace(xlb, xub);
    carReach.init_inputset(mu, ulb, uub);

    /* set precision */
    const double e[] = {0.2, 0.2, 0.2};

    /* set the specifications */
    double glb[] = {9, 0, -theta};  // goal area
    double gub[] = {9.5, 0.5, theta};
    double OBS[30][3] = {
			 {1, 0, -theta},
			 {1.2, 9, theta},

			 { 2.2, 0, -theta},
			 { 2.4, 5, theta},
			 
			 { 2.2, 6, -theta },
			 { 2.4, 10, theta },
			 
			 { 3.4, 0, -theta },
			 { 3.6, 9, theta },
			 
			 { 4.6, 1, -theta },
			 { 4.8, 10, theta },
			 
			 { 5.8, 0, -theta },
			 { 6, 6, theta },
			 
			 { 5.8, 7, -theta },
			 { 6, 10, theta },
			 
			 { 7, 1, -theta },
			 { 7.2, 10, theta },
			 
			 { 8.2, 0, -theta},
			 { 8.4, 8.5, theta},
			 
			 { 8.4, 8.3, -theta},
			 { 9.3, 8.5, theta},
			 
			 { 9.3, 7.1, -theta},
			 { 10,  7.3, theta},
			 
			 { 8.4, 5.9, -theta},
			 { 9.3, 6.1, theta},
			 
			 { 9.3, 4.7, -theta},
			 { 10,  4.9, theta},
			 
			 { 8.4, 3.5, -theta},
			 { 9.3, 3.7, theta},
			 
			 { 9.3, 2.3, -theta},
			 { 10, 2.5, theta}
    }; // obstacles
    
    /* initialize problem setting */
    rocs::CSolver solver(&carReach);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    double OBS1[30][3];
    for (size_t i = 0; i < 15; ++i) {
	OBS1[2*i][0] = OBS[2*i][0] - e[0]/2.0 < xlb[0]? OBS[2*i][0] : OBS[2*i][0] - e[0]/2.0;
	OBS1[2*i][1] = OBS[2*i][1] - e[1]/2.0 < xlb[1]? OBS[2*i][1] : OBS[2*i][1] - e[1]/2.0;
	OBS1[2*i][2] = OBS[2*i][2];
	OBS1[2*i+1][0] = OBS[2*i+1][0] + e[0]/2.0 > xub[0]? OBS[2*i+1][0] : OBS[2*i+1][0] + e[0]/2.0;
	OBS1[2*i+1][1] = OBS[2*i+1][1] + e[1]/2.0 > xub[1]? OBS[2*i+1][1] : OBS[2*i+1][1] + e[1]/2.0;
	OBS1[2*i+1][2] = OBS[2*i+1][2];
    	solver.init(rocs::AVOID, OBS1[2*i], OBS1[2*i+1]);
    }
    solver.init_avoid_area();
    solver.print_controller_info();
    
    /* solve the problem: use interval paver_test */
    solver.reachability_control(&carReach, e);
    solver.print_controller_info();
    
    /* save the problem data and the solution */
    // rocs::matWriter wtr("data_carFullmap.mat");
    // wtr.open();
    // wtr.write_problem_setting(carReach, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string datafile = "controller_carFullmap.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTCntlSys<carde> >(carReach);
    ctlrWtr.write_ivec_array(solver._goal, "G");
    ctlrWtr.write_ivec_array(solver._obs, "xobs");
    ctlrWtr.write_sptree_controller(solver);


    return 0;
}
