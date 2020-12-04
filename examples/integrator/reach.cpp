/**
 *  reach.cpp
 *
 *  Reachability control of a double integrator.
 *
 *  Created by Yinan Li on July 14, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include "src/system.hpp"
#include "src/csolver.h"
#include "src/matlabio.h"



struct integrator {
  static const int n = 2;  // system dimension
  static const int m = 1;  // control dimension
    
  /* template constructor
   * @param[out] dx
   * @param[in] x = [pos, vel]
   * @param u = [acc]
   */
  template<typename S>
  integrator(S &dx, const S &x, rocs::Rn u) {
    dx[0] = x[0] + x[1] + 0.5*u[0];
    dx[1] = x[1] + u[0];
  }
};


int main()
{
    /**
     * Set the state and control space. 
     */
    double xlb[] = {-10, -5};
    double xub[] = {1.85, 5};
    double ulb[] = {-2.0};
    double uub[] = {2.0};
    double mu[] = {0.3};

    rocs::DTCntlSys<integrator> reach("reach", 1, integrator::n, integrator::m);
    reach.init_workspace(xlb, xub);
    reach.init_inputset(mu, ulb, uub);

    rocs::CSolver solver(&reach, 0, rocs::RELMAX);
    double glb[]{-3.5, -0.5};
    double gub[]{-2.5, 0.5};
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    
    double eta[]{0.2, 0.2};
    solver.reachability_control(&reach, eta);
    solver.print_controller_info();

    rocs::matWriter wtr("data_reach.mat");
    wtr.open();
    wtr.write_problem_setting(reach, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();
    
    return 0;
}
