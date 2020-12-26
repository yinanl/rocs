/**
 *  dba_integrator.cpp
 *
 *  DBA control of the simplified SCARA manipulator dynamics (the double integrator model).
 *
 *  Created by Yinan Li on July 14, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include "src/system.hpp"
#include "src/csolver.h"

// #include "src/matlabio.h"
#include "src/hdf5io.h"

#include "scara.hpp"


const double h = 0.8*l1;
const double r = 0.5*l1;
double a1 = atan(h / r); // the upper bound for theta1
double a2 = asin(h / l1);
/**
 * Define constraints for theta2 (for current r, h) in the form of
 * f(x)\leq 0
 *
 */
template<typename T>
T collision1(const T &x) {
    T y{x[0]-a2, -x[0]-x[1]+M_PI-atan( (h-l1*sin(x[0]))/(l1*cos(x[0])-r) )};
    return y;
}

template<typename T>
T collision2(const T &x) {
    T y{x[0]-a1, a2-x[0],
	-x[0]-x[1]+M_PI+atan( (l1*sin(x[0])-h)/(l1*cos(x[0])) )};
    return y;
}


int main()
{
    /**
     * One double integrator.
     */
    // double xlb[] = {-M_PI, -1};
    // double xub[] = {M_PI, 1};
    // double ulb[] = {-2.0};
    // double uub[] = {2.0};
    // double mu[] = {0.1};
    // rocs::DTCntlSys<simple> reach("reach", ts, simple::n, simple::m);
    // reach.init_workspace(xlb, xub);
    // reach.init_inputset(mu, ulb, uub);

    // // double goal[2][2] = {{0.4080, -0.1}, {0.6172, 0.1}};
    // double goal[2][2] = {{1.5739, -0.1}, {1.7355, 0.1}};
    // rocs::CSolver solver(&reach, rocs::RELMAX, 100);
    // solver.init(rocs::GOAL, goal[0], goal[1]);
    // solver.init_goal_area();
    // double eta[]{0.05, 0.05};
    // solver.reachability_control(&reach, eta);
    // solver.print_controller_info();

    // rocs::matWriter wtr("data_simple_reach.mat");
    // wtr.open();
    // wtr.write_problem_setting(reach, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();

    /**
     * Set the state and control space.
     */
    double xlb[] = {0, -M_PI, -1, -1};
    double xub[] = {M_PI/2.0, M_PI, 1, 1};

    double ulb[] = {-5.0, -5.0};
    double uub[] = {5.0, 5.0};
    double mu[] = {0.5, 0.5};

    rocs::DTCntlSys<integrator> reach("reach", ts, integrator::n, integrator::m);
    reach.init_workspace(xlb, xub);
    reach.init_inputset(mu, ulb, uub);

    rocs::CSolver solver(&reach, 0, rocs::RELMAX, 80);
    double goal[2][4] = {{0.4980, 1.5739, -0.1, -0.1},
    			 {0.5772, 1.7055, 0.1, 0.1}}; // (0.05,0.2)=(0.5126,1.6264)
    // double goal[2][4] = {{0.4080, 1.5739, -0.1, -0.1},
    // 			 {0.6172, 1.7355, 0.1, 0.1}};
    // double goal[2][4] = {{0.4903, -0.9363, -0.1, -0.1},
    // 			 {0.6069, -0.8363, 0.1, 0.1}}; //(0.27,0.03)=(0.5488,-0.8763)
    solver.init(rocs::GOAL, goal[0], goal[1]);
    solver.init_goal_area();
    double obs[2][4] = {{a1, xlb[1], xlb[2], xlb[3]},
    			{M_PI/2.0, xub[1], xub[2], xub[3]}};
    solver.init(rocs::AVOID, obs[0], obs[1]);
    const double e[]{0.01, 0.01, 3, 3};  // only bisect x[0] and x[1].
    solver.init(rocs::AVOID, &collision1<rocs::ivec>, e);
    solver.init(rocs::AVOID, &collision2<rocs::ivec>, e);

    double eta[]{0.05, 0.05, 0.1, 0.1};
    solver.reachability_control(&reach, eta);
    solver.print_controller_info();

    // rocs::matWriter wtr("data_2integrators_reach1.mat");
    // wtr.open();
    // wtr.write_problem_setting(reach, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string datafile = "controller_2dbint_reach1_80.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTCntlSys<integrator> >(reach);
    ctlrWtr.write_ivec_array(solver._goal, "G");
    ctlrWtr.write_ivec_array(solver._obs, "xobs");
    ctlrWtr.write_sptree_controller(solver);

    return 0;
}
