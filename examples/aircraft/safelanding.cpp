/**
 *  safelanding.cpp
 *
 *  Safe landing of DC9-30.
 *
 *  Created by Yinan Li on Feb. 10, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>

#include "src/system.hpp"
#include "src/csolver.h"
#include "src/matlabio.h"


/* Parameters of the model */
double mg = 60000.0*9.81;
double mi = 1.0/60000; /* weight inverse: 1/m */



/* ODE of the longitudinal equation of motions for DC9-30 */
struct eomlong {
    static const int n = 3;  // state dimension
    static const int m = 2;  // control dimension
    
    /**
     * Constructors: 
     * @param[out] dx (dV, dgamma, dh)
     * @param[in] x (V,gamma,h).
     * @param[in] u (T, alpha) control input (thrust, angle of attack)
     */
    template<typename S>
    eomlong(S *dx, const S *x, rocs::Rn u) {
	double c = 1.25+4.2*u[1];
	dx[0] = mi*(u[0]*cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*sin(x[1]));
	dx[1] = (1.0/(60000*x[0]))*(u[0]*sin(u[1])+68.6*c*x[0]*x[0]-mg*cos(x[1]));
	dx[2] = x[0]*sin(x[1]);
    }
    
};

/** A target set in the form of f(x)<=0:
 * x(0)*sin(x(1)) >= -0.91
 *
 **/
template<typename T>
T target_area(const T &x) {
    T y(7);
    y[0] = -0.91-x[0]*sin(x[1]);
    y[1] = 63 - x[0];
    y[2] = x[0] - 75;
    y[3] = -3*M_PI/180 - x[1];
    y[4] = x[1];
    y[5] = -x[2];
    y[6] = x[2] - 2.5;
    return y;
}

// auto target = [](const double[] x) {
// 		  if(63 <= x[0] &&  x[0] <=  75 &&
// 		     -3*M_PI/180 <= x[1] && x[1] <=   0 &&
// 		     0 <= x[2] && x[2] <= 2.5) {
// 		      return true;
// 		  }
// 		  return false;
// 	      };


int main(int argc, char *argv[])
{   
    /**
     * Define the control system 
     **/
    /* Set sampling time and disturbance */
    double tau = 0.25;
    double delta = 10;
    /* Set parameters for computation */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    rocs::CTCntlSys<eomlong> aircraft("landing", tau, eomlong::n,
				      eomlong::m, delta, &controlparams);
    
    /* set the state space */
    double xlb[] = {58, -3*M_PI/180, 0};
    double xub[] = {83, 0, 56};
    double ulb[] = {0,0};
    double uub[] = {32000, 8*M_PI/180};
    double mu[] = {32000, 9.0/8.0*M_PI/180};
    aircraft.init_workspace(xlb, xub);
    aircraft.init_inputset(mu, ulb, uub);
    // std::cout << "# of u= " << aircraft._ugrid._nv
    // 	      << ", dimension= " << aircraft._ugrid._dim << '\n';
    // for(int i = 0; i < aircraft._ugrid._nv; ++i) {
    // 	std::cout << aircraft._ugrid._data[i][0] << ','
    // 		  << aircraft._ugrid._data[i][1] << '\n';
    // }
    aircraft.allocate_flows();

    /* Set the target set */
    const double eta[] = {25.0/362*2, 3*M_PI/180/66*2, 56.0/334*2};
    const double emin[] = {2.5/362,0.3*M_PI/180/66,5.6/334};
    double glb[] = {63, -3*M_PI/180, 0};
    double gub[] = {75, 0, 2.5};

    /* Solve the reachability problem */
    rocs::CSolver solver(&aircraft, rocs::RELMAX, 100);
    // solver.init(rocs::GOAL, glb, gub);
    solver.init(rocs::GOAL, &target_area<rocs::ivec>, eta);
    solver.init_goal_area(); /* save the goal area information */
    solver.reachability_control(&aircraft, eta);
    solver.print_controller_info();

    /* save the problem data and the solution */
    rocs::matWriter wtr("data_safe_landing.mat");
    wtr.open();
    wtr.write_problem_setting(aircraft, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    // /* Simulate the closed-loop control */
    // double x[3] = {81, -M_PI/180, 55};
    // double u[2] = {0, 0};
    // while (1) {
    // 	if (target(x)) {
    // 	    std::cout << "Arrived: " << x[0] <<  " "  << x[1] << " " << x[2] << std::endl;
    // 	    break;
    // 	} else {
    // 	    u[0] = solver._ctlr._cntl[0]
    // 	}
    // }

    // /**
    //  * DBA control synthesis
    //  **/
    // std::vector<rocs::CSolver*> w(nNodes);
    // for (rocs::UintSmall i = 0; i < nNodes; ++ i) {
    // 	w[i] = new rocs::CSolver(&aircraft,nProps);
    // 	w[i]->set_M(arrayM[i]);
    // 	for (rocs::UintSmall j = 1; j < nProps; ++j) {
    // 	    if (arrayM[i][j] < nNodes) {
    // 		w[i]->labeling(P[2*(j-1)], P[2*j+1], j);
    // 	    }
    // 	}
    // 	// std::cout << "Initial partition of w" << i <<":\n";
    // 	// w[i]->print_controller();
    // }
    // rocs::dba_control< rocs::CTCntlSys<eomlong> >(w, &aircraft, nNodes, acc, e);

    // /**
    //  * Display and save memoryless controllers.
    //  */
    // rocs::write_results_to_mat(aircraft, specfile, w);
    
    // /**
    //  * Release dynamic memory 
    //  **/
    // for (rocs::UintSmall i = 0; i < nNodes; ++i) {
    // 	delete w[i];
    // }

    aircraft.release_flows();
    return 0;
}
