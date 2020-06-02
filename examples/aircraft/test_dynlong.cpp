/**
 *  test_dyn_longitudinal.cpp
 *
 *  The longitudinal kinematics of DC9-30.
 *
 *  Created by Yinan Li on Mar.27, 2020.
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
    
    /* test the reachable set */
    // rocs::ivec x= {rocs::interval(60,60.05),
    // 		   rocs::interval(-M_PI/180,-21/22*M_PI/180),
    // 		   rocs::interval(16,16.2)};
    rocs::ivec x= {rocs::interval(80.0,81.0),
		   rocs::interval(-0.0196,-0.0131),
		   rocs::interval(54.25,56.00)};
    // rocs::Rn u = {0,0};
    std::vector<rocs::ivec> y(aircraft._ugrid._nv, rocs::ivec(3));
    std::cout << "The initial interval: " << x << '\n';
    aircraft.get_reach_set(y, x);
    std::cout << "The reachable set of x: \n";
    // std::cout << y[12] << '\n';
    for (int i = 0; i < y.size(); ++i)
    	std::cout << y[i] <<'\n';

    aircraft.release_flows();
    return 0;
}
