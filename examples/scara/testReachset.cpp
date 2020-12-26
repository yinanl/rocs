/**
 *  testReachset.cpp
 *
 *  Test the reachable set computation of a two-link scara manipulator.
 *
 *  Created by Yinan Li on July 05, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include <array>
#include <boost/numeric/odeint.hpp>

#include "src/system.hpp"
#include "src/csolver.h"

#include "scara.h"

typedef std::array<double, 4> state_type;
struct scara_vf {
    rocs::Rn _u;
    scara_vf(rocs::Rn u) : _u(u) {}
    void operator()(state_type &x, state_type &dxdt, double t) const {
	double z1 = I1+I2+m1*r1*r1+m2*(l1*l1+r2*r2);
	double z2 = m2*l1*r2;
	double z3 = I2+m2*r2*r2;
	double detM = z3*(z1-z3) - z2*z2*std::cos(x[1])*std::cos(x[1]);
	double a = z2*std::sin(x[1])*(2*x[2]+x[3])*x[3];
	double b = z2*std::cos(x[1]);
	double c = z2*x[2]*std::sin(x[1])-_u[1];
	
	dxdt[0] = x[2];
	dxdt[1] = x[3];
	dxdt[2] = (z3*_u[0] + z3*a + (z3+b)*c) / detM;
	dxdt[3] = ((z1+2*b)*(-c) - (z3+b)*(_u[0]+a)) / detM;
    }
};


int main() {

    /* set the state space */
    double xlb[4] = {0, -M_PI, -0.5, -0.5};
    double xub[4] = {M_PI, M_PI, 0.5, 0.5};
    
    /* set the control values */
    double ulb[2] = {-0.001, -0.001};
    double uub[2] = {0.001, 0.001};
    double mu[2] = {0.0002, 0.0002};

    /* set the sampling time and disturbance */
    double t = 0.05;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    
    /* define the control system */
    rocs::CTCntlSys<scaraode> scara("inverted pendulum", t,
				    scaraode::n, scaraode::nu,
				    delta, &controlparams);
    scara.init_workspace(xlb, xub);
    scara.init_inputset(mu, ulb, uub);
    scara.allocate_flows();

    /* test if reachable set covers the nominal trajectory */
    /* compute reachable set */
    // rocs::ivec x0 = {rocs::interval(0.09, 0.11),
    // 		     rocs::interval(-0.01, 0.01),
    // 		     rocs::interval(-0.01, 0.01),
    // 		     rocs::interval(-0.01, 0.01)};
    rocs::ivec x0 = {rocs::interval(0.6172, 0.63),
    		     rocs::interval(1.61, 1.63),
    		     rocs::interval(-0.005, 0.005),
		     rocs::interval(-0.005, 0.005)};
    std::vector<rocs::ivec> x(scara._ugrid._nv, rocs::ivec(4));
    std::cout << "The initial interval: " << x0 << '\n';
    std::cout << "The integrating time: " << t << '\n';
    scara.get_reach_set(x, x0);

    // double dt{0.001};
    // for (size_t i = 0; i < x.size(); ++i) {
    // 	/* integrate the nominal trajectory */
    // 	state_type y{0.1, 0, 0, 0};
    // 	rocs::Rn u(scara._ugrid._data[i]);
    // 	boost::numeric::odeint::runge_kutta_cash_karp54<state_type> rk45;
    // 	boost::numeric::odeint::integrate_const(rk45, scara_vf(u),
    // 						y, 0.0, t, dt);
    // 	std::cout << "x(t)= [" << y[0] << ','<< y[1] << ',' << y[2] << ',' << y[3] << "]\n";
    // 	std::cout << "R(t, x0, [" << scara._ugrid._data[i][0] << ','
    // 		  << scara._ugrid._data[i][1] << "])= ";
    // 	std::cout << x[i] <<'\n';
    // 	std::cout << "Check: ";
    // 	for (int j = 0; j < 4; ++j) {
    // 	    if (y[j] > x[i][j].getinf() && y[j] < x[i][j].getsup())
    // 		std::cout << true << ' ';
    // 	    else
    // 		std::cout << false << ' ';
    // 	}
    // 	std::cout << '\n';
    // }
    
    scara.release_flows();
    
    return 0;
}
