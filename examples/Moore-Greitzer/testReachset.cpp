/**
 *  engine.cpp
 *
 *  Reach and stay control of the ode model of Moore-Greitzer engine.
 *
 *  Created by Yinan Li on June 9, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include <array>
#include <boost/numeric/odeint.hpp>

#include "src/system.hpp"
#include "src/csolver.h"
#include "src/matlabio.h"

typedef std::array<double, 3> state_type;

double a = 1./3.5;
double B = 2.0;
double H = 0.18;
double W = 0.25;
double lc = 8.0;
double cx = 1.0/lc;
double cy = 1.0/(4*lc*B*B);
double aH = a+H;
double H2 = H/(2.0*W*W*W);
double W2 = 3*W*W;

/* user defined dynamics */
struct mgode {
    static const int n = 3;  // system dimension
    static const int nu = 2;  // control dimension

    /* template constructor
     * @param[out] dx
     * @param[in] x
     * @param u
     */
    template<typename S>
    mgode(S *dx, const S *x, rocs::Rn u) {
	dx[0] = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]) + u[0];
	dx[1] = cy * (x[0] - x[2]*sqrt(x[1]));
	dx[2] = u[1];
    }
};

struct engine_vf {
    rocs::Rn _u;
    engine_vf(rocs::Rn u) : _u(u) {}
    void operator()(state_type &x, state_type &dxdt, double t) const {
	dxdt[0] = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]) + _u[0];
	dxdt[1] = cy * (x[0] - x[2]*sqrt(x[1]));
	dxdt[2] = _u[1];
    }
};

int main() {

    /* set the state space */
    double xlb[3]{0, 0, 0.5};
    double xub[3]{1, 1, 0.8};
    
    /* set the control values */
    // double L = 0.0003;
    // double ulb[2]{-5, -L};
    // double uub[2]{5, L};
    // double mu[2]{1, L/10};
    double Lmu = 0.01;
    double Lu = 0.05;
    double ulb[2]{-Lu, -Lmu};
    double uub[2]{Lu, Lmu};
    double mu[2]{Lu/5, 2*Lmu/5};

    /* set the sampling time and disturbance */
    double t = 0.01;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 10;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    
    /* define the control system */
    rocs::CTCntlSys<mgode> engine("Moore-Greitzer", t, mgode::n, mgode::nu,
				  delta, &controlparams);
    engine.init_workspace(xlb, xub);
    engine.init_inputset(mu, ulb, uub);
    engine.allocate_flows();

    
    /* test if reachable set covers the nominal trajectory */
    /* compute reachable set */
    double dt{0.001};
    // state_type y{0.5039, 0.6605, 0.62};
    state_type y{1.526/2, 0.5, 0.65};
    // double e[3]{0.2, 0.5, 0.1};
    // rocs::ivec x0 = {rocs::interval(y[0]-e[0], y[0]+e[0]),
    // 		     rocs::interval(y[1]-e[1], y[1]+e[1]),
    // 		     rocs::interval(y[2]-e[2], y[2]+e[2])};
    rocs::ivec x0 = {rocs::interval(0.526, 1),
    		     rocs::interval(0.001, 1),
    		     rocs::interval(0.5, 0.8)};
    std::vector<rocs::ivec> x(engine._ugrid._nv, rocs::ivec(3));
    std::cout << "The initial interval: " << x0 << '\n';
    std::cout << "The integrating time: " << t << '\n';
    engine.get_reach_set(x, x0);
    
    for (size_t i = 0; i < x.size(); ++i) {
	/* integrate the nominal trajectory */
	state_type y{1.526/2, 0.5, 0.65};
	rocs::Rn u(engine._ugrid._data[i]);
	boost::numeric::odeint::runge_kutta_cash_karp54<state_type> rk45;
	boost::numeric::odeint::integrate_const(rk45, engine_vf(u),
						y, 0.0, t, dt);
	std::cout << "x(t)= [" << y[0] << ','<< y[1] << ',' << y[2] << "]\n";
	std::cout << "R(t, x0, [" << engine._ugrid._data[i][0] << ','
		  << engine._ugrid._data[i][1] << "])= ";
    	std::cout << x[i] <<'\n';
	std::cout << "Check: ";
	for (int j = 0; j < 3; ++j) {
	    if (y[j] > x[i][j].getinf() && y[j] < x[i][j].getsup())
		std::cout << true << ' ';
	    else
		std::cout << false << ' ';
	}
	std::cout << '\n';
    }
    
    engine.release_flows();
    
    return 0;
}
