/**
 *  testReachSet.cpp
 *
 *  Test abstraction by using a car kinematics model.
 *
 *  Created by Yinan Li on Nov. 25, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>

#include "src/grid.h"
#include "src/definitions.h"
#include "src/abstraction.hpp"
#include "src/hdf5io.h"


const double h = 0.3;  // sampling time
const double dt = 0.001; //integration step size for odeint

/* For reachable set computation */
struct twoagent {
    static const int n = 3;  // system dimension
    static const int nu = 2;  // control dimension
    rocs::ivec d{rocs::interval(-0.8, 0.8),
		 rocs::interval(-0.8, 0.8)};

    /* template constructor
     * @param[out] dx
     * @param[in] x = [xr, yr, psir]
     * @param u = [v, w]
     * @param d = [v', w']
     */
    template<typename S>
    twoagent(S *dx, const S *x, rocs::Rn u) {
	dx[0] = -u[0] + d[0]*cos(x[2]) + u[1]*x[1];
	dx[1] = d[0]*sin(x[2]) - u[1]*x[0];
	dx[2] = d[1] - u[1];
    }
};

/* For solving ode */
struct car_ode {
    rocs::Rn u;
    rocs::Rn d;
    car_ode (const rocs::Rn p1, const rocs::Rn p2): u(p1), d(p2){}
    /**
     * ODE model
     * @param x system state: [x,y,theta], n=3
     * @param dxdt vector field
     * @param t time
     */
    void operator() (rocs::Rn &x, rocs::Rn &dxdt, double t) const
    {
	dxdt[0] = -u[0] + d[0]*cos(x[2]) + u[1]*x[1];
	dxdt[1] = d[0]*sin(x[2]) - u[1]*x[0];
	dxdt[2] = d[1] - u[1];
    }
};


int main()
{
    /* Config */
    clock_t tb, te;
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;

    /**
     * Case I
     */
    /* Set the state and control space */
    const int xdim = 3;
    const int udim = 2;
    
    double xlb[] = {-3, -3, -M_PI};
    double xub[] = {3, 3, M_PI};
    double eta[] = {0.2, 0.2, 0.2};
    
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /**
     * Define the two-agent system
     */
    double t = 0.3;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    rocs::CTCntlSys<twoagent> safety("collision-free", t,
    				     twoagent::n, twoagent::nu,
    				     delta, &controlparams);

    safety.init_workspace(xlb, xub);
    safety.init_inputset(mu, ulb, uub);
    safety.allocate_flows();

    rocs::abstraction< rocs::CTCntlSys<twoagent> > abst(&safety);
    abst.init_state(eta, xlb, xub);
    std::cout << "# of in-domain nodes: " << abst._x._nv << '\n';


    /* Test reachable set computation */
    int suc = 1;
    size_t na = safety._ugrid._nv;
    rocs::Rn x_test{-1.4, 0.6, -1.9};
    // rocs::Rn x_test{1.4, 0.6, -1.9};
    rocs::ivec x0 = {rocs::interval(x_test[0]-eta[0]/2., x_test[0]+eta[0]/2.),
		     rocs::interval(x_test[1]-eta[1]/2., x_test[1]+eta[1]/2.),
		     rocs::interval(x_test[2]-eta[2]/2., x_test[2]+eta[2]/2.)};
    std::cout << "Checking rechable set computation for " << x0 << "...\n";
    rocs::Rn u(udim);
    rocs::Rn xpost(xdim);
    rocs::ivec box(xdim);
    std::vector<rocs::ivec> reachset(na, rocs::ivec(xdim));
    rocs::Rn corner(xdim), p0(xdim);
    int quo, rem;
    rocs::ivec margin{rocs::interval(-rocs::EPSIVAL, rocs::EPSIVAL),
    		      rocs::interval(-rocs::EPSIVAL, rocs::EPSIVAL),
    		      rocs::interval(-rocs::EPSIVAL, rocs::EPSIVAL)};
    rocs::ivec yt(xdim);
    
    safety.get_reach_set(reachset, x0);
    std::vector<rocs::Rn> d(4, rocs::Rn {0.7, 0.7});
    d[1][0] = -0.7;
    d[2][0] = -0.7;
    d[2][1] = -0.7;
    d[3][1] = -0.7;
    /* Test valid control inputs */
    for(size_t j = 0; j < na; ++j) {
	safety._ugrid.id_to_val(u, j); //get control values
	for(auto de:d) {
	    /* Test if the reachable set covers ode solutions of all corners */
	    for(int k = 0; k < std::pow(2, xdim); ++k) {
		quo = k;
		for(int l = 0; l < xdim; ++l) {
		    if(quo % 2) {
			corner[l] = x_test[l]+eta[l]/2.; //upper bound
			p0[l] = corner[l];
		    } else {
			corner[l] = x_test[l]-eta[l]/2.; //lower bound
			p0[l] = corner[l];
		    }
		    quo /= 2;
		}
		// std::cout << "Corner "
		// 	  << '(' << corner[0] << ',' << corner[1] << ',' << corner[2] << ")\n";
		boost::numeric::odeint::integrate_const(rk45, car_ode(u,de), corner, 0.0, h, dt);
		yt = reachset[j] + margin;
		if(!yt.isin(corner)) {
		    std::cout << "Incorrect at u=" << j
			      << '(' << u[0] << ',' << u[1] << ')'
			      << ", d=" << '(' << de[0] << ',' << de[1] << ") for "
			      << '(' << p0[0] << ',' << p0[1] << ',' << p0[2] << "):\n"
			      << '(' << corner[0] << ',' << corner[1] << ',' << corner[2] << ')'
			      << " is not in " << yt << '\n';
		    suc = 0;
		}
	    }
	}
    }//end for control values

    if(suc)
	std::cout << "The test point passes the test.\n";
    
    return suc;
}
