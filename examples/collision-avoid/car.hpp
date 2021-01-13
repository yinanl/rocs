/**
 *  car.hpp
 *
 *  Kinematics of a car-like mobile robot
 *
 *  Created by Yinan Li on Feb. 21, 2017.
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _robotcar_h
#define _robotcar_h


#include <vector>
#include <cmath>


/* user defined dynamics */
const double tau = 0.3;  // sampling time

struct carde {

    static const int n = 3;  // system dimension
    static const int m = 2;

    /**
     * Discrete-time dynamics
     * @param h sampling time
     * @param x system state: [x,y,theta], n=3
     * @param u control array (size of 2, velocity and steering angle)
     * @param nu the number of different control values
     */
    template<typename S>
    carde(S &dx, const S &x, rocs::Rn u) {
	if (std::fabs(u[0]) < 1e-6) { //v=0
	    dx = x;
	} else if (std::fabs(u[1]) < 1e-6) { //w=0
	    dx[0] = x[0] + u[0]* cos(x[2])*tau;
	    dx[1] = x[1] + u[0]* sin(x[2])*tau;
	    dx[2] = x[2];
	} else { //v,w not 0
	    dx[0] = x[0] + u[0]/u[1]*2*sin(u[1]*tau/2.)*cos(x[2]+u[1]*tau/2.);
	    dx[1] = x[1] + u[0]/u[1]*2*sin(u[1]*tau/2.)*sin(x[2]+u[1]*tau/2.);
	    dx[2] = x[2] + u[1] * tau;
	}
    }
    
}; // struct carde

#endif
