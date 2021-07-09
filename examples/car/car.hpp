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
const double h = 0.3;  // sampling time
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
	double alpha, at, r;
	alpha = atan(tan(u[1]) / 2);

	if (std::fabs(u[0]) < 1e-6) {
	    dx = x;
	} else if (std::fabs(u[1]) < 1e-6) {
	    dx[0] = x[0] + u[0]* cos(x[2])*h;
	    dx[1] = x[1] + u[0]* sin(x[2])*h;
	    dx[2] = x[2];
	} else {
	    at = alpha + u[0] * tan(u[1]) * h / 2;
	    r = 2 * sin(u[0]*tan(u[1])*h/2) / (cos(alpha) * tan(u[1]));

	    dx[0] = x[0] + r * cos(at + x[2]);
	    dx[1] = x[1] + r * sin(at + x[2]);
	    dx[2] = x[2] + u[0] * tan(u[1]) * h;
	}
    }

}; // struct carde

#endif
