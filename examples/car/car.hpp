/**
 *  car.hpp
 *
 *  Dynamics of 2d car
 *  - dynamics
 *  - constraints
 *  - datasets
 *
 *  Created by Yinan Li on Feb. 21, 2017.
 *  Revised by Yinan Li on May 12, 2018.
 *
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
	}
	else if (std::fabs(u[1]) < 1e-6) {

	    dx[0] = x[0] + u[0]* cos(x[2])*h;
	    dx[1] = x[1] + u[0]* sin(x[2])*h;
	    dx[2] = x[2];
	}
	else {
	    at = alpha + u[0] * tan(u[1]) * h / 2;
	    r = 2 * sin(u[0]*tan(u[1])*h/2) / (cos(alpha) * tan(u[1]));
      
	    dx[0] = x[0] + r * cos(at + x[2]);
	    dx[1] = x[1] + r * sin(at + x[2]);
	    dx[2] = x[2] + u[0] * tan(u[1]) * h;
	}
    }
    
}; // struct carde


/*** constraints ***/
const double H = 1.0; // car width
const double L = 2.0; // car length
const double d = 2.0; // marginal parking space
const double D = 0.5; // distance to the curb

/* rear and front car corners collide with the car body */
template<typename T>
T cst_v5(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}
template<typename T>
T cst_v6(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]-L) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]-L) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]-L) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]-L) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}
template<typename T>
T cst_v7(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]-(2*L+d)) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]-(2*L+d)) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]-(2*L+d)) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]-(2*L+d)) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}
template<typename T>
T cst_v8(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]-(3*L+d)) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]-(3*L+d)) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]-(3*L+d)) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]-(3*L+d)) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}

/* the car body collides with the rear car area: consider 4 corners */
template<typename T>
T cst_v1rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(H/2.0) - H;
    return y;
}
template<typename T>
T cst_v2rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(-H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(-H/2.0) - H;
    return y;
}
template<typename T>
T cst_v3rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(-H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(-H/2.0) - H;
    return y;
}
template<typename T>
T cst_v4rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(H/2.0) - H;
    return y;
}

/* the car body collides with the front car area: consider 4 corners */
template<typename T>
T cst_v1front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(H/2.0) - (3.0*L+d);
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(L/2.0) + sin(x[2])*(H/2.0) + (2.0*L+d);
    return y;
}
template<typename T>
T cst_v2front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(-H/2.0) - (3*L+d);
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(-H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(L/2.0) + sin(x[2])*(-H/2.0) + (2.0*L+d);
    return y;
}
template<typename T>
T cst_v3front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(-H/2.0) - (3*L+d);
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(-H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(-L/2.0) + sin(x[2])*(-H/2.0) + (2.0*L+d);
    return y;
}
template<typename T>
T cst_v4front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(H/2.0) - (3*L+d);
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(-L/2.0) + sin(x[2])*(H/2.0) + (2.0*L+d);
    return y;
}

/* the car body collides with the curb area */
template<typename T>
T cst_v1curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(H/2.0) + D;
    return y;
}
template<typename T>
T cst_v2curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(-H/2.0) + D;
    return y;
}
template<typename T>
T cst_v3curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(-H/2.0) + D;
    return y;
}
template<typename T>
T cst_v4curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(H/2.0) + D;
    return y;
}


#endif
