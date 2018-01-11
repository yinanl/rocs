/**
 *  car.hpp
 *
 *  Dynamics of 2d car [example from \cite{Gunther2016}]
 *  - dynamics
 *  - constraints
 *  - datasets
 *
 *  Created by yinan li on Feb. 21, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _robotcar_h
#define _robotcar_h


#include <vector>
#include <cmath>

#include "../../src_itvl/vectorfield.h"


const int XD = 3;
const int UD = 2;

/**
 * Discrete-time dynamics
 * @param h sampling time
 * @param x system state: [x,y,theta], n=3
 * @param u control array (size of 2, velocity and steering angle)
 * @param nu the number of different control values
 */
template <class XT, class UT>
std::vector<XT> carvf(const XT &x, const UT &u, const double h)
{
    double alpha, at, r;

    /* define output array of interval vectors */
    std::vector<XT> y(u.size());

    XT dx(3);
    for (int i = 0; i < u.size(); ++i) {

	alpha = atan(tan(u[i][1]) / 2);

	if (fabs(u[i][0]) < 1e-6) {
      
	    y[i] = x;
	}
	else if (fabs(u[i][1]) < 1e-6) {

	    dx[0] = x[0] + u[i][0]*cos(x[2])*h;
	    dx[1] = x[1] + u[i][0]*sin(x[2])*h;
	    dx[2] = x[2];

	    y[i] = dx;
	}
	else {

	    at = alpha + u[i][0] * tan(u[i][1]) * h / 2;
	    r = 2 * sin(u[i][0]*tan(u[i][1])*h/2) / (cos(alpha) * tan(u[i][1]));
      
	    dx[0] = x[0] + r * cos(at + x[2]);

	    dx[1] = x[1] + r * sin(at + x[2]);

	    dx[2] = x[2] + u[i][0] * tan(u[i][1]) * h;

	    y[i] = dx;
	}
    
    }

    return y;
}


/**
 * Derived functor for 2d car dynamics.
 */
class car : public VFunctor {
    
public:

    /* constructors */
    car() {}
    car(double h) : VFunctor(h) {}
    car(input_type &u, double h): VFunctor(u, h) {}
    car(grid ug, double h) : VFunctor(ug, h) {}
    car(const int dim, const double lb[], const double ub[],
	const double mu[], double h) : VFunctor(dim, lb, ub, mu, h) {}
    
    
    /* override operator () */
    virtual std::vector<ivec> operator()(const ivec &x) {

	if (!_uset.empty())
	    return carvf(x, _uset, _tau);
	else if (!_ugrid._data.empty())
	    return carvf(x, _ugrid._data, _tau);
	else
	    assert(false);
    }

    virtual std::vector<state_type> operator()(const state_type &x) {

	if (!_uset.empty())
	    return carvf(x, _uset, _tau);
	else if (!_ugrid._data.empty())
	    return carvf(x, _ugrid._data, _tau);
	else
	    assert(false);
    }


    virtual ~car() {}
};




#endif
