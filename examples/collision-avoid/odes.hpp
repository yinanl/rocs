/**
 *  odes.hpp
 *
 *  ODEs of a car-like mobile robot and a two-agent system.
 *
 *  Created by Yinan Li on Feb. 12, 2021.
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _odes_h
#define _odes_h


#include <vector>
#include <cmath>

#include "src/definitions.h"


/* The kinematics of a mobile robot in inertia frame */
struct car_dynamics {
    rocs::Rn u;
    car_dynamics (const rocs::Rn param): u (param) {}

    void operator() (rocs::Rn &x, rocs::Rn &dxdt, double t) const
    {
	dxdt[0] = u[0]*std::cos(x[2]);
	dxdt[1] = u[0]*std::sin(x[2]);
	dxdt[2] = u[1];
    }
};

/* Kinematics of the other robot in the body frame of an ego robot */
struct twoagent_ode {
    rocs::Rn u;
    rocs::Rn d;
    twoagent_ode (const rocs::Rn p1, const rocs::Rn p2): u(p1), d(p2){}
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


/* Frame conversions */
void inertia_to_body(rocs::Rn &z, rocs::Rn z0) {
    double x = std::cos(z0[2])*z[0] + std::sin(z0[2])*z[1];
    double y = -std::sin(z0[2])*z[0] + std::cos(z0[2])*z[1];
    z[0] = x;
    z[1] = y;
}
void body_to_inertia(rocs::Rn &z, rocs::Rn z0) {
    double x = std::cos(z0[2])*z[0] - std::sin(z0[2])*z[1];
    double y = std::sin(z0[2])*z[0] + std::cos(z0[2])*z[1];
    z[0] = x;
    z[1] = y;
}

#endif
