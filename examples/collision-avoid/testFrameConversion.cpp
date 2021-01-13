/**
 *  testFrameConversion.cpp
 *
 *  Test the frame conversion between inertial and body frames.
 *
 *  Created by Yinan Li on Nov. 27, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <cstdlib>
#include <ctime>


#include "src/definitions.h"



/* define dynamics */
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


int main() {
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;
    rocs::Rn z0{5.08561,10.1089,1.8}, z1{3.7, 10.8, 0.0};
    rocs::Rn u{0.9, -0.9}, d{0.6, 0.6/4.61};
    rocs::Rn zr(3);
    for(int i = 0; i < 3; ++i) {
	zr[i] = z1[i]-z0[i];
    }
    std::cout << "z1 position to z0 in inertial frame initially: ";
    for(int i = 0; i < 3; ++i) {
	std::cout << zr[i] << ' ';
    }
    std::cout << '\n';
    
    inertia_to_body(zr, z0);
    std::cout << "z1 position to z0 in z0 body frame initially: ";
    for(int i = 0; i < 3; ++i) {
	std::cout << zr[i] << ' ';
    }
    std::cout << '\n';
    
    const double h = 0.3;  // sampling time
    const double dt = 0.001; //integration step size for odeint
    /* Integrate z0 and z1 separately */
    boost::numeric::odeint::integrate_const(rk45, car_dynamics(u), z0, 0.0, h, dt);
    boost::numeric::odeint::integrate_const(rk45, car_dynamics(d), z1, 0.0, h, dt);
    std::cout << "z1-z0 = ";
    for(int i = 0; i < 3; ++i) {
	std::cout << z1[i]-z0[i] << ' ';
    }
    std::cout << '\n';
    /* Integrate by the local model */
    boost::numeric::odeint::integrate_const(rk45, twoagent_ode(u,d), zr, 0.0, h, dt);
    std::cout << "z1 position to z0 in z0 body frame after 0.3s: ";
    for(int i = 0; i < 3; ++i) {
	std::cout << zr[i] << ' ';
    }
    std::cout << '\n';
    
    body_to_inertia(zr, z0);
    std::cout << "z1 position to z0 in inertial frame after 0.3s: ";
    for(int i = 0; i < 3; ++i) {
	std::cout << zr[i] << ' ';
    }
    std::cout << '\n';
}
