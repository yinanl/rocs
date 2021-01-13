/**
 *  testReachset.cpp
 *
 *  Test reachable set computation for simplified car kinematics model:
 *  compare two ways of computing reachable set.
 *
 *  Created by Yinan Li on Oct. 29, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>
#include <array>
#include <boost/numeric/odeint.hpp>

#include "RungeKutta4.hh"
#include "src/system.hpp"
#include "src/csolver.h"

#include "car.hpp"

/* state space dim */
const int state_dim=3;
/* input space dim */
const int input_dim=2;
/* sampling time */
const double tau = 0.3;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* we integrate the vehicle ode by tau sec (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, const input_type &u) {
    /* the ode describing the vehicle */
    auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) {
        xx[0] = u[0]*std::cos(x[2]);
        xx[1] = u[0]*std::sin(x[2]);
        xx[2] = u[1];
    };
    /* simulate (use 10 intermediate steps in the ode solver) */
    scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,10);
};

/* we integrate the growth bound by 0.3 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, const state_type &, const input_type &u) {
    double c = std::abs(u[0]);
    r[0] = r[0]+c*r[2]*tau;
    r[1] = r[1]+c*r[2]*tau;
};


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


int main()
{
    // /* set the state space */
    // const double theta = 3.5;
    // double xlb[] = {0, 0, -theta};
    // double xub[] = {10, 10, theta};
    // const double eta[] = {.2,.2,.2}; /* set precision */
    
    // /* set the control values */
    // double ulb[] = {-1.0, -1.0};
    // double uub[] = {1.0, 1.0};
    // double mu[] = {0.3, 0.3};

    // /* define the control system */
    // rocs::DTCntlSys<carde> car("reach goal", h, carde::n, carde::m);
    // car.init_workspace(xlb, xub);
    // car.init_inputset(mu, ulb, uub);


    // /* Test */
    // state_type x_test{-1.4, 0.6, -1.9};

    // /* interval method */
    // rocs::ivec x0 = {rocs::interval(x_test[0]-eta[0]/2., x_test[0]+eta[0]/2.),
    // 		     rocs::interval(x_test[1]-eta[1]/2., x_test[1]+eta[1]/2.),
    // 		     rocs::interval(x_test[2]-eta[2]/2., x_test[2]+eta[2]/2.)};
    // std::vector<rocs::ivec> y0(car._ugrid._nv, rocs::ivec(3));
    // car.get_reach_set(y0, x0);

    // std::cout << "The initial interval: " << x0 << '\n';
    // std::cout << "The integrating time: " << tau << '\n';
    // for(size_t i = 0; i < car._ugrid._nv; ++i) {
    // 	/* nominal trajectory + growth bound */
    // 	state_type r{0.1, 0.1, 0.1};
    // 	state_type x = x_test;
    // 	input_type u{car._ugrid._data[i][0], car._ugrid._data[i][1]};
    // 	radius_post(r, x, u);
    // 	vehicle_post(x, u);
    // 	rocs::ivec y = {rocs::interval(x[0]-r[0], x[0]+r[0]),
    // 			rocs::interval(x[1]-r[1], x[1]+r[1]),
    // 			rocs::interval(x[2]-r[2], x[2]+r[2])};
    // 	std::cout << "R(t, x0, [" << car._ugrid._data[i][0] << ','
    // 		  << car._ugrid._data[i][1] << "])= ";
    // 	std::cout << y << "(gb), or " << y0[i] << "(interval).\n";
    // }


    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;
    rocs::Rn z0{5.08561,10.1089,1.8}, z1{3.7, 10.8, 0.0};
    rocs::Rn zr(3);
    for(int i = 0; i < 3; ++i) {
	zr[i] = z1[i]-z0[i];
    }
    rocs::Rn u{0.9, -0.9}, d{0.6, 0.6/4.61};

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
    std::cout << "zr = ";
    for(int i = 0; i < 3; ++i) {
	std::cout << zr[i] << ' ';
    }
    std::cout << '\n';

    
    return 0;
}
