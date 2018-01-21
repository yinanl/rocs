/**
 *  ipdl.hpp
 *
 *  Inverted pendulum [Michigan control tutorial](http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling)
 *
 *  Created by yinan li on April 26, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _invertpdl_h
#define _invertpdl_h


#include <cmath>
#include <cassert>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <armadillo>
#include "src_itvl/vectorfield.h"




const double m = 0.2;
const double b = 0.1;
const double J = 0.006;
const double g = 9.8;
const double l = 0.3;

const double a1 = (m*g*l) / (J + m*l*l);
const double a2 = b / (J + m*l*l);
const double a3 = l / (J + m*l*l);

const double abs_tol = 1E-12;
const double rel_tol = 1E-12;


/* continuous-time dynamics */
struct ipdl_vf {

    double _u;

    ipdl_vf(const double param): _u(param) {}

    void operator()(rocs::state_type &x, rocs::state_type &dxdt, double t) const
    {
	dxdt[0] = x[1];
	dxdt[1] = a1*std::sin(x[0]) - a2*x[1] + a3*std::cos(x[0])*_u;
    }

};


/** 
 * Point flow map
 * The original dynamics is given in continuous-time
 * we use time-T map to change the dynamics into discrete-time.
 * @param x [inout] an inital interval vector x0
 * @param u a list of control inputs
 * @param T integrate time horizon
 * @param dt integrate time step
 */
void ipdl_tmap(rocs::state_type &x, const double u, const double T, const double dt) {

    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::state_type> rk45;
  
    boost::numeric::odeint::integrate_const(rk45, ipdl_vf(u), x, 0.0, T, dt);
  
    /* rk45.do_step(ipdl(u), x, 0.0, T); */
}



/**
 * Interval flow map
 * @param x an inital interval vector x0
 * @param u a list of control inputs
 * @param h integrate time horizon
 * @param dt integrate time step
 * @return a list of end-state interval xt w.r.t. different u
 */
std::vector<rocs::ivec> ipdl_tmapbox(const rocs::ivec &x,
				     const rocs::input_type &u,
				     const double h, const double dt) {

    /* growth bound matrix L(u) \cite{Gunther2016} */
    arma::mat Lu = {{0, 1}, {std::sqrt(a1*a1 + a3*a3*u[0][0]*u[0][0]), -a2}};
    
    arma::vec r(x.radius());
    std::vector<double> xc = x.mid();
    
    // std::vector<ivec> yu(u.size());
    // ivec y(x.getdim());
    // arma::vec gb;
    
    // for (int i = 0; i < u.size(); ++i) {
	
    // 	Lu(1,0) = sqrt(a1*a1 + a3*a3*u[i][0]*u[i][0]);
    // 	gb = expmat(Lu*h) * r;
    // 	xc = x.mid();
    // 	ipdl_tmap(xc, u[i][0], h, dt);

    // 	y.setval(0, interval(xc[0]-gb[0], xc[0]+gb[0]));
    // 	y.setval(1, interval(xc[1]-gb[1], xc[1]+gb[1]));
    // 	yu[i] = y;
    // }

    std::vector<rocs::ivec> yu(u.size(), rocs::ivec(x.getdim()));
    arma::vec gb;
    
    for (int i = 0; i < u.size(); ++i) {
	
	Lu(1,0) = std::sqrt(a1*a1 + a3*a3*u[i][0]*u[i][0]);
	gb = arma::expmat(Lu*h) * r;
	xc = x.mid();
	ipdl_tmap(xc, u[i][0], h, dt);

	yu[i].setval(0, rocs::interval(xc[0]-gb[0], xc[0]+gb[0]));
	yu[i].setval(1, rocs::interval(xc[1]-gb[1], xc[1]+gb[1]));
    }
    
    return yu;
}


/**
 * Derived class for inverted pendulum.
 * 
 */
class ipdl : public rocs::VFunctor {
public:

    ipdl();
    ipdl(double T, double dt): VFunctor(T, dt) {}
    ipdl(rocs::input_type &U, double T, double dt): VFunctor(U, T, dt) {}
    ipdl(rocs::grid ug, double T, double dt) : VFunctor(ug, T, dt) {}
    ipdl(const int dim, const double lb[], const double ub[],
	 const double mu[], double T, double dt) : VFunctor(dim, lb, ub, mu, T, dt) {}

    /* override operator () */
    virtual std::vector<rocs::ivec> operator()(const rocs::ivec &x) {

	assert(_tau > 0);

	if (!_ugrid._data.empty())
	    return ipdl_tmapbox(x, _ugrid._data, _tau, _dt);
	else {
	    std::cout << "Callback of inverted pendulum vectorfield failed: no input data.\n";
	    assert(false);
	}
    }


    virtual ~ipdl() {}
};


#endif
