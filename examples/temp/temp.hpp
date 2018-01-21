/**
 *  temp.hpp
 *
 *  A header file defines temperature control system dynamics.
 *
 *  Created by yinan li on August 14, 2017.
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _temperature_h
#define _temperature_h


#include <cmath>
#include <cassert>
#include "src_itvl/vectorfield.h"


const double r1 = 0.002;
const double c = 16;
const double r2 = 0.1;


/*
 * 4-mode room temperature control:
 * discrete-time interval model
 */
std::vector<rocs::ivec> tpc_vf(const rocs::ivec &x, double tau) {

    int n = x.getdim();
    assert(n == 2);

    std::vector<rocs::ivec> y(4, rocs::ivec(2));
    
    /* mode 1: off */
    y[0].setval(0, std::exp(-r1*tau)*x[0] + (1.0 - std::exp(-r1*tau))*c);
    y[0].setval(1, x[1]);
    
    /* mode 2: heating */
    y[1].setval(0, std::exp(-r1*tau)*x[0] + (1.0 - std::exp(-r1*tau))*(x[1]-r2/r1) + r2*tau);
    y[1].setval(1, x[1] + r2*tau);

    /* mode 3: cooling */
    y[2].setval(0, std::exp(-r1*tau)*x[0] + (1.0 - std::exp(-r1*tau))*(x[1]+r2/r1) - r2*tau);
    y[2].setval(1, x[1] - r2*tau);

    /* mode 4: on */
    y[3].setval(0, std::exp(-r1*tau)*x[0] + (1.0 - std::exp(-r1*tau))*x[1]);
    y[3].setval(1, x[1]);
    
    return y;
    
}


/* functor for temperature vector field:
 *
 * to be consistent with the vf without controls
 * */

class TPC : public rocs::VFunctor {
    
public:

    /* constructors */
    TPC() {}
    
    TPC(rocs::input_type &modes, double h): VFunctor(modes, h) {}
    
    
    /* override operator () */
    virtual std::vector<rocs::ivec> operator()(const rocs::ivec &x) {

	return tpc_vf(x, _tau);
    }


    virtual ~TPC() {}
};


#endif
