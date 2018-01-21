/**
 *  dcdc.h
 *
 *  A header file of dcdc converter interval extensions.
 *
 *  Created by yinan li on Nov. 17, 2016.
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _dcdc_h
#define _dcdc_h


#include <cmath>
#include <cassert>
#include <armadillo>
#include "src_itvl/vectorfield.h"


const double xc = 70.0;
const double xl = 3.0;
const double rc = 0.005;
const double rl = 0.05;
const double r0 = 1.0;
const double vs = 1.0;

const double a11 = - rl / xl;
const double a14 = -1 / (xc * (rc + r0));
const double b1 = vs / xl;

const double a21 = (-1/xl) * (rl + r0*rc / (r0+rc));
const double a22 = (-1/xl) * (r0/(r0+rc));
const double a23 = (1/xc) * (r0/(r0+rc));
const double a24 = (-1/xc) * (1/(r0+rc));

const double TS = 0.5;
const size_t M = 2;


/**
 * 2-mode boost DC-DC converter:
 * discrete-time interval model
 */
std::vector<rocs::ivec> dcdc(const rocs::ivec &x) {

    int n = x.getdim();
    assert(n == 2);

    std::vector<rocs::ivec> y(M);
    arma::mat I = arma::eye<arma::mat>(n, n);
    arma::vec b = {b1, 0};
    
    /* mode 1 (y[0]): decoupled */
    // x1(t+h)=e^(a11 h)x1(t) + b1/a11(e^(a11 h)-1)
    // x2(t+h)=e^(a14 h)x2(t)
    arma::mat A1 = {{a11, 0}, {0, a14}};
    arma::mat F1 = arma::expmat(A1 * TS);
    arma::vec g1 = arma::inv(A1) * (F1 - I) * b;
    y[0] = linmap(F1, g1, x);

    
    /* mode 2: use the solution form and zonotope operations */
    // x(t+h) = e^(Ah)x(t) + A^-1(e^(Ah)-I)b
    arma::mat A2 = { {a21, a22}, {a23, a24} };
    arma::mat F2 = arma::expmat(A2 * TS);
    arma::vec g2 = arma::inv(A2) * (F2 - I) * b;

    y[1] = linmap(F2, g2, x);
    
    return y;
}


// /* dummy implementation */
// std::vector<ivec> dcdc_dummy(const ivec &x) {

//     int n = x.getdim();
//     assert(n == 2);

//     std::vector<ivec> y(M);
//     arma::mat I = arma::eye<arma::mat>(n, n);
    
//     /* mode 1 (y[0]): decoupled */
//     // x1(t+h)=e^(a11 h)x1(t) + b1/a11(e^(a11 h)-1)
//     // x2(t+h)=e^(a14 h)x2(t)
//     ivec y1(n);
//     y1[0] = interval(exp(a11*TS)*x[0].getinf()+b1/a11*(exp(a11*TS)-1),
// 		     exp(a11*TS)*x[0].getsup()+b1/a11*(exp(a11*TS)-1) );
//     y1[1] = interval(exp(a14*TS)*x[1].getinf(), exp(a14*TS)*x[1].getsup());
//     y[0] = y1;

//     /* mode 2: use the solution form and zonotope operations */
//     // x(t+h) = e^(Ah)x(t) + A^-1(e^(Ah)-I)b
//     arma::mat A2 = { {a21, a22}, {a23, a24} };
//     arma::mat F2 = expmat(A2 * TS);
//     arma::vec b = {b1, 0};
//     arma::vec g2 = inv(A2) * (F2 - I) * b;
    
//     std::vector<double> c = x.mid();
//     std::vector<double> r = x.radius();

//     double yr1, yr2;
//     yr1 = fabs( F2(0,0)*r[0] ) + fabs( F2(0,1)*r[1] );
//     yr2 = fabs( F2(1,0)*r[0] ) + fabs( F2(1,1)*r[1]);
    
//     ivec y2(n);
//     y2[0] = interval(F2(0,0)*c[0]+F2(0,1)*c[1] - yr1 + g2[0], F2(0,0)*c[0]+F2(0,1)*c[1] + yr1 + g2[0]);
//     y2[1] = interval(F2(1,0)*c[0]+F2(1,1)*c[1] - yr2 + g2[1], F2(1,0)*c[0]+F2(1,1)*c[1] + yr2 + g2[1]);

//     y[1] = y2;

//     return y;
// }




std::vector<rocs::ivec> dcdc_vf(const rocs::ivec &x, double tau) {

    int n = x.getdim();
    assert(n == 2);

    std::vector<rocs::ivec> y(2);
    arma::mat I = arma::eye<arma::mat>(n, n);
    arma::vec b = {b1, 0};
    
    /* mode 1 (y[0]): decoupled */
    // x1(t+h)=e^(a11 h)x1(t) + b1/a11(e^(a11 h)-1)
    // x2(t+h)=e^(a14 h)x2(t)
    arma::mat A1 = {{a11, 0}, {0, a14}};
    arma::mat F1 = arma::expmat(A1 * tau);
    arma::vec g1 = arma::inv(A1) * (F1 - I) * b;
    y[0] = linmap(F1, g1, x);

    
    /* mode 2: use the solution form and zonotope operations */
    // x(t+h) = e^(Ah)x(t) + A^-1(e^(Ah)-I)b
    arma::mat A2 = { {a21, a22}, {a23, a24} };
    arma::mat F2 = arma::expmat(A2 * tau);
    arma::vec g2 = arma::inv(A2) * (F2 - I) * b;

    y[1] = linmap(F2, g2, x);
    
    return y;
    
}


/**
 * Functor for dcdc converter vector field:
 *
 * to be consistent with the vf without controls
 * */

class DCDC : public rocs::VFunctor {
    
public:

    /* constructors */
    DCDC() {}
    
    DCDC(rocs::input_type &modes, double h): VFunctor(modes, h) {}
    
    
    /* override operator () */
    virtual std::vector<rocs::ivec> operator()(const rocs::ivec &x) {

	return dcdc_vf(x, _tau);
    }


    virtual ~DCDC() {}
};


#endif
