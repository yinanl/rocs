/**
 *  mgode.hpp
 *
 *  ODEs of Moore-Greitzer engine
 *
 *  Created by Yinan Li on June 11, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _mgode_h
#define _mgode_h

double a = 1.67;
double B = 2.0;
double H = 0.18;
double W = 0.25;
double lc = 8.0;
double cx = 1.0/lc;
double cy = 1.0/(4*lc*B*B);
double aH = a+H;
double H2 = H/(2.0*W*W*W);
double W2 = 3*W*W;

/* user defined dynamics */
struct mgode {

    static const int n = 3;  // system dimension
    static const int nu = 2;  // control dimension
    
    /* template constructor
     * @param[out] dx
     * @param[in] x
     * @param u
     */
    template<typename S>
    mgode(S *dx, const S *x, rocs::Rn u) {
	dx[0] = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]) + u[0];
	dx[1] = cy * (x[0] - x[2]*sqrt(x[1]));
	dx[2] = u[1];
    }
};

#endif
