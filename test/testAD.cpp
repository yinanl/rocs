/**
 * Example of error controlled Taylor model flowpipe computation:
 * Reversed Van de Pol system
 */

#include <iostream>
#include <cmath>

#include "src/flow_taylor.hpp"
#include "src/timer.hpp"


/* user defined dynamics (a template functor) */
struct vdpode {
    static const int n = 2;
    
    /* template constructor
     * @param[out] dx
     * @param[in] x
     * @param u
     */
    template<typename S>
    vdpode(S *dx, const S *x) {
	dx[0] = -x[1];
	dx[1] = x[0] + x[1]*(x[0]*x[0]-1);
    }
};



int main() {

    /* Set sampling time and disturbance */
    double t = 0.05;  // fixed sampling time
    double delta = 10;

    /* Set parameters for computation */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    
    rocs::flowTaylor<vdpode, double, bool> vdpreal(&controlparams, t, delta);

    /* initial state and interval */
    double rx0[2] = {0.5, 0.3};
    
    vdpreal.reset_taylor_coeffs();
    vdpreal.init_taylor_coeffs(rx0);
    vdpreal.eval_taylor_coeffs(10);
    vdpreal.print_taylor_coeffs();

    return 0;
}
