/**
 * Example of error controlled Taylor model flowpipe computation:
 * Reversed Van de Pol system
 */

#include <fstream>
#include <iomanip>
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

    std::ofstream file;
    file.open("reachset.txt", std::ios::out | std::ios::app);
    
    /* Set sampling time and disturbance */
    double t = 0.05;  // fixed sampling time
    double delta = 10;

    /* Set parameters for computation */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    
    /* create the Taylor model */
    rocs::flowTaylor<vdpode, rocs::interval, bool> vdp(&controlparams, t, delta);
    std::cout << "log(1-alpha)delta= " << vdp._logdel1 << '\n';
    std::cout << "log(tau)= " << vdp._logdel2 << '\n';

    double epsilon = 0.01;
    /* Case 1:
     * small [x0] 
     */
    double w = 0.015;
    
    double rx0[2] = {0.5, 0.3};
    rocs::ivec ix0(vdpode::n);
    ix0[0] = rocs::interval(rx0[0]-w/2,rx0[0]+w/2);
    ix0[1] = rocs::interval(rx0[1]-w/2,rx0[1]+w/2);
    std::cout << "Compute valid reachable set for " << ix0 << ":\n";
    // std::cout << "Computing bound for Taylor terms: " << '\n'
    // 	      << vdp.eval_taylorterm_bound(ix0);
    rocs::ivec ixt(vdpode::n);  // define output reachable set
    timer T;
    T.reset();
    std::cout << vdp.reachset_robust(ixt, ix0, epsilon) << ": ";
    std::cout << "Time taken: " << T.elapsed() << '\n';
    std::cout << ixt << std::endl;
    std::cout << vdp._kbar << std::endl;
    std::cout << rocs::eval_epsilon(vdp) << '\n';

    file << std::fixed << std::setprecision(8) << ixt << '\n';
    
    double rx1[2] = {-1.2, -0.6};
    w = 0.0035;
    rocs::ivec ix1(vdpode::n);
    ix1[0] = rocs::interval(rx1[0]-w/2,rx1[0]+w/2);
    ix1[1] = rocs::interval(rx1[1]-w/2,rx1[1]+w/2);
    std::cout << "Compute valid reachable set for " << ix1 << ":\n";
    std::cout << vdp.reachset_robust(ixt, ix1, epsilon) << ": ";
    std::cout << ixt << std::endl;
    std::cout << vdp._kbar << std::endl;
    std::cout << rocs::eval_epsilon(vdp) << '\n';
    
    file << std::fixed << std::setprecision(8) << ixt << '\n';

    double rx2[2] = {-2.3, 1.7};
    w = 0.0002;
    rocs::ivec ix2(vdpode::n);
    ix2[0] = rocs::interval(rx2[0]-w/2,rx2[0]+w/2);
    ix2[1] = rocs::interval(rx2[1]-w/2,rx2[1]+w/2);
    std::cout << "Compute valid reachable set for " << ix2 << ":\n";
    std::cout << vdp.reachset_robust(ixt, ix2, epsilon) << ": ";
    std::cout << ixt << std::endl;
    std::cout << vdp._kbar << std::endl;
    std::cout << rocs::eval_epsilon(vdp) << '\n';
    
    file << std::fixed << std::setprecision(8) << ixt << '\n';
    

    // /* Case 2: 
    //  * [x0]=X (the state space) 
    //  */
    // rocs::ivec x0(vdpode::n);
    // x0[0] = rocs::interval(-3,3);
    // x0[1] = rocs::interval(-3,3);
    
    // std::cout << "Compute valid reachable set for " << x0 << ":\n";
    // std::cout << "Computing bound for Taylor terms: " << '\n'
    // 	      << vdp.eval_taylorterm_bound(x0);
    // rocs::ivec xt(vdpode::n);
    // T.reset();
    // std::cout << vdp.compute_reachset_valid(xt, x0) << ": ";
    // std::cout << "Time taken: " << T.elapsed() << '\n';
    // std::cout << xt << std::endl;


    file.close();

    return 0;
}
