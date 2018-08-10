/**
 *  roa.cpp
 *
 *  Compute the region of attraction (ROA) for reversed Van de Pol
 *
 *  Created by Yinan Li on April 27, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */



#include <iostream>

// #include "timer.hpp"

#include "src/system.hpp"
#include "src/csolver.h"
#include "src/matlabio.h"


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

/** constraint function f(x)<=0
 * xT P x -c <= 0, P = [1.5 -0.5; -0.5 1], c = 1.43
 *
 **/
template<typename T>
T roa_core(const T &x) {
  const double c = 1.43;
  T y(1);
  y[0] = x[0]*(1.5*x[0]-x[1]) + x[1]*x[1] - c;
  return y;
}


int main()
{
    /* set the state space */
    double xlb[] = {-4, -4};
    double xub[] = {4, 4};
    
    /* set the sampling time and disturbance */
    double tau = 0.05;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    
    /* define the control system */
    rocs::CTSys<vdpode> vdproa("VanDePol", tau, vdpode::n, delta, &controlparams);
    vdproa.init_workspace(xlb, xub);
    vdproa.allocate_flows();

    /* compute ROA */
    /* eps = 0.03 */
    rocs::CSolver solver1(&vdproa);
    solver1.init(rocs::GOAL, &roa_core<rocs::ivec>, 0.01);
    solver1.cobuchi(&vdproa, 0.03, rocs::ABSMAX, 0.03, rocs::ABSMAX);
    solver1.print_controller_info();
    /* eps = 0.01 */
    rocs::CSolver solver2(&vdproa);
    solver2.init(rocs::GOAL, &roa_core<rocs::ivec>, 0.01);
    solver2.cobuchi(&vdproa, 0.01, rocs::ABSMAX, 0.01, rocs::ABSMAX);
    solver2.print_controller_info();
    /* eps = 0.005 */
    rocs::CSolver solver3(&vdproa);
    solver3.init(rocs::GOAL, &roa_core<rocs::ivec>, 0.01);
    solver3.cobuchi(&vdproa, 0.005, rocs::ABSMAX, 0.005, rocs::ABSMAX);
    solver3.print_controller_info();
    
    
    /* save the problem data and the solutions */
    rocs::matWriter wtr("data_roavdp.mat");
    wtr.open();
    wtr.write_problem_setting(vdproa, solver1);
    // wtr.write_controller(solver1);
    double eta[] = {0.001, 0.001};
    wtr.write_winset_boundary(eta, solver1, "bd01");
    wtr.write_winset_boundary(eta, solver2, "bd02");
    wtr.write_winset_boundary(eta, solver3, "bd03");
    wtr.close();


    vdproa.release_flows();
    return 0;
}
