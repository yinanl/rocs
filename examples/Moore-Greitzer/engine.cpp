/**
 *  engine.cpp
 *
 *  Reach and stay control of the ode model of Moore-Greitzer engine.
 *
 *  Created by Yinan Li on June 9, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"

const double a = 1./3.5;
const double B = 2.0;
const double H = 0.18;
const double W = 0.25;
const double lc = 8.0;
const double cx = 1.0/lc;
const double cy = 1.0/(4*lc*B*B);
const double aH = a+H;
const double H2 = H/(2.0*W*W*W);
const double W2 = 3*W*W;

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

/** Target ball f(x)<=0
 * x^2+y^2<=r^2 <==> (x^2+y^2)-r^2<=0
 *
 **/
template<typename T>
T target_ball(const T &x) {
  const double r = 0.003;
  T y(1);
  y[0] = (x[0]-0.4519)*(x[0]-0.4519)+(x[1]-0.6513)*(x[1]-0.6513)-r*r;
  return y;
}


int main(int argc, char *argv[]) {

    /* set the state space */
    double xlb[3]{0.44, 0.6, 0.5};
    double xub[3]{0.54, 0.7, 0.8};

    /* set the control values */
    double Lmu = 0.01;
    double Lu = 0.05;
    double ulb[2]{-Lu, -Lmu};
    double uub[2]{Lu, Lmu};
    double mu[2]{Lu/5, 2*Lmu/5};

    /* set the sampling time and disturbance */
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 10;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);

    /* define the control system */
    double t = 0.1;
    rocs::CTCntlSys<mgode> engine("Moore-Greitzer", t, mgode::n, mgode::nu,
    				  delta, &controlparams);
    engine.init_workspace(xlb, xub);
    engine.init_inputset(mu, ulb, uub);
    engine.allocate_flows();

     /* set the specifications */
    /* case I */
    // std::stringstream convert(argv[1]);
    // double e;
    // convert >> e;
    // double e = 0.003;
    // double glb[]{0.5039-e, 0.6605-e, 0.5};
    // double gub[]{0.5039+e, 0.6605+e, 0.8};
    // double olb[]{0.520, 0.658, 0.5};
    // double oub[]{0.526, 0.// 664, 0.8};
    
    /* case II */
    double e = 0.003;
    double glb[]{0.4519-e, 0.6513-e, 0.5};
    double gub[]{0.4519+e, 0.6513+e, 0.8};
    double olb[]{0.497, 0.650, 0.5};
    double oub[]{0.503, 0.656, 0.8};

    /* solve the problem */
    rocs::CSolver solver(&engine, 0, rocs::RELMAX, 60);
    // solver.init(rocs::GOAL, glb, gub);
    double et[]{0.001, 0.001, 1};
    solver.init(rocs::GOAL, &target_ball<rocs::ivec>, et);
    solver.init_goal_area();
    solver.init(rocs::AVOID, olb, oub);
    solver.init_avoid_area();
    solver.print_controller_info();
    // solver.print_controller();

    double ei[]{0.0002, 0.0002, 0.004};
    double er[]{0.005, 0.005, 0.005};
    double ermin[]{0.0002, 0.0002, 0.002};
    // solver.cobuchi(&engine, ei, er);
    solver.reachstay_control(&engine, ei, er, true, ermin);
    solver.print_controller_info();

    /* save the problem data and the solution */
    rocs::matWriter wtr("data_caseIIReachstay.mat");
    wtr.open();
    wtr.write_problem_setting(engine, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    engine.release_flows();

    
    // double tr = 1;
    // rocs::CTCntlSys<mgode> engine_r("Moore-Greitzer", tr, mgode::n, mgode::nu,
    // 				    delta, &controlparams);
    // engine_r.init_workspace(xlb, xub);
    // engine_r.init_inputset(mu, ulb, uub);
    // engine_r.allocate_flows();

    // rocs::CSolver solver(&engine_r, 0, rocs::RELMAX, 60);
    // // solver.init(rocs::GOAL, glb, gub);
    // double et[]{0.001, 0.001, 1};
    // solver.init(rocs::GOAL, &target_ball<rocs::ivec>, et);
    // solver.init_goal_area();
    // solver.init(rocs::AVOID, olb, oub);
    // solver.init_avoid_area();
    // solver.print_controller_info();

    // double er[]{0.005, 0.005, 0.005};
    // double ermin[]{0.0002, 0.0002, 0.002};
    // solver.reachability_control(&engine_r, er, true, ermin);

    // rocs::matWriter wtr("data_caseIIReach_1s.mat");
    // wtr.open();
    // wtr.write_problem_setting(engine_r, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();

    // engine_r.release_flows();

    return 0;
}
