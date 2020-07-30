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

double a = 1./3.5;
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


int main(int argc, char *argv[]) {

    /* set the state space */
    double xlb[3]{0.45, 0.6, 0.5};
    double xub[3]{0.55, 0.7, 0.8};

    /* set the control values */
    double Lmu = 0.01;
    double Lu = 0.05;
    double ulb[2]{-Lu, -Lmu};
    double uub[2]{Lu, Lmu};
    double mu[2]{Lu/5, 2*Lmu/5};

    /* set the sampling time and disturbance */
    double t = 0.1;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 10;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);

    /* define the control system */
    rocs::CTCntlSys<mgode> engine("Moore-Greitzer", t, mgode::n, mgode::nu,
				  delta, &controlparams);
	double tr = 1;
	rocs::CTCntlSys<mgode> engine_r("Moore-Greitzer", tr, mgode::n, mgode::nu,
				  delta, &controlparams);
    engine.init_workspace(xlb, xub);
    engine.init_inputset(mu, ulb, uub);
    engine.allocate_flows();
	engine_r.init_workspace(xlb, xub);
    engine_r.init_inputset(mu, ulb, uub);
    engine_r.allocate_flows();

    /* set the specifications */
    /* case I */
    // std::stringstream convert(argv[1]);
    // double e;
    // convert >> e;
    double e = 0.003;
    double glb1[]{0.5039-e, 0.6605-e, 0.5};
    double gub1[]{0.5039+e, 0.6605+e, 0.8};
    double olb1[]{0.520, 0.658, 0.5};
    double oub1[]{0.526, 0.664, 0.8};
    // /* case II */
    // // double e = 0.003;
    // double glb2[]{0.4519-e, 0.6513-e, 0.5};
    // double gub2[]{0.4519+e, 0.6513+e, 0.8};
    // double olb2[]{0.497, 0.650, 0.5};
    // double oub2[]{0.503, 0.656, 0.8};

    /* solve the problem */
    rocs::CSolver solver(&engine, rocs::RELMAX, 100);
    solver.init(rocs::GOAL, glb1, gub1);
    solver.init_goal_area();
    solver.init(rocs::AVOID, olb1, oub1);
    solver.init_avoid_area();
    solver.print_controller_info();
    // solver.print_controller();

    double ei[]{0.0002, 0.0002, 0.004};
    double er[]{0.005, 0.005, 0.005};
    double ermin[]{0.0002, 0.0002, 0.002};
    // solver.cobuchi(&engine, ei, er);
    // solver.reachstay_control(&engine, ei, er, true, ermin);
    solver.invariance_control(&engine, ei);
	solver.reachability_control(&engine_r, er, true, ermin);
    solver.print_controller_info();

    /* save the problem data and the solution */
    // rocs::matWriter wtr("data_caseICobuchi.mat");
    rocs::matWriter wtr("data_caseIReachstay.mat");
    // rocs::matWriter wtr("data_caseIReach.mat");
    wtr.open();
    wtr.write_problem_setting(engine, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    engine.release_flows();
    engine_r.release_flows();

    return 0;
}
