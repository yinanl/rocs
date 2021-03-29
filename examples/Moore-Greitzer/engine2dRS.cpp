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
// #include "src/matlabio.h"
#include "src/hdf5io.h"


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
struct mgode2 {
    static const int n = 2;  // system dimension
    static const int nu = 2;  // control dimension

    /* template constructor
     * @param[out] dx
     * @param[in] x
     * @param u
     */
    template<typename S>
    mgode2(S *dx, const S *x, rocs::Rn u) {
	dx[0] = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]) + u[0];
	dx[1] = cy * (x[0] - u[1]*sqrt(x[1]));
    }
};


int main(int argc, char *argv[]) {

    /* set the state space */
    double xlb[]{0.44, 0.6};
    double xub[]{0.54, 0.7};

    /* set the control values */
    double Lmu = 0.01;
    double Lu = 0.05;
    double ulb[]{-Lu, 0.5};
    double uub[]{Lu, 0.8};
    double mu[]{Lu/5, 0.01};

    /* set the sampling time and disturbance */
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 10;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);

    /* define the control system */
    double t= 0.1;
    rocs::CTCntlSys<mgode2> engine("Moore-Greitzer", t, mgode2::n, mgode2::nu,
				  delta, &controlparams);
    engine.init_workspace(xlb, xub);
    engine.init_inputset(mu, ulb, uub);
    engine.allocate_flows();

    /* set the specifications */
    // std::stringstream convert(argv[1]);
    // double e;
    // convert >> e;

    /* case I */
    double e = 0.003;
    double glb[]{0.5039-e, 0.6605-e};
    double gub[]{0.5039+e, 0.6605+e};
    double olb[]{0.520, 0.658};
    double oub[]{0.526, 0.664};

    // /* case II */
    // double e = 0.003;
    // double glb[]{0.4519-e, 0.6513-e};
    // double gub[]{0.4519+e, 0.6513+e};
    // double olb[]{0.497, 0.650};
    // double oub[]{0.503, 0.656};

    /* solve the problem */
    rocs::CSolver solver(&engine, 0, rocs::RELMAX, 70);
    solver.init(rocs::AVOID, olb, oub);
    solver.init_avoid_area();
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    solver.print_controller_info();
    solver.print_controller();

    double ei[]{0.00018, 0.00018};
    double er[]{0.005, 0.005};
    double ermin[]{0.00018, 0.00018};
    // solver.cobuchi(&engine, ei, ermin);
    solver.reachstay_control(&engine, ei, ermin);
    // solver.reachstay_control(&engine, ei, er, true, ermin);
    // solver.reachability_control(&engine, ermin);
    solver.print_controller_info();

    /* save the problem data and the solution */
    // rocs::matWriter wtr("data_2d_caseIReachstay.mat");
    // wtr.open();
    // wtr.write_problem_setting(engine, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string datafile = "controller_2d_caseIReachstay.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::CTCntlSys<mgode2> >(engine);
    ctlrWtr.write_ivec_array(solver._goal, "G");
    ctlrWtr.write_ivec_array(solver._obs, "xobs");
    ctlrWtr.write_sptree_controller(solver);


    engine.release_flows();

    return 0;
}
