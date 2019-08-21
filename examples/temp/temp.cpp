/**
 *  temp.cpp
 *
 *  Temperature control example
 *
 *  Created by yinan li on August 14, 2017.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"


const double r1 = 0.002;
const double c = 16;
const double r2 = 0.1;

const double tau = 10;  // 50 secs
const double exp1 = std::exp(-r1*tau);
const double b1 = (1.0 - exp1)*c;


/*
 * 4-mode room temperature control:
 * discrete-time interval model
 */
struct tpc {
    static const int n = 2;
    static const int m = 4;
    
    template<typename S>
    tpc(S &y, const S &x, const int m) {
	switch(m) {
	case 1: /* mode 1: off */
	    y[0] = exp1*x[0] + b1;
	    y[1] = x[1];
	    break;
	case 2: /* mode 2: heating */
	    y[0] = exp1*x[0] + (1.0 - exp1)*(x[1]-r2/r1) + r2*tau;
	    y[1] = x[1]+r2*tau;
	    break;
	case 3: /* mode 3: cooling */
	    y[0] = exp1*x[0] + (1.0 - exp1)*(x[1]+r2/r1) - r2*tau;
	    y[1] = x[1]-r2*tau;
	    break;
	case 4: /* mode 4: on */
	    y[0] = exp1*x[0] + (1.0 - exp1)*x[1];
	    y[1] = x[1];
	    break;
	default:
	    break;
	}
    }
};


int main()
{
    /* set the state space */
    /* x[0]= room temperature, x[1]=heater temperature */
    double xlb[] = {10, 12};
    double xub[] = {28, 30};

    /* define the control system */
    rocs::DTSwSys<tpc> tpcReachstay("tpc", tau, tpc::n, tpc::m);
    tpcReachstay.init_workspace(xlb, xub);

    /* set the specification */
    double glb[] = {18, 20};
    double gub[] = {20, 22};

    double olb[] = {22, 26};
    double oub[] = {28, 30};

    /* solve the problem */
    rocs::CSolver solver(&tpcReachstay);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    solver.init(rocs::AVOID, olb, oub);
    solver.init_avoid_area();
    
    // solver.reach_stay(0.1, ABSMAX, 0.1, ABSMAX);
    // solver.invariance_control(0.5, ABSMAX);
    solver.cobuchi(&tpcReachstay, 0.15, rocs::ABSMAX, 0.15, rocs::ABSMAX);
    solver.print_controller_info();

    /* save the problem data and the solution */
    rocs::matWriter wtr("data_tpc_reachstay.mat");
    wtr.open();
    wtr.write_problem_setting(tpcReachstay, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    return 0;
}
