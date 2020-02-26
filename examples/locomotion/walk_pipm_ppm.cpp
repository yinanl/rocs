/**
 *  walk_pipm_ppm.cpp
 *  Contact switch from PIPM to PPM.
 *
 *  Created by Yinan Li on Sept. 23, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <time.h>
#include <set>

#include <vector>
#include <iostream>

#include "src/abstraction.hpp"
#include "src/dsolver.h"
#include "src/matlabio.h"

#include "walkstep.hpp"


/***** KEYFRAME STATES *****/
const double Dx= 0.0002;
const double z0 = 1e-5;

const double xfoot1 = 0.0;
const double vxapex1 = 0.5;
const double w1 = 3.1321;

// const double xfoot2 = 0.5;
// const double vxapex2 = 0.6;
// const double w2 = 2.8592;

const double xarm1 = 0.6;
const double vxapex3 = 1.7;
const double w3 = 3.7436;

// const double xarm2 = 0.68;
// const double vxapex4 = 1.7;
// const double w4 = 3.7436;
/***** KEYFRAME STATES *****/


/***** ROBUST MARGINS *****/
const double ds1[] = {-0.002, 0.002};
const double dz1[] = {-0.05, 0.05};
const double ds2[] = {-0.06, 0.06};
const double dz2[] = {-0.005, 0.005};

const double zeta_r1[] = {0.0, 10}; // zeta boundaries for robust tube intersection
const double zeta_r2[] = {-1000000, -0.1}; 
/***** ROBUST MARGINS *****/




int main()
{
    /* Make sure the subfolder pipm2ppm is created. */
    const char *filename = "pipm2ppm/data_walk_pipm2ppm_D1_005.mat";
    /* Sampling time */
    double tau = 0.02; //<=0.05

    /* Dimensions */
    double n = 2;
    double m = 1;

    /* State space */
    double xlb[] = {-0.15, 0.2};
    double xub[] = {0.8, 1.8};
    
    /* Control values */
    double ulb[] = {2}; // the first variable denotes mode
    double uub[] = {4};
    double mu[] = {0.02}; // dw=0.02~0.04

    /* Robustness margins */
    double e1[] = {0,0};
    double e2[] = {0,0};

    /* Discretization parameter */
    double eta[] = {0.005, 0.005};

    
    /* Control system */
    WalkStep step1("semi1", tau, n, m, PIPM, xfoot1);
    step1.init_workspace(xlb, xub);
    step1.init_inputset(mu,ulb,uub);
    
    WalkStep step2("semi2", tau, n, m, PPM, xarm1);
    step2.init_workspace(xlb, xub);
    step2.init_inputset(mu,ulb,uub);
    
    
    /* Abstraction */
    rocs::abstraction<WalkStep> abst1(&step1);
    abst1.init_state(eta, xlb, xub);
    rocs::abstraction<WalkStep> abst2(&step2);
    abst2.init_state(eta, xlb, xub);


    /* Assign transitions */
    clock_t tb, te;
    tb = clock();
    abst1.assign_transitions(e1, e2);
    abst2.assign_transitions(e1, e2);
    te = clock();
    std::cout << "Time of computing abstraction: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';
    // rocs::matWriter tswtr("pipm2ppm/ts_pipm2ppm_D5_002.mat");
    // tswtr.open();
    // tswtr.write_transitions(abst2._ts);
    // tswtr.write_transitions(abst1._ts);
    // tswtr.close();


    /*** Solve the one-step walking problem by two semi-steps ***/
    RobustSet Rq2(xarm1, vxapex3, -w3*w3, z0, Dx, ds2, dz2);
    RobustSet Rq1(xfoot1, vxapex1, w1*w1, z0, Dx, ds1, dz1);
    InterSetPIPMPPM Rinter(&Rq1, &Rq2, zeta_r1, zeta_r2);
    
    tb = clock();
    
    /* Get the target set */
    std::function<bool(const rocs::ivec &)> in_target =
    	std::bind(&RobustSet::in_robust_set, &Rq2, std::placeholders::_1);
    std::vector<size_t> target = abst2.get_discrete_states(in_target);
    /* Solve the second semi-step */
    rocs::DSolver solver2(&(abst2._ts));
    solver2.reachability(target);
    /* Determine the guard set */
    std::vector<size_t> guard;
    rocs::ivec x(n);
    for (int i = 0; i < solver2._win.size(); ++i) {
    	if (solver2._win[i]) {
	    for (int j = 0; j < n; ++j) {
		x.setval(j, rocs::interval(abst1._x._data[i][j]-abst1._x._gw[j]/2,
					   abst1._x._data[i][j]+abst1._x._gw[j]/2));
	    }
	    
	    if (Rinter.in_interset_sigma(x) && Rinter.in_interset_zeta(x.mid())) {
		guard.push_back(i);
	    }
    	}
    }    
    /* Get the initial set */
    std::function<bool(const rocs::ivec &)> in_initial =
    	std::bind(&RobustSet::in_robust_set, &Rq1, std::placeholders::_1);
    std::vector<size_t> initial = abst1.get_discrete_states(in_initial);
    /* Solve the first semi-step */
    rocs::DSolver solver1(&(abst1._ts));
    solver1.reachability(guard);
    te = clock();
    std::cout << "Time of control synthesis: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';

    
    /*** Save data to file ***/
    rocs::matWriter wtr(filename);
    wtr.open();
    wtr.write_uniform_grids(abst1._x, "xgrid");
    wtr.write_uniform_grids(abst1._ptrsys->_ugrid, "ugrid");
    wtr.write_real_array(target, "goalset");
    wtr.write_real_array(initial, "initset");
    wtr.write_real_array(guard, "guardset");
    wtr.write_discrete_controller(solver2, "leastctlr2", "optctlr2");
    wtr.write_discrete_controller(solver1, "leastctlr1", "optctlr1");
    wtr.close();
    
    return 0;
}
