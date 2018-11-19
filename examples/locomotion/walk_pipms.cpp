/**
 *  walk_pipms.cpp
 *  one walking step in PIPM mode between two keyframe states.
 *
 *  Created by Yinan Li on Sept 22, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <time.h>
#include <set>
#include <vector>
#include <functional>
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

const double xfoot2 = 0.5;
const double vxapex2 = 0.6;
const double w2 = 2.8592;
/***** KEYFRAME STATES *****/


/***** ROBUST MARGINS *****/
const double ds1[] = {-0.002, 0.002};
const double dz1[] = {-0.05, 0.05};
const double ds2[] = {-0.006, 0.006};
const double dz2[] = {-0.05, 0.05};
/***** ROBUST MARGINS *****/


int main()
{
    const char *filename = "pipms/data_walk_pipms_D2_004.mat";
    
    /* Sampling time */
    double tau = 0.02; //<=0.05

    /* Dimensions */
    double n = 2;
    double m = 1;

    /* State space */
    double xlb[] = {-0.15, 0.2};
    double xub[] = {0.7, 1.1};
    
    /* Control values */
    double ulb[] = {2}; // the first variable denotes mode
    double uub[] = {4};
    double mu[] = {0.02}; // dw=0.02~0.04

    /* Robustness margins */
    double e1[] = {0,0};
    double e2[] = {0,0};

    
    /* Control system */
    WalkStep step1("semi1", tau, n, m, PIPM, xfoot1);
    step1.init_workspace(xlb, xub);
    step1.init_inputset(mu,ulb,uub);
    WalkStep step2("semi2", tau, n, m, PIPM, xfoot2);
    step2.init_workspace(xlb, xub);
    step2.init_inputset(mu,ulb,uub);
    
    
    /* Abstraction */
    // double eta[] = {0.005, 0.005};
    double eta[] = {0.004, 0.004};
    rocs::abstraction<WalkStep> abst1(&step1);
    abst1.init_state(n, eta, xlb, xub);
    rocs::abstraction<WalkStep> abst2(&step2);
    abst2.init_state(n, eta, xlb, xub);


    /* Assign transitions */
    clock_t tb, te;
    tb = clock();
    abst1.assign_transitions(e1, e2);
    abst2.assign_transitions(e1, e2);
    te = clock();
    std::cout << "Time of computing abstraction: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';


    /*** Solve the one-step walking problem by two semi-steps ***/
    /* Define two robust sets */
    RobustSet Rq2(xfoot2, vxapex2, w2*w2, z0, Dx, ds2, dz2);
    RobustSet Rq1(xfoot1, vxapex1, w1*w1, z0, Dx, ds1, dz1);
    //InterSetPIPMs Rinter(&Rq1, &Rq2);
    
    tb = clock();
    /* Solve the second semi-step */
    std::function<bool(const rocs::ivec &)> in_target =
    	std::bind(&RobustSet::in_robust_set, &Rq2, std::placeholders::_1);
    std::vector<size_t> target = abst2.get_discrete_states(in_target);
    rocs::DSolver solver2(&(abst2._ts));
    solver2.reachability(target);
    /* Determine the guard set */
    std::vector<size_t> guard;
    rocs::ivec x(n);
    for (int i = 0; i < solver2._winset.size(); ++i) {
    	if (solver2._winset[i]) {
	    for (int j = 0; j < n; ++j) {
		x.setval(j, rocs::interval(abst1._x._data[i][j]-abst1._x._gw[j]/2,
					   abst1._x._data[i][j]+abst1._x._gw[j]/2));
	    }
	    // if (Rinter.in_interset(x)) {
	    // 	guard.push_back(i);
	    // }
	    if (Rq1.in_robust_tube(x)) {
		guard.push_back(i);
	    }
    	}
    }
    /* Solve the first semi-step */
    std::function<bool(const rocs::ivec &)> in_initial =
    	std::bind(&RobustSet::in_robust_set, &Rq1, std::placeholders::_1);
    std::vector<size_t> initial = abst1.get_discrete_states(in_initial);
    rocs::DSolver solver1(&(abst1._ts));
    solver1.reachability(guard); // solve the first semi step
    te = clock();
    std::cout << "Time of control synthesis: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';

    
    /* Save data to file */
    rocs::matWriter wtr(filename);
    wtr.open();
    wtr.write_uniform_grids(abst1._x, "xgrid");
    wtr.write_uniform_grids(abst1._ptrsys->_ugrid, "ugrid");
    wtr.write_real_array(target, "goalset");
    wtr.write_real_array(guard, "guardset");
    wtr.write_real_array(initial, "initset");
    //wtr.write_transitions(abst2._ts);
    wtr.write_discrete_controller(solver2, "leastctlr2", "optctlr2");
    //wtr.write_transitions(abst1._ts);
    wtr.write_discrete_controller(solver1, "leastctlr1", "optctlr1");
    wtr.close();
    
    return 0;
}
