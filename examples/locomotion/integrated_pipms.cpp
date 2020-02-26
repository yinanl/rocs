/**
 *  integrated_pipms.cpp
 *
 *  One walking step from PIPM to PIPM: 
 *  25 different cells for each robutst set.
 *
 *  Created by Yinan Li on Oct. 14, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <time.h>
#include <set>

#include <vector>
#include <iostream>

#include <stdlib.h>
#include <string>

#include "src/abstraction.hpp"
#include "src/dsolver.h"
#include "src/matlabio.h"

#include "walkstep.hpp"


/***** KEYFRAME STATES *****/
const double Dx= 0.0002;
const double z0 = 1e-5;

/* the initial keyframe state */
const double xfoot1 = 0.0;
const double v1 = 0.55;
const double w1 = 2.6592;

/* the final keyframe state */
const double xfoot2 = 0.5;
const double v2 = 0.55;
const double w2 = 2.6592;
/***** KEYFRAME STATES *****/


/***** ROBUST MARGINS *****/
const int N1 = 5;
const int N2 = 5;

const double ds1 = 0.002;
const double dz1 = 0.002;
const double ds2 = 0.002;
const double dz2 = 0.002;
/***** ROBUST MARGINS *****/


void reachability_control(rocs::abstraction<WalkStep> &abst1, const double s1[], const double z1[],
			  rocs::abstraction<WalkStep> &abst2, const double s2[], const double z2[],
			  const char *filename) {

    clock_t tb, te;

    /*** Solve the one-step walking problem by two semi-steps ***/
    RobustSet Rq2(xfoot2, v2, w2*w2, z0, Dx, s2, z2);
    RobustSet Rq1(xfoot1, v1, w1*w1, z0, Dx, s1, z1);
    
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
    rocs::ivec x(abst1._x._dim);
    for (int i = 0; i < solver2._win.size(); ++i) {
    	if (solver2._win[i]) {
    	    for (int j = 0; j < abst1._x._dim; ++j) {
    		x.setval(j, rocs::interval(abst1._x._data[i][j]-abst1._x._gw[j]/2,
    					   abst1._x._data[i][j]+abst1._x._gw[j]/2));
    	    }
    	    if (Rq1.in_robust_tube(x)) {
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
    solver1.reachability(guard); // solve the first semi step
    te = clock();
    std::cout << "Time of control synthesis: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';

    
    /* Save data to file */
    rocs::matWriter wtr(filename);
    wtr.open();
    wtr.write_real_array(target, "goalset");
    wtr.write_real_array(initial, "initset");
    wtr.write_real_array(guard, "guardset");
    wtr.write_discrete_controller(solver2, "leastctlr2", "optctlr2");
    wtr.write_discrete_controller(solver1, "leastctlr1", "optctlr1");
    wtr.close();
}



int main()
{
    /* File name components */
    char filename[100];
    std::string strq1, strq2;
    const std::string prefix = "integrated/pipms_2/data_integrate_pipms_q";
    const std::string suffix = ".mat";
    const std::string connect = "_p";
    
    double loc1 = xfoot1;
    double loc2 = xfoot2;
    
    /* Sampling time */
    double tau = 0.02; //<=0.05

    /* Dimensions */
    double n = 2;
    double m = 1;

    /* State space */
    double xlb[] = {-0.2, 0.2};
    double xub[] = {0.7, 1.1};
    
    /* Control values */
    double ulb[] = {2}; // the first variable denotes mode
    double uub[] = {4};
    double mu[] = {0.02}; // dw=0.02~0.04

    /* Robustness margins */
    double e1[] = {0,0};
    double e2[] = {0,0};

    /* Control system */
    WalkStep step1("semi1", tau, n, m, PIPM, loc1);
    step1.init_workspace(xlb, xub);
    step1.init_inputset(mu,ulb,uub);
    WalkStep step2("semi2", tau, n, m, PIPM, loc2);
    step2.init_workspace(xlb, xub);
    step2.init_inputset(mu,ulb,uub);
    
    
    /* Abstraction */
    double eta[] = {0.003, 0.003};
    rocs::abstraction<WalkStep> abst1(&step1);
    abst1.init_state(eta, xlb, xub);
    rocs::abstraction<WalkStep> abst2(&step2);
    abst2.init_state(eta, xlb, xub);
    /* Save grids */
    rocs::matWriter wtr("integrated/pipms_2/data_grids_pipms.mat");
    wtr.open();
    wtr.write_uniform_grids(abst1._x, "xgrid");
    wtr.write_uniform_grids(abst1._ptrsys->_ugrid, "ugrid1");
    wtr.write_uniform_grids(abst2._ptrsys->_ugrid, "ugrid2");
    wtr.close();


    /* Assign transitions by computation */
    clock_t tb, te;
    tb = clock();
    abst1.assign_transitions(e1, e2);
    abst2.assign_transitions(e1, e2);
    te = clock();
    std::cout << "Time of computing abstraction: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';
    
    // /* Save transitions (have some problem in writing becuase of the huge size) */
    // wtr.openfile("integrated/data_trans_pipm2ppm_ppm.mat");
    // wtr.write_transitions(abst2._ts, "trans_post", "post", "postptr", "trans_pre", "pre", "preptr");
    // wtr.close();
    // wtr.openfile("integrated/data_trans_pipm2ppm_pipm.mat");
    // wtr.write_transitions(abst1._ts, "trans_post", "post", "postptr", "trans_pre", "pre", "preptr");
    // wtr.close();

    /* Read transitions from file */
    // rocs::matReader rdr("integrated/data_trans_pipm2ppm_pipm.mat");
    // rdr.open();
    // rdr.read_transitions(abst1._ts);
    // rdr.close();
    // rdr.openfile("integrated/data_trans_pipm2ppm_ppm.mat");
    // rdr.open();
    // rdr.read_transitions(abst2._ts);
    // rdr.close();


    /*** Iterate over 9x9 keyframe states ***/
    double rbs1[2], rbz1[2], rbs2[2], rbz2[2];
    /* Loop initial keyframe states q */
    for (int i1 = 0; i1 < N1; ++ i1) {
	for (int j1 = 0; j1 < N2; ++ j1) {
	    rbs1[0] = (-N1 + 2*i1) * ds1;
	    rbs1[1] = (-N1+2 + 2*i1) * ds1;
	    rbz1[0] = (-N2 + 2*j1) * dz1;
	    rbz1[1] = (-N2+2 + 2*j1) * dz1;

	    /* Compute the index of the initial state */
	    strq1 = std::to_string(i1*N1 + j1+1);
	    // std::cout << "The initial state " << i1*N1 + j1+1 << '\n';
	    
	    /* Loop final keyframe states q' */
	    for (int i2 = 0; i2 < N1; ++ i2) {
		for (int j2 = 0; j2 < N2; ++ j2) {
		    rbs2[0] = (-N1 + 2*i2) * ds2;
		    rbs2[1] = (-N1+2 + 2*i2) * ds2;
		    rbz2[0] = (-N2 + 2*j2) * dz2;
		    rbz2[1] = (-N2+2 + 2*j2) * dz2;

		    /* Compute the index of the final state */
		    strq2 = std::to_string(i2*N1 + j2+1);
		    // std::cout << "The final state " << i2*N2 + j2+1 << '\n';

		    /* Compose the filename */
		    std::string sf = prefix + strq1 + connect + strq2 + suffix;
		    strcpy(filename, sf.c_str());
		    std::cout << "Writing to " << filename << '\n';

		    std::cout << "rbs1=[" << rbs1[0] << ", " << rbs1[1] << "]"
		    	      << ", rbz1=[" << rbz1[0] << ", " << rbz1[1] << "]" << '\n';
		    std::cout << "rbs2=[" << rbs2[0] << ", " << rbs2[1] << "]"
		    	      << ", rbz2=[" << rbz2[0] << ", " << rbz2[1] << "]" << '\n';
		    /* Reachability control and save results */
		    reachability_control(abst1, rbs1, rbz1, abst2, rbs2, rbz2, filename);
			
		    // if (i1==1 & j1==1 & i2==1 & j2==1) {
		    // 	/* Compute the index of the final state */
		    // 	strq2 = std::to_string(i2*N1 + j2+1);
		    	
		    // 	/* Compose the filename */
		    // 	std::string sf = prefix + strq1 + connect + strq2 + suffix;
		    // 	strcpy(filename, sf.c_str());
		    	
		    // 	std::cout << "rbs1=[" << rbs1[0] << ", " << rbs1[1] << "]"
		    // 		  << ", rbz1=[" << rbz1[0] << ", " << rbz1[1] << "]" << '\n';
		    // 	/* Reachability control and save results */
		    // 	reachability_control(abst1, rbs1, rbz1, abst2, rbs2, rbz2, filename);
		    // }
		}
	    }
	}
    }
    
    
    return 0;
}
