/**
 *  Safe landing control of DC9-30 with the abstraction-based engine.
 *
 *  Created by Yinan Li on Jan 11, 2021.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>
#include <boost/algorithm/string.hpp>
#include <cmath>
#include <sys/stat.h>

#include "src/system.hpp"
#include "src/abstraction.hpp"
#include "src/DBAparser.h"
#include "src/bsolver.hpp"
#include "src/hdf5io.h"


/* Parameters of the model */
double mg = 60000.0*9.81;
double mi = 1.0/60000; /* weight inverse: 1/m */



/* ODE of the longitudinal equation of motions for DC9-30 */
struct eomlong {
    static const int n = 3;  // state dimension
    static const int m = 2;  // control dimension
    
    /**
     * Constructors: 
     * @param[out] dx (dV, dgamma, dh)
     * @param[in] x (V,gamma,h).
     * @param[in] u (T, alpha) control input (thrust, angle of attack)
     */
    template<typename S>
    eomlong(S *dx, const S *x, rocs::Rn u) {
	double c = 1.25+4.2*u[1];
	dx[0] = mi*(u[0]*cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*sin(x[1]));
	dx[1] = (1.0/(60000*x[0]))*(u[0]*sin(u[1])+68.6*c*x[0]*x[0]-mg*cos(x[1]));
	dx[2] = x[0]*sin(x[1]);
    }
    
};

/** A target set in the form of f(x)<=0:
 * x(0)*sin(x(1)) >= -0.91
 *
 **/
template<typename T>
T target_area(const T &x) {
    T y(7);
    y[0] = -0.91-x[0]*sin(x[1]);
    y[1] = 63 - x[0];
    y[2] = x[0] - 75;
    y[3] = -3*M_PI/180 - x[1];
    y[4] = x[1];
    y[5] = -x[2];
    y[6] = x[2] - 2.5;
    return y;
}


int main(int argc, char *argv[])
{   
    /**
     * Define the control system 
     **/
    /* Set sampling time and disturbance */
    double tau = 0.25;
    double delta = 10;
    /* Set parameters for computation */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    rocs::CTCntlSys<eomlong> aircraft("landing", tau, eomlong::n,
				      eomlong::m, delta, &controlparams);
    /* set the state space */
    double xlb[] = {58, -3*M_PI/180, 0};
    double xub[] = {83, 0, 56};
    double ulb[] = {0,0};
    double uub[] = {32000, 8*M_PI/180};
    double mu[] = {32000, 9.0/8.0*M_PI/180};
    aircraft.init_workspace(xlb, xub);
    aircraft.init_inputset(mu, ulb, uub);
    // std::cout << "# of u= " << aircraft._ugrid._nv
    // 	      << ", dimension= " << aircraft._ugrid._dim << '\n';
    // for(int i = 0; i < aircraft._ugrid._nv; ++i) {
    // 	std::cout << aircraft._ugrid._data[i][0] << ','
    // 		  << aircraft._ugrid._data[i][1] << '\n';
    // }
    aircraft.allocate_flows();


    /**
     * Compute abstraction
     */
    /* Set the target set */
    const double eta[] = {25.0/362, 3*M_PI/180/66, 56.0/334};
    double glb[] = {63, -3*M_PI/180, 0};
    double gub[] = {75, 0, 2.5};
    rocs::abstraction< rocs::CTCntlSys<eomlong> > abst(&aircraft);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    auto target = [&abst, &eta, &glb, &gub](size_t i) {
		      std::vector<double> x(abst._x._dim);
		      abst._x.id_to_val(x, i);
		      double c[]{eta[0]/2.0, eta[1]/2.0, eta[2]/2.0}; //+1e-10;
		      if(x[0]-c[0] >= glb[0] && x[0]+c[0] <= gub[0] &&
			 x[1]-c[1] >= glb[1] && x[1]+c[1] <= gub[1] &&
			 x[2]-c[2] >= glb[2] && x[2]+c[2] <= gub[2])
			  return 1;
		      else
			  return 0;
		  };
    abst.assign_labels(target);
    abst.assign_label_outofdomain(0);
    /* Compute abstraction */
    clock_t tb, te;
    std::string transfile = "abstraction.h5";
    struct stat buffer;
    float tabst;
    if(stat(transfile.c_str(), &buffer) == 0) {
	/* Read from a file */
	std::cout << "Transitions have been computed. Reading transitions...\n";
	rocs::h5FileHandler transRdr(transfile, H5F_ACC_RDONLY);
	tb = clock();
	transRdr.read_transitions(abst._ts);
	te = clock();
	tabst = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of reading abstraction: " << tabst << '\n';
    } else {
	std::cout << "Transitions haven't been computed. Computing transitions...\n";
	/* Robustness margins */
	double e1[] = {0,0,0};
	double e2[] = {0,0,0};
	tb = clock();
	abst.assign_transitions(e1, e2);
	te = clock();
	tabst = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of computing abstraction: " << tabst << '\n';
	/* Write transitions to file */
	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
	transWtr.write_transitions(abst._ts);
    }
    std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
    
    std::vector<size_t> targetIDs;
    std::vector< std::vector<double> > targetPts;   //initial invariant set
    std::vector<double> x(abst._x._dim);
    for(size_t i = 0; i < abst._labels.size(); ++i) {
    	if(abst._labels[i] > 0) {
    	    targetIDs.push_back(i);
    	    abst._x.id_to_val(x, i);
    	    targetPts.push_back(x);
    	}
    }

    
    /**
     * Reachability control synthesis
     */
    rocs::DSolver solver(&abst._ts);
    tb = clock();
    solver.reachability(targetIDs);
    te = clock();
    float tsyn = (float)(te - tb)/CLOCKS_PER_SEC;
    std::cout << "Time of control synthesis: " << tsyn << '\n';
    std::cout << "The size of winning set: " << solver._nw << '\n';

    
    /**
     * Write the control synthesis result into .h5 file.
     */
    std::string datafile = "controller_abst_safelanding.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::CTCntlSys<eomlong> >(aircraft);
    ctlrWtr.write_2d_array<double>(targetPts, "G");
    ctlrWtr.write_array<double>(eta, 3, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(solver);


    aircraft.release_flows();
    return 0;
}
