/**
 *  safesetAbst.cpp
 *
 *  Compute the safety (collision-free) set of two agents
 *  by using abstraction-based methods.
 *
 *  Created by Yinan Li on July 13, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>
#include <sys/stat.h>

#include "src/system.hpp"
#include "src/abstraction.hpp"
#include "src/DBAparser.h"
#include "src/dsolver.h"
#include "src/bsolver.hpp"

#include "src/hdf5io.h"


struct twoagent {
    static const int n = 3;  // system dimension
    static const int nu = 2;  // control dimension
    rocs::ivec d{rocs::interval(-0.8, 0.8),
		 rocs::interval(-0.8, 0.8)};

    /* template constructor
     * @param[out] dx
     * @param[in] x = [xr, yr, psir]
     * @param u = [v, w]
     * @param d = [v', w']
     */
    template<typename S>
    twoagent(S *dx, const S *x, rocs::Rn u) {
	dx[0] = -u[0] + d[0]*cos(x[2]) + u[1]*x[1];
	dx[1] = d[0]*sin(x[2]) - u[1]*x[0];
	dx[2] = d[1] - u[1];
    }
};


int main(int argc, char *argv[])
{
    std::string specfile = "Ga.txt";
    double eta[]{0.2, 0.2, 0.2}; /* partition precision */

    /* Input arguments:
     * carAbst [dbafile precision(e.g. 0.2 0.2 0.2)]
     */
    if (argc > 5) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }
    if(argc > 1) {
	specfile = std::string(argv[1]);
	if(argc >= 5) {
	    for(int i = 2; i < 5; ++i)
		eta[i-2] = std::atof(argv[i]);
	} else {
	    std::cout << "Input precision should be of 3-dim, e.g. [0.2, 0.2, 0.2,].\n";
	    std::exit(1);
	}
    }

    std::cout << "Partition precision: "
	      << eta[0] << ' ' << eta[1] << ' ' << eta[2]
	      << '\n';


    /**
     * Set the state and control space.
     */
    const double theta = 3.5;
    double xlb[] = {-3, -3, -theta};
    double xub[] = {3, 3, theta};
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /**
     * Define the two-agent system
     */
    double t = 0.3;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    rocs::CTCntlSys<twoagent> safety("collision-free", t,
				     twoagent::n, twoagent::nu,
				     delta, &controlparams);

    safety.init_workspace(xlb, xub);
    safety.init_inputset(mu, ulb, uub);
    safety.allocate_flows();

    rocs::abstraction< rocs::CTCntlSys<twoagent> > abst(&safety);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    /**
     * Assign 1 to the target invariant set and 0 to others.
     * Mark 0 for any box intersect or inside the cylinder: x^2+y^2<=rmin^2, any phi.
     * The invariant set is the region outside of the cylinder.
     */
    auto inv_set = [&abst, &eta](size_t i) {
		       const double rmin = 1.21;
		       std::vector<double> x(abst._x._dim);
		       abst._x.id_to_val(x, i);
		       double xl = x[0] - eta[0]/2.;
		       double xr = x[0] + eta[0]/2.;
		       double yl = x[1] - eta[1]/2.;
		       double yr = x[1] + eta[1]/2.;
		       double xsqr = (xr*xr) > (xl*xl) ? (xl*xl) : (xr*xr);
		       double ysqr = (yr*yr) > (yl*yl) ? (yl*yl) : (yr*yr);
		       if(xsqr + ysqr < rmin*rmin)
			   // The closest corner to the circle is inside the circle.
			   return 0;
		       // else if(xr>-rmin && xl<rmin && yr>-rmin && yl<rmin)
		       // 	   // The distance between the center of the box and the circle
		       // 	   // in both x, y direction should be less than (rmin + box radius)
		       // 	   return 0;
		       else
			   return 1;
		   };
    abst.assign_labels(inv_set);
    abst.assign_label_outofdomain(1); //out of domain is safe
    /* Compute abstraction */
    clock_t tb, te;
    std::string transfile = "abstca_0.8-1.2-0.3-0.2.h5";
    /* https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c */
    struct stat buffer;
    if(stat(transfile.c_str(), &buffer) == 0) {
	/* Read from a file */
	std::cout << "Transitions have been computed. Reading transitions...\n";
	rocs::h5FileHandler transRdr(transfile, H5F_ACC_RDONLY);
	tb = clock();
	transRdr.read_transitions(abst._ts);
	te = clock();
	float time = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of reading abstraction: " << time << '\n';
	std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
    } else {
	std::cout << "Transitions haven't been computed. Computing transitions...\n";
	/* Robustness margins */
	double e1[] = {0,0,0};
	double e2[] = {0,0,0};
	tb = clock();
	abst.assign_transitions(e1, e2);
	te = clock();
	float tabst = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of computing abstraction: " << tabst << '\n';
	std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
	/* Write transitions to file */
	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
	transWtr.write_transitions(abst._ts);
    }

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
     * Invariance control synthesis
     */
    rocs::DSolver solver(&abst._ts);
    tb = clock();
    solver.invariance(targetIDs);
    te = clock();
    float tsyn = (float)(te - tb)/CLOCKS_PER_SEC;
    std::cout << "Time of control synthesis: " << tsyn << '\n';
    std::cout << "The size of winning set: " << solver._nw << '\n';
    /**
     * Write the control synthesis result into .h5 file.
     */
    std::string datafile = "controller_safety_abst_0.8-1.2-0.3-0.2.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::CTCntlSys<twoagent> >(safety);
    ctlrWtr.write_2d_array<double>(targetPts, "G");
    ctlrWtr.write_array<double>(eta, 3, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(solver);


    // /**
    //  * Read DBA from dba*.txt file
    //  */
    // std::cout << "Reading the specification...\n";
    // rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    // std::vector<rocs::UintSmall> acc;
    // std::vector<std::vector<rocs::UintSmall>> arrayM;
    // if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc))
    // 	std::exit(1);

    
    // /**
    //  * Solve a Buchi game on the product of NTS and DBA.
    //  */
    // std::cout << "Start solving a Buchi game on the product of the abstraction and DBA...\n";
    // rocs::BSolver solver; // memories will be allocated for psolver
    // solver.construct_dba((int)nAP, (int)nNodes, (int)q0, acc, arrayM);
    // tb = clock();
    // solver.load_abstraction(abst);
    // solver.generate_product(abst);
    // solver.solve_buchigame_on_product();
    // te = clock();
    // float tsyn = (float)(te - tb)/CLOCKS_PER_SEC;
    // std::cout << "Time of synthesizing controller: " << tsyn << '\n';
    // /**
    //  * Write the optimal controller to a txt file.
    //  */
    // std::string datafile = "controller_safety_abst_0.8-1.2-0.3-0.2.txt";
    // std::cout << "Writing the controller...\n";
    // solver.write_controller_to_txt(const_cast<char*>(datafile.c_str()));
    // std::cout << "Controller writing is done.\n";


    safety.release_flows();
    return 0;
}
