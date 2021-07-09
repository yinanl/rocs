/**
 *  Moore-Greitzer engine DBA control with the specification-guided engine.
 *  
 *  Created by Yinan Li on July 3, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>

#include <string>
#include <cmath>
#include <boost/numeric/odeint.hpp>

#include "src/DBAparser.h"
#include "src/abstraction.hpp"
#include "src/bsolver.hpp"
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


int main(int argc, char *argv[])
{
    /* Input arguments: 
     * engine2dAbstI dbafile 
     */
    if (argc != 2) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }

    
    clock_t tb, te;
    /* set the state space */
    double xlb[]{0.44, 0.6};
    double xub[]{0.54, 0.7};
    
    /* set the control values */
    // double Lmu = 0.01;
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


    /**
     * Abstraction 
     */
    rocs::abstraction< rocs::CTCntlSys<mgode2> > abst(&engine);
    const double eta[]{0.00018, 0.00018};
    // const double eta[]{0.002, 0.002};
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    double obs[][2]{{0.520, 0.526}, {0.658, 0.664}};
    auto avoid = [&obs, abst, eta](size_t& id) {
		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, id);
    		     double c1 = eta[0]/2.0; //+1e-10;
    		     double c2 = eta[1]/2.0; //+1e-10;
		     if ((obs[0][0]-c1) <= x[0] && x[0] <= (obs[0][1]+c1) &&
			 (obs[1][0]-c2) <= x[1] && x[1] <= (obs[1][1]+c2))
			 return -1;
		     return 0;
		 };
    abst.assign_labels(avoid);
    abst.assign_label_outofdomain(0);
    std::vector<size_t> obstacles;
    for (size_t i = 0; i < abst._x._nv; ++i) {
    	if (abst._labels[i] < 0)
    	    obstacles.push_back(i);
    }

    std::string transfile = "abstI_0.00018.h5";
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
	double e1[] = {0.0, 0.0};
	double e2[] = {0.0, 0.0};
	tb = clock();
	abst.assign_transitions(e1, e2);
	te = clock();
	tabst = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of computing abstraction: " << tabst << '\n';
	/* Write abstraction to file */
	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
	transWtr.write_transitions(abst._ts);
    }
    std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
    

    /**
     * Read DBA from spec*.txt file
     */
    std::cout << "Reading the specification...\n";
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    std::string specfile = std::string(argv[1]);
    if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc)) 
	std::exit(1);
    

    /* 
     * Assign labels to states: has to be consistent with the dba file.
     */
    double e = 0.003;
    double goal[][2]{{0.5039-e, 0.5039+e}, {0.6605-e, 0.6605+e}};
    
    auto label_target = [&goal, &abst, &eta](size_t i) {
			    std::vector<double> x(abst._x._dim);
			    abst._x.id_to_val(x, i);
			    double c1= eta[0]/2.0; //+1e-10;
			    double c2= eta[1]/2.0; //+1e-10;
			    
			    return (goal[0][0] <= (x[0]-c1) && (x[0]+c1) <= goal[0][1] && 
				    goal[1][0] <= (x[1]-c2) && (x[1]+c2) <= goal[1][1]) ?
				1: abst._labels[i];
			};
    abst.assign_labels(label_target);
    std::cout << "Specification assignment is done.\n";
    std::vector<size_t> targetIDs;
    std::vector<rocs::Rn> targetPts;   //initial invariant set
    rocs::Rn x(abst._x._dim);
    for(size_t i = 0; i < abst._labels.size(); ++i) {
    	if(abst._labels[i] > 0) {
    	    targetIDs.push_back(i);
    	    abst._x.id_to_val(x, i);
    	    targetPts.push_back(x);
    	}
    }
    

    /**
     * Solve a Buchi game on the product of NTS and DBA.
     */
    std::cout << "Start solving a Buchi game on the product of the abstraction and DBA...\n";
    rocs::BSolver solver; // memories will be allocated for psolver
    solver.construct_dba((int)nAP, (int)nNodes, (int)q0, acc, arrayM);
    tb = clock();
    solver.load_abstraction(abst);
    solver.generate_product(abst);
    solver.solve_buchigame_on_product();
    te = clock();
    float tsyn = (float)(te - tb)/CLOCKS_PER_SEC;
    std::cout << "Time of synthesizing controller: " << tsyn << '\n';

    /**
     * Display and save memoryless controllers.
     */
    std::cout << "Writing the controller...\n";
    // std::string datafile = "controller_abstI_0.00018.txt";
    // solver.write_controller_to_txt(const_cast<char*>(datafile.c_str()));
    std::string datafile = "controller_abstI_0.00018.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::CTCntlSys<mgode2> >(engine);
    ctlrWtr.write_2d_array<double>(targetPts, "G");
    ctlrWtr.write_array<double>(eta, mgode2::n, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(&(solver._sol));
    std::cout << "Controller writing is done.\n";

    std::cout << "Total time of used (abstraction+synthesis): " << tabst+tsyn << '\n';
    

    engine.release_flows();
    return 0;
}
