/**
 *  dcdcAbstRS.cpp
 *
 *  Abstraction-based method for reach-stay control of DC-DC converter.
 *
 *  Created by Yinan Li on Jan 11, 2021.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>
#include <string>
#include <sys/stat.h>

#include "src/DBAparser.h"
#include "src/system.hpp"
#include "src/abstraction.hpp"
#include "src/bsolver.hpp"
#include "src/hdf5io.h"
// #include "src/matlabio.h"

#include "dcdc.hpp"


int main(int argc, char *argv[])
{
    /* Input arguments: 
     * dcdcAbstRS dbafile 
     */
    if (argc != 2) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }
    std::vector<std::string> tokens;
    boost::split(tokens, argv[1], boost::is_any_of("."));
    
    
    /* set the state space */
    double xlb[] = {0.649, 0.9898};
    double xub[] = {1.65, 1.19};
    /* set the target */
    double glb[] = {1.15, 1.09};
    double gub[] = {1.55, 1.17};

    /* define the control system */
    rocs::DTSwSys<dcde> dcdcInv("dcdc", tau, dcde::n, dcde::m);
    dcdcInv.init_workspace(xlb, xub);
    
    
    /**
     * Abstraction 
     */
    const double eta[]{0.001, 0.0002};
    rocs::abstraction< rocs::DTSwSys<dcde> > abst(&dcdcInv);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    auto target = [&abst, &eta, &glb, &gub](size_t& id) {
		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, id);
    		     double c[]{eta[0]/2.0, eta[1]/2.0}; //+1e-10;
		      if(x[0]-c[0] >= glb[0] && x[0]+c[0] <= gub[0] &&
			 x[1]-c[1] >= glb[1] && x[1]+c[1] <= gub[1])
			  return 1;
		      else
			  return 0;
		 };
    abst.assign_labels(target);
    abst.assign_label_outofdomain(0);
    /* Compute abstraction */
    clock_t tb, te;
    std::string transfile = "abst_" + tokens[0] + ".h5";
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
     * Read DBA from spec*.txt file
     */
    std::cout << "Reading the specification...\n";
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    std::string specfile = std::string(argv[1]);
    if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc)) 
	std::exit(1);

    
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
    std::string datafile = "controller_abst_" + tokens[0] + ".h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTSwSys<dcde> >(dcdcInv);
    ctlrWtr.write_2d_array<double>(targetPts, "G");
    ctlrWtr.write_array<double>(eta, 2, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(&(solver._sol));
    std::cout << "Controller writing is done.\n";

    std::cout << "Total time of used (abstraction+synthesis): " << tabst+tsyn << '\n';

    
    return 0;
}
