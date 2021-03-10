/**
 *  staticplan.cpp
 *
 *  Generate a static plan for a planar vehicle to satisfy a DBA
 *  by using abstraction-based control synthesis.
 *
 *  Created by Yinan Li on Feb. 8, 2021.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>

#include "src/DBAparser.h"
#include "src/abstraction.hpp"

#include "src/bsolver.hpp"
#include "src/patcher.h"

#include "src/hdf5io.h"

#include "car.hpp"


int main(int argc, char *argv[])
{
    /**
     * Default values
     **/
    std::string specfile;
    double eta[]{0.2, 0.2, 0.2}; /* partition precision */

    /* Input arguments:
     * carAbst dbafile precision(e.g. 0.2 0.2 0.2)
     */
    if (argc < 2 || argc > 5) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }
    specfile = std::string(argv[1]);
    if (argc > 2 && argc < 4) {
	std::cout << "Input precision should be of 3-dim, e.g. 0.2 0.2 0.2.\n";
	std::exit(1);
    }
    if (argc == 5) {
	for(int i = 2; i < 5; ++i)
	    eta[i-2] = std::atof(argv[i]);
    }
    std::cout << "Partition precision: " << eta[0] << ' '
	      << eta[1] << ' ' << eta[2] << '\n';

    clock_t tb, te;
    /* set the state space */
    const double theta = 3.5;
    double xlb[] = {0, 0, -theta};
    double xub[] = {10, 10, theta};
    /* set the control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};
    /* define the control system */
    rocs::DTCntlSys<carde> car("DBA", tau, carde::n, carde::m);
    car.init_workspace(xlb, xub);
    car.init_inputset(mu, ulb, uub);


    /**
     * Construct and save the abstraction
     */
    // const double eta[] = {0.2, 0.2, 0.2}; /* set precision */
    rocs::abstraction<rocs::DTCntlSys<carde>> abst(&car);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    /* Assign the label of avoid area to -1 */
    rocs::UintSmall nAvoid = 6;
    double obs[6][4] = {
	{0.0, 1.0, 4.6, 6.4},
	{2.0, 3.0, 0.0, 3.6},
	{3.5, 4.5, 8.5, 10.0}, //{3.5, 4.5, 7.5, 10.0},
	{5.5, 6.5, 0.0, 1.0},
	{5.5, 6.5, 3.4, 5.5}, //{5.5, 6.5, 3.4, 6.5},
	{6.5, 10.0, 4.5, 5.5} //{6.5, 10.0, 5.5, 6.5}
    };
    auto label_avoid = [&obs, &nAvoid, &abst, &eta](size_t i) {
    		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, i);
    		     double c1= eta[0]/2.0; //+1e-10;
    		     double c2= eta[1]/2.0; //+1e-10;
    		     for(size_t i = 0; i < nAvoid; ++i) {
    			 if ((obs[i][0]-c1) <= x[0] && x[0] <= (obs[i][1]+c1) &&
    			     (obs[i][2]-c2) <= x[1] && x[1] <= (obs[i][3]+c2))
    			     return -1;
    		     }
    		     return 0;
    		 };
    abst.assign_labels(label_avoid);
    abst.assign_label_outofdomain(-1); // out of domain is banned
    std::vector<size_t> obstacles;
    for (size_t i = 0; i < abst._x._nv; ++i) {
    	if (abst._labels[i] < 0)
    	    obstacles.push_back(i);
    }
    /* Compute abstraction */
    std::string suffix;
    for(auto &item : eta) {
	std::stringstream ss;
    	ss << std::setprecision(1);
	ss << item;
	suffix += '-';
	suffix += ss.str();
    }
    std::string transfile = "abstfull" + suffix + ".h5";
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
	std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
    } else {
	std::cout << "Transitions haven't been computed. Computing transitions...\n";
	/* Robustness margins */
	double e1[] = {0,0,0};
	double e2[] = {0,0,0};
	tb = clock();
	abst.assign_transitions(e1, e2);
	// abst.assign_transitions();
	te = clock();
	tabst = (float)(te - tb)/CLOCKS_PER_SEC;
	std::cout << "Time of computing abstraction: " << tabst << '\n';
	std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
	/* Write abstraction to file */
	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
	transWtr.write_transitions(abst._ts);
	transWtr.write_array<size_t>(obstacles, "obs");
	transWtr.write_array<double>(eta, carde::n, "eta");
	transWtr.write_2d_array<double>(abst._x._data, "xgrid");
	transWtr.write_problem_setting< rocs::DTCntlSys<carde> >(car);
    }


    /**
     * Read DBA from dba*.txt file
     */
    std::cout << "Reading the specification...\n";
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc))
	std::exit(1);
    std::vector<std::string> tokens;
    boost::split(tokens, specfile, boost::is_any_of("."));


    /*
     * Assign labels to states: has to be consistent with the dba file.
     */
    rocs::UintSmall nGoal = 3;
    double e = tau/2.0;
    double goal[3][4] = {{0.5+e, 2.0-e, 7.5+e, 9.5-e},
			 {8.0, 0.8-e, 8.0, 0.8-e}, // {5.5, 0.5-e, 8.5, 0.5-e},
			 {7.5+e, 9.5-e, 0.8+e, 3.0-e}}; //{7.5+e, 9.5-e, 0.8+e, 4.0-e}
    auto label_target = [&goal, &nGoal, &abst, &eta](size_t i) {
		      std::vector<double> x(abst._x._dim);
		      abst._x.id_to_val(x, i);
		      // double c1= eta[0]/2.0; //+1e-10;
		      // double c2= eta[1]/2.0; //+1e-10;
		      double xl = x[0] - eta[0]/2.;
		      double xr = x[0] + eta[0]/2.;
		      double yl = x[1] - eta[1]/2.;
		      double yr = x[1] + eta[1]/2.;
		      boost::dynamic_bitset<> label(nGoal, false); // n is the number of goals
		      for(rocs::UintSmall i = 0; i < nGoal; ++i) {
			  if(i != 1) {
			  label[2-i] = (goal[i][0] <= xl && xr <= goal[i][1] &&
					goal[i][2] <= yl && yr <= goal[i][3])
			      ? true: false;
			  } else {
			      double rxr = xr - goal[i][0];
			      double rxl = xl - goal[i][0];
			      double ryr = yr - goal[i][2];
			      double ryl = yl - goal[i][2];
			      double xsqr = (rxr*rxr) < (rxl*rxl) ? (rxl*rxl) : (rxr*rxr);
			      double ysqr = (ryr*ryr) < (ryl*ryl) ? (ryl*ryl) : (ryr*ryr);
			      label[2-i] = (xsqr+ysqr)<goal[i][1]*goal[i][3] ? true : false;
			  }
		      }
		      return label.to_ulong();
		  };
    abst.assign_labels(label_target);
    
    std::cout << "Save labels to file.\n";
    std::string labelfile = "labels_" + tokens[0] + "_" + transfile;
    rocs::h5FileHandler labelWtr(labelfile, H5F_ACC_TRUNC);
    labelWtr.write_array<int>(abst._labels, "labels");
    std::cout << "Specification assignment is done.\n";


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
    std::string datafile = "controller_" + tokens[0] + "_";
    for(int i = 0; i < 3; ++i) {
    	std::stringstream ss;
    	ss << std::setprecision(1);
    	ss << eta[i];
    	datafile += ss.str();
    	if (i < 2)
    	    datafile += "-";
    }
    datafile += ".h5";
    std::cout << "Writing the controller...\n";
    // solver.write_controller_to_txt(const_cast<char*>(datafile.c_str()));
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::DTCntlSys<carde> >(car);
    // ctlrWtr.write_2d_array<double>(targetPts, "G");
    ctlrWtr.write_array<double>(eta, carde::n, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(&(solver._sol));
    std::cout << "Controller writing is done.\n";

    std::cout << "Total time of used (abstraction+synthesis): " << tabst+tsyn << '\n';


    /**
     * Create a patcher and save the winning graph
     */
    std::string winfile = "gwin.h5";
    rocs::Patcher local;
    std::cout << "Extracting the winning graph...\n";
    tb = clock();
    local.initialize_winning_graph(solver._sol);
    te = clock();
    std::cout << "Time of extracting the winning graph: " << (float)(te - tb)/CLOCKS_PER_SEC << '\n';
    
    /* Save the winning graph to a file */
    std::cout << "Writing the winning graph to file...\n";
    rocs::h5FileHandler graphWtr(winfile, H5F_ACC_TRUNC);
    tb = clock();
    if(graphWtr.write_winning_graph(local)) {
    	std::cout << "Error in saving the winning graph to file.\n";
    	return 1;
    }
    te = clock();
    std::cout << "Time of writing the winning graph: " << (float)(te - tb)/CLOCKS_PER_SEC << '\n';

    return 0;
}
