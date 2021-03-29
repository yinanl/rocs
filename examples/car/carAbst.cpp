/**
 *  carAbst.cpp
 *
 *  The vehicle example from scots using abstraction-based control synthesis.
 *
 *  Created by Yinan Li on Feb. 18, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <boost/algorithm/string.hpp>

#include "src/abstraction.hpp"
#include "src/DBAparser.h"
#include "src/bsolver.hpp"

#include "src/hdf5io.h"

#include "car.hpp"


int main(int argc, char *argv[])
{
    /**
     * Default values
     **/
    std::string specfile;
    double eta[]{0.2, 0.2, 0.2}; /* partition precision */
    bool preprocess = false;

    /* Input arguments:
     * carAbst dbafile precision(e.g. 0.2 0.2 0.2) -preprocess
     */
    if (argc < 2 || argc > 6) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }
    specfile = std::string(argv[1]);
    if (argc > 2 && argc < 5) {
	std::cout << "Input precision should be of 3-dim, e.g. 0.2 0.2 0.2.\n";
	std::exit(1);
    }
    if (argc >= 5) {

	for(int i = 2; i < 5; ++i)
	    eta[i-2] = std::atof(argv[i]);

    	if (argc > 5) {
    	    std::string param = std::string(argv[5]);
    	    if (param == "-preprocess")
    		preprocess = true;
    	    else {
    		std::cout << "Wrong argument for preprocessing.\n";
    		std::exit(1);
    	    }
    	}
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
    rocs::DTCntlSys<carde> car("reach goal", h, carde::n, carde::m);
    car.init_workspace(xlb, xub);
    car.init_inputset(mu, ulb, uub);


    /**
     * Construct and save abstraction
     */
    // const double eta[] = {0.2, 0.2, 0.2}; /* set precision */
    rocs::abstraction<rocs::DTCntlSys<carde>> abst(&car);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    /* Assign the label of avoid area to -1 */
    rocs::UintSmall nAvoid = 4;
    double obs[4][4] = {
	{1.6, 5.7, 4.0, 5.0},
	{3.0, 5.0, 5.0, 8.0},
	{4.3, 5.7, 1.8, 4.0},
	{5.7, 8.5, 1.8, 2.5}
    };
    auto label_avoid = [&obs, &nAvoid, &abst, &eta](size_t i) {
    		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, i);
    		     double c1= eta[0]/2.0+1e-10;
    		     double c2= eta[1]/2.0+1e-10;
    		     for(size_t i = 0; i < nAvoid; ++i) {
    			 if ((obs[i][0]-c1) <= x[0] && x[0] <= (obs[i][1]+c1) &&
    			     (obs[i][2]-c2) <= x[1] && x[1] <= (obs[i][3]+c2))
    			     return -1;
    		     }
    		     return 0;
    		 };
    abst.assign_labels(label_avoid);
    std::vector<size_t> obstacles;
    for (size_t i = 0; i < abst._x._nv; ++i) {
    	if (abst._labels[i] < 0)
    	    obstacles.push_back(i);
    }
    /* Compute abstraction */
    /* Robustness margins */
    double e1[] = {0,0,0};
    double e2[] = {0,0,0};
    tb = clock();
    abst.assign_transitions(e1, e2);
    // abst.assign_transitions();
    te = clock();
    float tabst = (float)(te - tb)/CLOCKS_PER_SEC;
    std::cout << "Time of computing abstraction: " << tabst << '\n';
    std::cout << "# of transitions: " << abst._ts._ntrans << '\n';


    /**
     * Read DBA from dba*.txt file
     */
    std::cout << "Reading the specification...\n";
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc))
	std::exit(1);


    /*
     * Assign labels to states: has to be consistent with the dba file.
     */
    double goal[3][4] = {{1.0, 2.0, 0.5, 2.0},
			 {0.5, 2.5, 7.5, 8.5},
			 {7.1, 9.1, 4.6, 6.4}};
    rocs::UintSmall nGoal = 3;
    auto label_target = [&goal, &nGoal, &abst, &eta](size_t i) {
		      std::vector<double> x(abst._x._dim);
		      abst._x.id_to_val(x, i);
		      double c1= eta[0]/2.0; //+1e-10;
		      double c2= eta[1]/2.0; //+1e-10;
		      boost::dynamic_bitset<> label(nGoal, false); // n is the number of goals
		      for(rocs::UintSmall i = 0; i < nGoal; ++i) {
			  label[2-i] = (goal[i][0] <= (x[0]-c1) && (x[0]+c1) <= goal[i][1] &&
				      goal[i][2] <= (x[1]-c2) && (x[1]+c2) <= goal[i][3])
			      ? true: false;
		      }
		      return label.to_ulong();
		  };
    abst.assign_labels(label_target);
    std::cout << "Specification assignment is done.\n";


    /**
     * Solve a Buchi game on the product of NTS and DBA.
     */
    std::cout << "Start solving a Buchi game on the product of the abstraction and DBA...\n";
    rocs::BSolver solver;
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
    std::vector<std::string> tokens;
    boost::split(tokens, specfile, boost::is_any_of("."));
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
    ctlrWtr.write_array<double>(eta, carde::n, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(&(solver._sol));
    std::cout << "Controller writing is done.\n";

    std::cout << "Total time of used (abstraction+synthesis): " << tabst+tsyn << '\n';

    return 0;
}
