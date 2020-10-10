/**
 *  carDBAall.cpp
 *
 *  A general main program performing DBA control synthesis.
 *  Based on the T operator: U_{j} pre(W_i | O_ij)
 *
 *  The input spec files are: p1,1_1,2,3,4,4_preprocess.txt,
 *  in which the propositions are simplified.
 *  
 *  Created by Yinan Li on Jan. 30, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"
#include "src/DBAparser.h"

#include "car.hpp"


int main(int argc, char *argv[])
{
    /**
     * Default values
     **/
    double e[]{0.2, 0.2, 0.2}; /* partition precision */
    bool preprocess = false;
    
    /* Input arguments: 
     * carDBAall dbafile precision -preprocess 
     */
    if (argc > 3 | argc < 2) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }
    // if (argc > 2) {
    // 	std::stringstream convert(argv[2]);
    // 	if (!(convert >> e)) {
    // 	    std::cout << "Wrong argument for precision.\n";
    // 	    std::exit(1);
    // 	}
    // 	if (argc > 3) {
    // 	    std::string param = std::string(argv[3]);
    // 	    if (param == "-preprocess")
    // 		preprocess = true;
    // 	    else {
    // 		std::cout << "Wrong argument for preprocessing.\n";
    // 		std::exit(1);
    // 	    }
    // 	}
    // }
    // std::cout << "Partition precision: " << e << '\n';
    if (preprocess)
	std::cout << "Using preprocessing.\n";
    else
	std::cout << "Not using preprocessing.\n";
    

    /**
     * Read specification file
     **/
    std::cout << "Loading the specification...\n";
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    std::string specfile = std::string(argv[1]);
    rocs::UintSmall nProps = 0, nNodes = 0;
    if (!rocs::read_spec(specfile, nNodes, nProps, arrayM, acc)) 
    	std::exit(1);
    
    /* Mark accepting states */
    boost::dynamic_bitset<> isacc(nNodes, false);
    for (rocs::UintSmall i = 0; i < acc.size(); ++i)
	isacc[acc[i]] = true;
    std::cout << "isacc= ";
    for (boost::dynamic_bitset<>::size_type i = 0; 
         i < isacc.size(); ++i) {
        std::cout << isacc[i] << ' ';
    }
    std::cout << '\n';
    
    
    /**
     * Workspace Setup
     **/
    /* State space */
    const double theta = 3.5;  // M_PI
    double xlb[] = {0, 0, -theta};
    double xub[] = {10, 10, theta};
    /* Control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};
    /* Control system */
    rocs::DTCntlSys<carde> carvf("DBA", h, carde::n, carde::m);
    carvf.init_workspace(xlb, xub);
    carvf.init_inputset(mu, ulb, uub);
    /* Obstacles */
    std::vector<rocs::ivec> obs(4);
    rocs::ivec A(3);
    A.setval(0, rocs::interval(1.6, 5.7));
    A.setval(1, rocs::interval(4.0, 5.0));
    A.setval(2, rocs::interval(-theta, theta));
    obs[0] = A;
    A.setval(0, rocs::interval(3.0, 5.0));
    A.setval(1, rocs::interval(5.0, 8.0));
    A.setval(2, rocs::interval(-theta, theta));
    obs[1] = A;
    A.setval(0, rocs::interval(4.3, 5.7));
    A.setval(1, rocs::interval(1.8, 4.0));
    A.setval(2, rocs::interval(-theta, theta));
    obs[2] = A;
    A.setval(0, rocs::interval(5.7, 8.5));
    A.setval(1, rocs::interval(1.8, 2.5));
    A.setval(2, rocs::interval(-theta, theta));
    obs[3] = A;
    /* Goals */
    int nG = 3; // # of goals
    std::vector<rocs::ivec> bset(3);
    rocs::ivec B(3);
    B.setval(0, rocs::interval(1.0, 2.0));
    B.setval(1, rocs::interval(0.5, 2.0));
    B.setval(2, rocs::interval(-theta, theta));
    bset[0] = B;
    B.setval(0, rocs::interval(0.5, 2.5));
    B.setval(1, rocs::interval(7.5, 8.5));
    B.setval(2, rocs::interval(-theta, theta));
    bset[1] = B;
    B.setval(0, rocs::interval(7.1, 9.1));
    B.setval(1, rocs::interval(4.6, 6.4));
    B.setval(2, rocs::interval(-theta, theta));
    bset[2] = B;
    

    /**
     * DBA control synthesis
     */
    /* Initialize the set of S-domains */
    std::vector<rocs::CSolver*> w(nNodes);
    std::vector<rocs::SPtree*> sdoms(nNodes);
    /* Set avoid area */
    for (rocs::UintSmall i = 0; i < nNodes; ++i) {
	w[i] = new rocs::CSolver(&carvf, nProps, rocs::RELMAX);
	for (size_t j = 0; j < obs.size(); ++j) {
	    w[i]->init(rocs::AVOID, obs[j]);
	}
    }
    /* Assign labels */
    // for (rocs::UintSmall i = 0; i < nNodes; ++ i) {
    // 	w[i]->set_M(arrayM[i]);
    // 	std::cout << "Outgoig transitions of w" << i <<":\n";
    // 	for (auto &m : w[i]->_M)
    // 	    std::cout << m << ' ';
    // 	std::cout << '\n';
    // 	for (rocs::UintSmall j = 1; j < nProps; ++j) {
    // 	    if (arrayM[i][j] < nNodes) {
    // 		w[i]->labeling(bset[j-1], j);
    // 	    }
    // 	}
    // 	sdoms[i] = &(w[i]->_ctlr);
    // 	// std::cout << "Initial partition of w" << i <<":\n";
    // 	// w[i]->print_controller();
    // }

    rocs::UintSmall labels[]{3, 2, 1, 3, 0};
    for (rocs::UintSmall i = 0; i < nNodes; ++ i) {
    	w[i]->set_M(arrayM[i]);
	std::cout << "Outgoig transitions of w" << i <<":\n";
	for (auto &m : w[i]->_M)
	    std::cout << m << ' ';
	std::cout << '\n';
    	for (rocs::UintSmall j = 1; j < nProps; ++j) {
    	    if (arrayM[i][j] < nNodes) {
    		w[i]->labeling(bset[labels[j]], j);
    	    }
    	}
    	sdoms[i] = &(w[i]->_ctlr);
    	std::cout << "Initial partition of w" << i <<":\n";
    	w[i]->print_controller();
    }
    
    rocs::dba_control< rocs::DTCntlSys<carde> >(w, &carvf, sdoms, nNodes, isacc, e);
    

    /**
     * Display and save memoryless controllers.
     */
    for (rocs::UintSmall i = 0; i < w.size(); ++i) {
	std::cout << std::endl << "S-domain for q" << std::to_string(i) << '\n';
	w[i]->print_controller_info();
    }
    rocs::write_results_to_mat(carvf, specfile, w);
    

    /**
     * Release dynamic memory 
     **/
    for (rocs::UintSmall i = 0; i < nNodes; ++i)
	delete w[i];
    
    return 0;
}
