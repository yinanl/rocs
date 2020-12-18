/**
 *  carDBAall.cpp
 *
 *  A general main program performing DBA control synthesis.
 *  Based on the T operator: U_{j} pre(W_i | O_ij)
 *
 *  The input spec files are: dba1,2,3,4.txt,
 *  in which the propositions are not simplified.
 *  
 *  Created by Yinan Li on Jan. 30, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "src/DBAparser.h"
#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"
#include "src/hdf5io.h"

#include "car.hpp"


int main(int argc, char *argv[])
{
    /**
     * Default values
     **/
    std::string specfile;
    double e[]{0.2, 0.2, 0.2}; /* partition precision */
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
	    e[i-2] = std::atof(argv[i]);
	
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
    std::cout << "Partition precision: " << e[0] << ' '
	      << e[1] << ' ' << e[2] << '\n';
    
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
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc)) 
	std::exit(1);
    size_t nProps = arrayM.size();
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
    rocs::UintSmall nA = 4;
    double olb[4][3] = {{1.6, 4.0, -theta},
			{3.0, 5.0, -theta},
			{4.3, 1.8, -theta},
			{5.7, 1.8, -theta}};
    double oub[4][3] = {{5.7, 5.0, theta},
			{5.0, 8.0, theta},
			{5.7, 4.0, theta},
			{8.5, 2.5, theta}};
    /* Goals */
    rocs::UintSmall nG = 3; // # of goals
    double glb[3][3]= {{1.0, 0.5, -theta},
		       {0.5, 7.5, -theta},
		       {7.1, 4.6, -theta}};
    double gub[3][3]= {{2.0, 2.0, theta},
		       {2.5, 8.5, theta},
		       {9.1, 6.4, theta}};
    

    /**
     * DBA control synthesis
     */
    /* Initialize the set of S-domains */
    std::vector<rocs::CSolver*> w(nNodes);
    std::vector<rocs::SPtree*> sdoms(nNodes);
    /* Set avoid area */
    for (rocs::UintSmall i = 0; i < nNodes; ++i) {
	w[i] = new rocs::CSolver(&carvf, nProps, rocs::RELMAX);
	for (rocs::UintSmall j = 0; j < nA; ++j)
	    w[i]->init(rocs::AVOID, olb[j], oub[j]);
	// for (size_t j = 0; j < obs.size(); ++j)
	//     w[i]->init(rocs::AVOID, obs[j]);
    }
    rocs::UintSmall labels[]{4, 2, 1}; // corresponding to goal1,2,3.
    for (rocs::UintSmall i = 0; i < nNodes; ++ i) {
	w[i]->set_M(arrayM[i]);
	std::cout << "Outgoig transitions of w" << i <<":\n";
	for (auto &m : w[i]->_M)
	    std::cout << m << ' ';
	std::cout << '\n';
	for (rocs::UintSmall j = 0; j < nG; ++j) {
	    // w[i]->labeling(bset[j], labels[j]);
	    w[i]->labeling(glb[j], gub[j], labels[j]);
	}
	sdoms[i] = &(w[i]->_ctlr);
	std::cout << "Initial partition of w" << i <<":\n";
	w[i]->print_controller();
    }

    std::cout << "Start control synthesis...\n";
    std::vector<std::string> tokens;
    boost::split(tokens, specfile, boost::is_any_of("."));
    if(preprocess && tokens[0]=="dba1") {/* TO BE AUTOMATED */
	/* order: q4 --> q1,0 --> q3 --> q2 */
	/* W 4 */
	rocs::UintSmall n0 = 1;
	std::vector<rocs::CSolver*> wsub0(n0);
	wsub0[0] = w[4];
	boost::dynamic_bitset<> accsub0(n0, false);
	accsub0[0] = isacc[4];
	rocs::dba_control< rocs::DTCntlSys<carde> >(wsub0, &carvf, sdoms, n0, accsub0, e);
	/* W 2, 3 */
	rocs::UintSmall n1 = 2;
	std::vector<rocs::CSolver*> wsub1(n1);
	wsub1[0] = w[0];
	wsub1[1] = w[1];
	boost::dynamic_bitset<> accsub1(n1, false);
	accsub1[0] = isacc[0];
	accsub1[1] = isacc[1];
	rocs::dba_control< rocs::DTCntlSys<carde> >(wsub1, &carvf, sdoms, n1, accsub1, e);
	/* W 1 */
	rocs::UintSmall n2 = 1;
	std::vector<rocs::CSolver*> wsub2(n2);
	wsub2[0] = w[3];
	boost::dynamic_bitset<> accsub2(n2, false);
	accsub2[0] = isacc[3];
	rocs::dba_control< rocs::DTCntlSys<carde> >(wsub2, &carvf, sdoms, n2, accsub2, e);
	/* W 0 */
	rocs::UintSmall n3 = 1;
	std::vector<rocs::CSolver*> wsub3(n3);
	wsub3[0] = w[2];
	boost::dynamic_bitset<> accsub3(n3, false);
	accsub3[0] = isacc[2];
	rocs::dba_control< rocs::DTCntlSys<carde> >(wsub3, &carvf, sdoms, n3, accsub3, e);
    } else {
	rocs::dba_control< rocs::DTCntlSys<carde> >(w, &carvf, sdoms, nNodes, isacc, e);
    }
    

    /**
     * Display and save memoryless controllers.
     */
    for (rocs::UintSmall i = 0; i < w.size(); ++i) {
	std::cout << std::endl << "S-domain for q" << std::to_string(i) << '\n';
	w[i]->print_controller_info();
    }
    // rocs::write_results_to_mat(carvf, specfile, w);
    rocs::write_csolvers_to_h5(carvf, specfile, w);
    

    /**
     * Release dynamic memory 
     **/
    for (rocs::UintSmall i = 0; i < nNodes; ++i)
	delete w[i];
    
    return 0;
}
