/**
 *  carDBA.cpp
 *
 *  A general main program performing DBA control synthesis.
 *  Based on the T operator: U_{j} pre(W_i | O_ij)
 *
 *  The input spec files are: dba1,2,3.txt,
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

// #include "src/matlabio.h"
#include "src/hdf5io.h"

#include "car.hpp"


int main(int argc, char *argv[])
{
    /**
     * Default values
     **/
    std::string specfile;
    double e[]{0.2, 0.2, 0.2}; /* partition precision */
    
    /**
     * Input arguments: (2<= argc <= 5)
     * carAbst dbafile precision(e.g. 0.2 0.2 0.2)
     */
    if (argc < 2 || argc > 5) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }
    specfile = std::string(argv[1]);
    if (argc > 2 && argc < 5) {
	std::cout << "Input precision should be of 3-dim, e.g. 0.2 0.2 0.2.\n";
	std::exit(1);
    }
    if (argc > 2) {
	for(int i = 2; i < 5; ++i)
	    e[i-2] = std::atof(argv[i]);
    }
    std::cout << "Partition precision: " << e[0] << ' '
	      << e[1] << ' ' << e[2] << '\n';
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
    rocs::UintSmall labels[]{4, 2, 1}; // corresponding to goal1,2,3.
    auto init_w = [&carvf, &nProps, &labels, &arrayM,
		   &olb, &oub, &nA,
		   &glb, &gub, &nG](std::vector<rocs::CSolver*> &w,
				    rocs::UintSmall i,
				    rocs::UintSmall oid[]) {
		      w[i] = new rocs::CSolver(&carvf, nProps, rocs::RELMAX);
		      for (rocs::UintSmall j = 0; j < nA; ++j)
			  w[i]->init(rocs::AVOID, olb[j], oub[j]); // avoid areas
		      for (rocs::UintSmall j = 0; j < nG; ++j)
			  w[i]->labeling(glb[j], gub[j], labels[j]); // labeled areas
		      w[i]->set_M(arrayM[oid[i]]);
		  };

    std::cout << "Start control synthesis...\n";
    // std::vector<std::string> tokens;
    // boost::split(tokens, specfile, boost::is_any_of("."));
    rocs::UintSmall oid[nNodes];
    for(int i = 0; i < nNodes; ++i)
	oid[i] = i;
    rocs::dba_control< rocs::DTCntlSys<carde> >(w, &carvf, sdoms, nNodes,
						isacc, init_w, oid, e);
    

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
