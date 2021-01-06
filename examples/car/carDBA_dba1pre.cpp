/**
 *  carDBA.cpp
 *
 *  An optimized DBA control synthesis for general reachability.
 *  
 *  Created by Yinan Li on Jan. 4, 2021.
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

//#include "src/matlabio.h"
#include "src/hdf5io.h"

#include "car.hpp"


int main(int argc, char *argv[])
{
    /**
     * Default values
     **/
    std::string specfile = "dba1.txt";
    double e[]{0.2, 0.2, 0.2}; /* partition precision */
    
    /**
     * Input arguments: (argc = 1 or 4)
     * carAbst precision(e.g. 0.2 0.2 0.2)
     */
    if (argc > 1 && argc < 4) {
	std::cout << "Input precision should be of 3-dim, e.g. 0.2 0.2 0.2.\n";
	std::exit(1);
    }
    
    if (argc > 1) {
	for(int i = 1; i < 4; ++i)
	    e[i-1] = std::atof(argv[i]);
    }
    std::cout << "Partition precision: " << e[0] << ' '
	      << e[1] << ' ' << e[2] << '\n';
    std::cout << "Using preprocessing.\n";
    

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
    auto init_w = [&carvf, nProps, &olb, &oub, nA,
    		   &glb, &gub, nG, &labels, &arrayM](std::vector<rocs::CSolver*> &w,
						     rocs::UintSmall i,
						     rocs::UintSmall oid[]) {
    		      w[i] = new rocs::CSolver(&carvf, nProps, rocs::RELMAX);
    		      for (rocs::UintSmall j = 0; j < nA; ++j)
    			  w[i]->init(rocs::AVOID, olb[j], oub[j]); // avoid areas
    		      for (rocs::UintSmall j = 0; j < nG; ++j)
    			  w[i]->labeling(glb[j], gub[j], labels[j]); // labeled areas
    		      w[i]->set_M(arrayM[oid[i]]);
		      // /***** LOGGING  *****/
		      // w[i]->print_controller();
		      // /***** LOGGING  *****/
    		  };
    
    std::cout << "Start control synthesis...\n";
    std::vector<std::string> tokens;
    boost::split(tokens, specfile, boost::is_any_of("."));
    // if(preprocess && tokens[0]=="dba1") {/* TO BE AUTOMATED */
    /* order: q4 --> q1,0 --> q3 --> q2 */
    /* W 4 */
    rocs::UintSmall n0 = 1;
    std::vector<rocs::CSolver*> w0(n0);
    boost::dynamic_bitset<> acc0(n0, false);
    acc0[0] = isacc[4];
    rocs::UintSmall oid0[]{4};
    // rocs::dba_control< rocs::DTCntlSys<carde> >(w0, &carvf, sdoms, n0, acc0,
    // 						init_w, oid0, e);
    init_w(w0, 0, oid0);
    sdoms[oid0[0]] = &(w0[0]->_ctlr);
    w0[0]->init_winset();
    w[4] = w0[0];
    /* W 0, 1 */
    rocs::UintSmall n1 = 2;
    std::vector<rocs::CSolver*> w1(n1);
    boost::dynamic_bitset<> acc1(n1, false);
    acc1[0] = isacc[0];
    acc1[1] = isacc[1];
    rocs::UintSmall oid1[]{0, 1};
    rocs::dba_control< rocs::DTCntlSys<carde> >(w1, &carvf, sdoms, n1, acc1,
						init_w, oid1, e);
    w[0] = w1[0];
    w[1] = w1[1];
    /* W 3 */
    rocs::UintSmall n2 = 1;
    std::vector<rocs::CSolver*> w2(n2);
    boost::dynamic_bitset<> acc2(n2, false);
    acc2[0] = isacc[3];
    rocs::UintSmall oid2[]{3};
    rocs::dba_control< rocs::DTCntlSys<carde> >(w2, &carvf, sdoms, n2, acc2,
						init_w, oid2, e);
    w[3] = w2[0];
    /* W 2 */
    rocs::UintSmall n3 = 1;
    std::vector<rocs::CSolver*> w3(n3);
    boost::dynamic_bitset<> acc3(n3, false);
    acc3[0] = isacc[2];
    rocs::UintSmall oid3[]{2};
    rocs::dba_control< rocs::DTCntlSys<carde> >(w3, &carvf, sdoms, n3, acc3,
						init_w, oid3, e);
    w[2] = w3[0];

    

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
