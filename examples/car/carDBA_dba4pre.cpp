/**
 *  carDBA_dba4pre.cpp
 *
 *  Control synthesis for dba4 with preprocessing. (TO BE AUTOMATED)
 *
 *  Created by Yinan Li on Jan. 5, 2020.
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
    std::string specfile = "dba4.txt";
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
    rocs::UintSmall nA = 3;
    double olb[3][3] = {{0.0, 4.0, -theta},
			{5.4, 5.0, -theta},
			{4.5, 0.0, -theta}};
    double oub[3][3] = {{3.2, 5.0, theta},
			{6.0, 10.0, theta},
			{5.2, 2.5, theta}};
    /* Goals */
    rocs::UintSmall nG = 6;
    double glb[6][3]= {{1.0, 0.5, -theta}, //a1
		       {3.8, 3.1, -theta}, //a4
		       {0.0, 6.0, -theta}, //c
		       {6.0, 0.0, -theta}, //d
		       {0.5, 7.5, -theta}, //a2
		       {7.1, 1.9, -theta}};//a3
    double gub[6][3]= {{2.0, 2.0, theta}, //a1
		       {4.6, 4.5, theta}, //a4
		       {4.0, 10.0, theta}, //c
		       {10.0, 4.0, theta}, //d
		       {2.5, 8.5, theta}, //a2
		       {9.1, 2.9, theta}};//a3


    /**
     * DBA control synthesis
     */
    /* Initialize the set of S-domains */
    std::vector<rocs::CSolver*> w(nNodes);
    std::vector<rocs::SPtree*> sdoms(nNodes);
    rocs::UintSmall labels[]{32, 1, 8, 4, 10, 20}; //a1(32),a3(20,=a3&d),c(8),d(4),a2(10,=a2&c),a4(1)
    auto init_w = [&carvf, nProps, &labels, &arrayM,
		   &olb, &oub, nA,
		   &glb, &gub, nG](std::vector<rocs::CSolver*> &w,
				   rocs::UintSmall i,
				   rocs::UintSmall oid[]) {
		      w[i] = new rocs::CSolver(&carvf, nProps, rocs::RELMAX);
		      for (rocs::UintSmall j = 0; j < nA; ++j)
			  w[i]->init(rocs::AVOID, olb[j], oub[j]); // avoid areas
		      for (rocs::UintSmall j = 0; j < nG; ++j)
			  w[i]->labeling(glb[j], gub[j], labels[j]); // labeled areas
		      w[i]->set_M(arrayM[oid[i]]);
		  };


    /* Perform synthesis */
    /* order: q2,3,4 --> q5 --> q1,2 */
    /* W 2,3,4 */
    rocs::UintSmall n1 = 3;
    std::vector<rocs::CSolver*> w1(n1);
    boost::dynamic_bitset<> acc1(n1, false);
    acc1[0] = isacc[2];
    acc1[1] = isacc[3];
    acc1[2] = isacc[4];
    rocs::UintSmall oid1[]{2, 3, 4};
    rocs::dba_control< rocs::DTCntlSys<carde> >(w1, &carvf, sdoms, n1, acc1,
						init_w, oid1, e);
    w[2] = w1[0];
    w[3] = w1[1];
    w[4] = w1[2];
    /* W 5 */
    rocs::UintSmall n2 = 1;
    std::vector<rocs::CSolver*> w2(n2);
    boost::dynamic_bitset<> acc2(n2, false);
    acc2[0] = isacc[5];
    rocs::UintSmall oid2[]{5};
    rocs::dba_control< rocs::DTCntlSys<carde> >(w2, &carvf, sdoms, n2, acc2,
						init_w, oid2, e);
    w[5] = w2[0];
    /* W 0 1 */
    rocs::UintSmall n3 = 2;
    std::vector<rocs::CSolver*> w3(n3);
    boost::dynamic_bitset<> acc3(n3, false);
    acc3[0] = isacc[0];
    acc3[1] = isacc[1];
    rocs::UintSmall oid3[]{0, 1};
    rocs::dba_control< rocs::DTCntlSys<carde> >(w3, &carvf, sdoms, n3, acc3,
						init_w, oid3, e);
    w[0] = w3[0];
    w[1] = w3[1];


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
