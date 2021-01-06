/**
 *  engine2dDBA.cpp
 *
 *  Reach and stay control of the ode model of Moore-Greitzer engine.
 *
 *  Created by Yinan Li on June 9, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "src/system.hpp"
#include "src/csolver.h"
#include "src/DBAparser.h"

// #include "src/matlabio.h"
#include "src/hdf5io.h"


double a = 1./3.5;
double B = 2.0;
double H = 0.18;
double W = 0.25;
double lc = 8.0;
double cx = 1.0/lc;
double cy = 1.0/(4*lc*B*B);
double aH = a+H;
double H2 = H/(2.0*W*W*W);
double W2 = 3*W*W;

/* user defined dynamics */
struct mgode2 {
    static const int n = 2;  // system dimension
    static const int nu = 2;  // control dimension

    // double a{1.67};
    // double B{2.0};
    // double H{0.18};
    // double W{0.25};
    // double lc{8.0};
    // double cx{1.0/lc};
    // double cy{1.0/(4*lc*B*B)};
    // double aH{a+H};
    // double H2{H/(2.0*W*W*W)};
    // double W2{3*W*W};
    
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
     * carDBAall dbafile
     */
    if (argc != 2) {
	std::cout << "Improper number of arguments.\n";
	std::exit(1);
    }
    
    
    /**
     * Read specification file
     **/
    std::cout << "Loading the specification...\n";
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    std::string specfile = std::string(argv[1]);
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
     * Control system setup
     **/
    /* set the state space */
    double xlb[]{0.44, 0.6};
    double xub[]{0.54, 0.7};
    
    /* set the control values */
    double Lmu = 0.01;
    double Lu = 0.05;
    double ulb[2]{-Lu, 0.5};
    double uub[2]{Lu, 0.8};
    double mu[2]{Lu/5, 0.01};

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

    /* define the target areas */
    /* case I */
    double e = 0.003;
    double glb[][2]{{0.5039-e, 0.6605-e}};
    double gub[][2]{{0.5039+e, 0.6605+e}};
    double olb[][2]{{0.520, 0.658}};
    double oub[][2]{{0.526, 0.664}};
    // /* case II */
    // double e = 0.003;
    // double glb[][2] {{0.4519-e, 0.6513-e}};
    // double gub[][2] {{0.4519+e, 0.6513+e}};
    // double olb[][2] {{0.497, 0.650}};
    // double oub[][2] {{0.503, 0.656}};
    rocs::UintSmall nG = 1, nA = 1;


    /**
     * DBA control synthesis
     */
    
    /* Initialize the set of S-domains */
    std::vector<rocs::CSolver*> w(nNodes);
    std::vector<rocs::SPtree*> sdoms(nNodes);
    rocs::UintSmall labels[]{1}; // label the goal area
    auto init_w = [&engine, nProps, &labels, &arrayM,
		   &olb, &oub, nA,
		   &glb, &gub, nG](std::vector<rocs::CSolver*> &w,
				   rocs::UintSmall i,
				   rocs::UintSmall oid[]) {
		      w[i] = new rocs::CSolver(&engine, nProps, rocs::RELMAX);
		      for (rocs::UintSmall j = 0; j < nA; ++j)
			  w[i]->init(rocs::AVOID, olb[j], oub[j]); // avoid areas
		      for (rocs::UintSmall j = 0; j < nG; ++j)
			  w[i]->labeling(glb[j], gub[j], labels[j]); // labeled areas
		      w[i]->set_M(arrayM[oid[i]]);
		  };
    /* Perform synthesis */
    // double ei[]{0.0002, 0.0002};
    double er[]{0.00018, 0.00018};
    rocs::UintSmall oid[nNodes];
    for(int i = 0; i < nNodes; ++i)
	oid[i] = i;
    rocs::dba_control< rocs::CTCntlSys<mgode2> >(w,&engine,sdoms,nNodes,isacc,
						 init_w, oid, er);
    

    /**
     * Display and save memoryless controllers.
     */
    for (rocs::UintSmall i = 0; i < w.size(); ++i) {
	std::cout << std::endl << "S-domain for q" << std::to_string(i) << '\n';
	w[i]->print_controller_info();
    }
    // rocs::write_results_to_mat(engine, specfile, w);
    rocs::write_csolvers_to_h5(engine, specfile, w);
    

    /**
     * Release dynamic memory 
     **/
    for (rocs::UintSmall i = 0; i < nNodes; ++i)
	delete w[i];
    
    
    engine.release_flows();
    
    return 0;
}
