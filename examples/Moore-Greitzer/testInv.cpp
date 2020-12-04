/**
 *  testInv.cpp
 *
 *  Compare abstraction-based and interval methods on solving invariance problems.
 *  
 *  Created by Yinan Li on Dec. 02, 2020.
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

#include "src/txtfileio.h"
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
     * carAbst dbafile 
     */
    if (argc != 2) {
    	std::cout << "Improper number of arguments.\n";
    	std::exit(1);
    }
    std::string specfile = std::string(argv[1]);
    
    clock_t tb, te;
    /* set the state space */
    double e0 = 0.004;
    double xlb[]{0.5039-e0, 0.6605-e0};
    double xub[]{0.5039+e0, 0.6605+e0};
    
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

    /* The goal area */
    double e = 0.003;
    double goal[][2]{{0.5039-e, 0.5039+e}, {0.6605-e, 0.6605+e}};
    

    /**
     * Read DBA from spec *.txt file
     */
    std::cout << "Reading the specification...\n";
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
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


    double eta[]{0.00018, 0.00018};

    /**
     * Abstraction-based control
     */
    std::string transfile = "abst_goal_0.00018.h5";
    // std::string ctlrabst = "controller_Gb_abst.txt";
    std::string h5ctlrabst = "controller_Gb_abst_0.00018.h5";
    
    rocs::abstraction< rocs::CTCntlSys<mgode2> > abst(&engine);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    abst.assign_label_outofdomain(0);
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
    	double e1[] = {0.0, 0.0};
    	double e2[] = {0.0, 0.0};
    	tb = clock();
    	abst.assign_transitions(e1, e2);
    	te = clock();
    	tabst = (float)(te - tb)/CLOCKS_PER_SEC;
    	std::cout << "Time of computing abstraction: " << tabst << '\n';
    	std::cout << "# of transitions: " << abst._ts._ntrans << '\n';
    	/* Write abstraction to file */
    	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
    	transWtr.write_transitions(abst._ts);
    }
    /* Assign labels to states: has to be consistent with the dba file */
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
    /* Solve a Buchi game on the product of NTS and DBA */
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
    /* Display and save memoryless controllers */
    std::cout << "Writing the controller...\n";
    // solver.write_controller_to_txt(const_cast<char*>(ctlrabst.c_str()));
    rocs::h5FileHandler ctlrWtr(h5ctlrabst, H5F_ACC_TRUNC);
    ctlrWtr.write_problem_setting< rocs::CTCntlSys<mgode2> >(engine);
    ctlrWtr.write_2d_array<double>(targetPts, "G");
    ctlrWtr.write_array<double>(eta, mgode2::n, "eta");
    ctlrWtr.write_2d_array<double>(abst._x._data, "xgrid");
    ctlrWtr.write_discrete_controller(&(solver._sol));
    std::cout << "Controller writing is done.\n";

    std::cout << "Total time of used (abstraction+synthesis): " << tabst+tsyn << '\n';

    /* Test read discrete controller */
    std::vector<long long> w_x0;
    std::vector<long long> encode3;
    std::vector<NODE_POST> nts_ctrlr;
    std::vector<CTRL> ctrl;
    std::vector<int> q_prime;
    rocs::h5FileHandler Rdr(h5ctlrabst, H5F_ACC_RDONLY);
    Rdr.read_discrete_controller(w_x0, encode3, nts_ctrlr, ctrl, q_prime);
    std::cout << "w_x0: " << w_x0.size() << '\n';
    std::cout << "encode3: " << encode3.size() << '\n';
    std::cout << "nts_ctrlr: " << nts_ctrlr.size() << '\n';
    std::cout << "ctrl: " << ctrl.size() << '\n';
    std::cout << "q_prime: " << q_prime.size() << '\n';



    /**
     * DBA control synthesis
     */
    // double glb[][2]{{0.5039-e, 0.6605-e}};
    // double gub[][2]{{0.5039+e, 0.6605+e}};
    // double glb[][2]{{0.50102, 0.65762}};
    // double gub[][2]{{0.50678, 0.66338}};
    // rocs::UintSmall nG = 1, nA = 1;
    // /* Initialize the set of S-domains */
    // std::vector<rocs::CSolver*> w(nNodes);
    // std::vector<rocs::SPtree*> sdoms(nNodes);
    // rocs::UintSmall labels[]{1}; // corresponding to goal1,2,3.
    // for (rocs::UintSmall i = 0; i < nNodes; ++ i) {
    // 	w[i] = new rocs::CSolver(&engine, nProps, rocs::RELMAX, 50);
    // 	w[i]->set_M(arrayM[i]);
    // 	std::cout << "Outgoig transitions of w" << i <<":\n";
    // 	for (auto &m : w[i]->_M)
    // 	    std::cout << m << ' ';
    // 	std::cout << '\n';
    // 	for (rocs::UintSmall j = 0; j < nG; ++j)
    // 	    w[i]->labeling(glb[j], gub[j], labels[j]);
    // 	sdoms[i] = &(w[i]->_ctlr);
    // 	std::cout << "Initial partition of w" << i <<":\n";
    // 	w[i]->print_controller();
    // }
    // // double ei[]{0.0002, 0.0002};
    // double er[]{0.00018, 0.00018};
    // rocs::dba_control< rocs::CTCntlSys<mgode2> >(w,&engine,sdoms,nNodes,isacc,er);
    // /* Save controllers to .h5 files */
    // rocs::write_csolvers_to_h5(engine, specfile, w);

    // /* Use invariance algorithm */
    // rocs::CSolver solver(&engine, 0, rocs::RELMAX);
    // solver.init(rocs::GOAL, glb[0], gub[0]);
    // solver.init_goal_area();
    // solver.invariance_control(&engine, eta);
    // solver.print_controller_info();
    // std::string filename = "controller_Gb_itvl.h5";
    // rocs::h5FileHandler wtr(filename, H5F_ACC_TRUNC);
    // wtr.write_problem_setting(engine);
    // wtr.write_sptree_controller(solver);
    

    // /* Test reachable set */
    // rocs::ivec G{rocs::interval(0.5039-e, 0.6605-e),
    // 		 rocs::interval(0.5039+e, 0.6605+e)};
    // // rocs::ivec G{rocs::interval(0.50102, 0.50678),
    // // 		 rocs::interval(0.65762, 0.66338)};
    // rocs::ivec g0{rocs::interval(0.50102, 0.5039),rocs::interval(0.65762, 0.6605)};
    // rocs::ivec g1{rocs::interval(0.5039, 0.50678),rocs::interval(0.65762, 0.6605)};
    // rocs::ivec g2{rocs::interval(0.50102, 0.5039),rocs::interval(0.6605, 0.66338)};
    // rocs::ivec g3{rocs::interval(0.5039, 0.50678),rocs::interval(0.6605, 0.66338)};
    // std::vector<rocs::ivec> x0;
    // x0.push_back(g0);
    // x0.push_back(g1);
    // x0.push_back(g2);
    // x0.push_back(g3);

    // std::cout << "Boxes in goal that stay in goal:\n";
    // std::vector<rocs::ivec> y(engine._ugrid._nv, rocs::ivec(2));
    // for(auto &x : x0) {
    // 	engine.get_reach_set(y, x);
    // 	for(size_t i = 0; i < y.size(); ++i) {
    // 	    if(G.isin(y[i])) {
    // 		std::cout << "x=" << x
    // 			  << ", y=" << y[i]
    // 			  << ", u=" << i << '\n';
    // 	    }
    // 	}
    // }

    engine.release_flows();
    return 0;
}
