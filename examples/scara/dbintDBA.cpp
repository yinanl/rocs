/**
 *  dba_integrator.cpp
 *
 *  DBA control of the simplified SCARA manipulator dynamics (the double integrator model).
 *
 *  Created by Yinan Li on July 14, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <cmath>

#include "src/DBAparser.h"
#include "src/system.hpp"
#include "src/csolver.h"

// #include "src/matlabio.h"
#include "src/hdf5io.h"

#include "scara.hpp"


const double h = 0.8*l1;
const double r = 0.5*l1;
double a1 = atan(h / r); // the upper bound for theta1
double a2 = asin(h / l1);
/**
 * Define constraints for theta2 (for current r, h) in the form of
 * f(x)\leq 0
 *
 */
template<typename T>
T collision1(const T &x) {
    // T y(2);
    // y[1] = x[0] + x[1] - M_PI + atan((h - l1*sin(x[0])) / (l1*cos(x[0]) - r));
    // y[0] = x[0] - a2;
    T y{x[0]-a2, -x[0]-x[1]+M_PI-atan( (h-l1*sin(x[0]))/(l1*cos(x[0])-r) )};
    return y;
}

template<typename T>
T collision2(const T &x) {
    // T y(3);
    // y[0] = x[0] - a1;
    // y[1] = a2 - x[0];
    // y[2] = x[0] + x[1] - M_PI + atan((l1*sin(x[0])-h) / (l1*cos(x[0])));
    T y{x[0]-a1, a2-x[0],
	-x[0]-x[1]+M_PI+atan( (l1*sin(x[0])-h)/(l1*cos(x[0])) )};
    return y;
}


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
     * Setup the control system.
     */
    double xlb[] = {0, -M_PI, -1, -1};
    double xub[] = {M_PI/2.0, M_PI, 1, 1};
    double ulb[] = {-5.0, -5.0};
    double uub[] = {5.0, 5.0};
    double mu[] = {0.5, 0.5};

    rocs::DTCntlSys<integrator> scara("dba", ts, integrator::n, integrator::m);
    scara.init_workspace(xlb, xub);
    scara.init_inputset(mu, ulb, uub);


    // /**
    //  * Convert from the operational space to the joint space
    //  */
    // rocs::ivec xy1{rocs::interval(0.03, 0.07), rocs::interval(0.18, 0.22)};
    // rocs::ivec xy2{rocs::interval(0.25, 0.28), rocs::interval(0.02, 0.05)};
    // rocs::ivec theta1(4), theta2(4);
    // xy2theta(theta1, xy1);
    // xy2theta(theta2, xy2);
    // std::cout << "Target region 2: ";
    // for(int i = 0; i < 4; ++i)
    // 	std::cout << theta1[i] << ' ';
    // std::cout << '\n';
    // std::cout << "Target region 2: ";
    // for(int i = 0; i < 4; ++i)
    // 	std::cout << theta2[i] << ' ';
    // std::cout << '\n';


    /**
     * DBA control synthesis
     */
    /* Initialize the set of S-domains */
    std::vector<rocs::CSolver*> w(nNodes);
    std::vector<rocs::SPtree*> sdoms(nNodes);
    /* Set avoid area */
    double obs[2][4] = {{a1, xlb[1], xlb[2], xlb[3]},
    			{M_PI/2.0, xub[1], xub[2], xub[3]}};

    const double e[]{0.01, 0.01, 3, 3};  // only bisect x[0] and x[1].
    for (rocs::UintSmall i = 0; i < nNodes; ++i) {
	w[i] = new rocs::CSolver(&scara, nProps, rocs::RELMAX);
	w[i]->init(rocs::AVOID, obs[0], obs[1]);
	w[i]->init(rocs::AVOID, &collision1<rocs::ivec>, e);
	w[i]->init(rocs::AVOID, &collision2<rocs::ivec>, e);
	w[i]->init_avoid_area();
    }
    /* Assign labels */
    rocs::UintSmall nG = 2;
    /* Two goals:
     * (x=0.05, y=0.2), the joint space (0.5126, 1.6264)
     * (x=0.27, y=0.03), the joint space (0.5488, -0.8763)
     */
    // double goal[2][4] = {{0.4080, 1.6172, -0.1, -0.1},
    // 			 {0.6172, 1.6355, 0.1, 0.1}};
    double glb[][4] = {{0.4980, 1.5739, -0.1, -0.1},
		       {0.4903, -0.9363, -0.1, -0.1}};
    double gub[][4] = {{0.5772, 1.7055, 0.1, 0.1},
		       {0.6069, -0.8363, 0.1, 0.1}};
    rocs::UintSmall labels[]{1, 2}; // corresponding to goal1,2,3.
    for (rocs::UintSmall i = 0; i < nNodes; ++ i) {
	w[i]->set_M(arrayM[i]);
	for (rocs::UintSmall j = 0; j < nG; ++j)
	    w[i]->labeling(glb[j], gub[j], labels[j]);
	sdoms[i] = &(w[i]->_ctlr);
	w[i]->init_goal_area();
    }

    double eta[]{0.05, 0.05, 0.1, 0.1};
    rocs::dba_control< rocs::DTCntlSys<integrator> >(w,&scara,sdoms,nNodes,isacc,eta);


    /**
     * Display and save memoryless controllers.
     */
    for (rocs::UintSmall i = 0; i < w.size(); ++i) {
	std::cout << std::endl << "S-domain for q" << std::to_string(i) << '\n';
	w[i]->print_controller_info();
    }
    //rocs::write_results_to_mat(scara, specfile, w);
    rocs::write_csolvers_to_h5(scara, specfile, w);


    /**
     * Release dynamic memory
     **/
    for (rocs::UintSmall i = 0; i < nNodes; ++i)
	delete w[i];


    return 0;
}
