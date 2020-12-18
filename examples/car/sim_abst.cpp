/**
 *  sim_carAbst.cpp
 *
 *  To simulate the vehicle example using abstraction-based control.
 *
 *  Authors: Yinan Li, Zhibing Sun, Jun Liu
 *  Created: May 27, 2020
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include <string>
#include <math.h>
#include <boost/numeric/odeint.hpp>

#include "src/abstraction.hpp"
#include "src/DBAparser.h"
#include "src/bsolver.hpp"
#include "src/hdf5io.h"


/* define dynamics */
struct car_dynamics {
	rocs::Rn u;
	car_dynamics (const rocs::Rn param): u (param) {}

	void operator() (rocs::Rn &x, rocs::Rn &dxdt, double t) const
	{
	      double alpha = std::atan(std::tan(u[1])/2.0);
	      dxdt[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
	      dxdt[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
	      dxdt[2] = u[0]*std::tan(u[1]);
	}
};


int main(int argc, char *argv[])
{
    /**
     * Default arguments
     */
    std::string specfile{"dba1.txt"};
    double eta[]{0.2, 0.2, 0.2};

    /* Input arguments:
     * sim_abst specfile precision(e.g. 0.2 0.2 0.2)
     */
    if (argc > 2 && argc < 5) {
	std::cout << "Improper number of arguments. Input arguments are:\n";
	std::cout << "./sim_abst specfile precision(e.g. 0.2 0.2 0.2)\n";
	std::exit(1);
    }
    if(argc > 1) {
	specfile = std::string(argv[1]);
    }
    if (argc == 5) {
	for(int i = 2; i < 5; ++i)
	    eta[i-2] = std::atof(argv[i]);
    }
    std::vector<std::string> tokens;
    boost::split(tokens, specfile, boost::is_any_of("."));
    std::string ctlrfile = "controller_" + tokens[0] + "_";
    for(int i = 0; i < 3; ++i) {
	std::stringstream ss;
	ss << std::setprecision(1);
	ss << eta[i];
	ctlrfile += ss.str();
	if (i < 2)
	    ctlrfile += "-";
    }
    ctlrfile += ".h5";


    /**
     * Load specification
     **/
    std::cout << "\nReading the specification...\n";
    rocs::UintSmall nAP = 0, nNodes = 0, q0 = 0;
    std::vector<rocs::UintSmall> acc;
    std::vector<std::vector<rocs::UintSmall>> arrayM;
    if (!rocs::read_spec(specfile, nNodes, nAP, q0, arrayM, acc))
	std::exit(1);
    boost::dynamic_bitset<> isacc(nNodes, false);
    for (rocs::UintSmall i = 0; i < acc.size(); ++i)
	isacc[acc[i]] = true;


    /**
     * Setup the motion planning workspace
     **/
    /* set the state space */
    const int x_dim = 3;
    const double theta = 3.5;
    double xlb[] = {0, 0, -theta};
    double xub[] = {10, 10, theta};
    /* set the control values */
    const int u_dim = 2;
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};
    /* Generate grids */
    rocs::grid x_grid(x_dim,eta,xlb,xub);
    x_grid.gridding();
    rocs::grid u_grid(u_dim,mu,ulb,uub);
    std::cout << "Number of discrete states: " << x_grid._nv << "\n";
    std::cout << "Number of discrete inputs: " << u_grid._nv << "\n";


    /**
     * Load global controller
     **/
    std::cout << "\nLoading the controller from " << ctlrfile << "...\n";
    std::vector<long long> w_x0, encode3;
    std::vector<NODE_POST> nts_ctrlr;
    std::vector<CTRL> ctrl;
    std::vector<int> q_prime;
    rocs::h5FileHandler planRdr(ctlrfile, H5F_ACC_RDONLY);
    planRdr.read_discrete_controller(w_x0, encode3, nts_ctrlr, ctrl, q_prime);


    /**
     * Simulation
     **/
    const double tsim = 0.3;  //simulation time step
    /* Open files for writing results and logs */
    std::string simfile = "sim_traj_"+tokens[0]+".txt";
    std::ofstream ctlrWtr(simfile);
    if(!ctlrWtr.is_open())
	ctlrWtr.open(simfile, std::ios::out);
	/********** Logging **********/
    std::ofstream logger;
    logger.open("logs_"+tokens[0]+".txt", std::ios::out);
    /********** Logging **********/

    /* ode solver */
    const double dt = 0.001; //integration step size for odeint
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;

    /* Initial condition */
    rocs::Rn x{3, 2, M_PI/2.};
    rocs::Rn u{0, 0};
    long long x_index = x_grid.val_to_id(x);
    std::vector<long long>::iterator p1;
    p1 = std::lower_bound(w_x0.begin(), w_x0.end(), x_index);
    if(*p1 != x_index) {
	std::cout << x[0] <<  " "  << x[1] << " " << x[2];
	std::cout << " is not in the winning set, change a new initial state." << std::endl;
	return 1;
    }

    float tpat = 0;
    int max_num_achieve_acc=5, max_num_iteration=500; //3000000;
    rocs::UintSmall q;
    int i, j;
    std::vector<CTRL>::iterator p7;
    NODE_POST p5;
    std::cout << "\nLaunching simulation...\n";
    for(q = q0, i = j = 0; i<max_num_achieve_acc && j<max_num_iteration; ++j) {
	/* Determine the control u by the product state (x_index, q) */
	p5 = nts_ctrlr[encode3[x_index]];
	p7 = std::lower_bound(ctrl.begin()+p5.pos, ctrl.begin()+p5.pos+p5.num_a,
			      q, [](const CTRL &item, const int val) {
				     return item.q < val;});
	if(p7->q != q) {
	    std::cout << j << "th iteration: " << "Reach accepting state " << i << " times.\n";
	    std::cout << "Automaton state: " << q << '\n';
	    std::cout << "Position in ctrl list: " << p5.pos << ", number of actions: " << p5.num_a << '\n';
	    std::cout << "Error in ctrl" << std::endl;
	    break;
	}
	u_grid.id_to_val(u, p7->u);

	/* Print to screen */
	std::cout << "# " << j << ":\n";
	std::cout << "x: [" << x[0] <<  ','  << x[1] << ',' << x[2] << "]\n";
	std::cout << "Nominal control: " << p7->u << '('
		  << '[' << u[0] << ',' << u[1] << "])\n";
	std::cout << "current dba state: " << q << "\n";

	/* Write to txt file */
	// ctlrWtr << j << ':';
	for(int d = 0; d < x_dim; ++d) {
	    ctlrWtr << x[d];
	    if(d<x_dim-1) {
		ctlrWtr << ',';
	    }
	}
	ctlrWtr << ';';
	for(int d = 0; d < u_dim; ++d) {
	    ctlrWtr << u[d];
	    if(d<u_dim-1) {
		ctlrWtr << ',';
	    }
	}
	ctlrWtr << '\n';

	/* Integrate the dynamics */
	boost::numeric::odeint::integrate_const(rk45, car_dynamics(u), x, 0.0, tsim, dt);
	if(x[2] > M_PI) //convert angle into [-pi, pi]
	    x[2] -= 2*M_PI;
	if(x[2] < -M_PI)
	    x[2] += 2*M_PI;

	/********** Logging **********/
    logger << "x,q,u=" << x_index << ' ' << q << ' ' << p7->u << "; encode3[x]=" << encode3[x_index];
	logger << "\np5=" << p5.num_a << ' ' << p5.label << ' ' << p5.pos;
	logger << "\np7=";
	for(std::vector<CTRL>::iterator p = ctrl.begin()+p5.pos;
		p!=ctrl.begin()+p5.pos+p5.num_a; ++p) {
		logger << p->q << ' ' << p->u << ';';
	}
	logger << "\n\n";
    /********** Logging **********/

	/* Update automaton state q */
	x_index = x_grid.val_to_id(x);
	q = q_prime[p5.label*nNodes + q];
	if(isacc[q])
	    std::cout << "******** achieve acc \'" << q << "\' " << ++i << " time ********" << std::endl;
    }

    ctlrWtr.close();
	logger.close();


    return 0;
}
