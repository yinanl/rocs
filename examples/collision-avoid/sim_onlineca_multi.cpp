/**
 *  sim_onlineca_multi.cpp
 *
 *  Simulate the LTL motion planning with collision avoidance
 *  with multiple moving obstacles.
 *
 *  Created by Yinan Li on Feb. 27, 2021.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <utility>
#include <cmath>
#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <cstdlib>
#include <ctime>

#include "src/grid.h"
#include "src/definitions.h"
#include "src/abstraction.hpp"
#include "src/DBAparser.h"
#include "src/bsolver.hpp"
#include "src/patcher.h"
#include "src/hdf5io.h"

#include "car.hpp"
#include "odes.hpp"


int main(int argc, char *argv[])
{
    std::string specfile, ctlrfile, cafile, graphfile;
    specfile = "dba1.txt";
    ctlrfile = "controller_dba1_0.2-0.2-0.2.h5";
    graphfile = "gwin.h5";
    const double eta[] = {0.2, 0.2, 0.2};

    cafile = "controller_safety_abst-0.8-0.9-0.1-0.1.h5";
    const double eta_r[] = {0.1, 0.1, 0.1};

    /* Input arguments:
     * ./simMulti dbafile ctlrfile cafile graphfile eta[0] eta_r[0]
     */
    if(argc!=1) {
	if(argc==7) {
	    specfile = std::string(argv[1]);
	    ctlrfile = std::string(argv[2]);
	    cafile = std::string(argv[3]);
	    graphfile = std::string(argv[4]);
	    eta[0] = std::atof(argv[5]);
	    eta[1] = std::atof(argv[5]);
	    eta[2] = std::atof(argv[5]);
	    eta_r[0] = std::atof(argv[6]);
	    eta_r[1] = std::atof(argv[6]);
	    eta_r[2] = std::atof(argv[6]);
	} else {
	    std::cout << "Improper number of arguments.\n";
	    std::exit(1);
	}
    }
    std::cout << "Simulate " << specfile << " with controller " << ctlrfile << '\n';


    /**
     * Setup the motion planning workspace
     **/
    const int x_dim = 3;
    const double theta = 3.5;
    const double xlb[] = {0, 0, -theta};
    const double xub[] = {10, 10, theta};
    /* set the control values */
    const int u_dim = 2;
    const double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    /* discretization precision */
    double mu[] = {0.3, 0.3};
    /* generate grid */
    rocs::grid x_grid(x_dim,eta,xlb,xub);
    x_grid.gridding();
    rocs::grid u_grid(u_dim,mu,ulb,uub);
    std::cout << "Number of discrete states: " << x_grid._nv << "\n";
    std::cout << "Number of discrete inputs: " << u_grid._nv << "\n";

    clock_t tb, te;


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
     * Load global controller
     **/
    std::cout << "\nLoading global controller...\n";
    std::vector<long long> w_x0, encode3;
    std::vector<NODE_POST> nts_ctrlr;
    std::vector<CTRL> ctrl;
    std::vector<int> q_prime;
    rocs::h5FileHandler planRdr(ctlrfile, H5F_ACC_RDONLY);
    planRdr.read_discrete_controller(w_x0, encode3, nts_ctrlr, ctrl, q_prime);


    /**
     * Load local safety controller
     **/
    std::cout << "\nLoading local safety controller...\n";
    rocs::h5FileHandler reader(cafile, H5F_ACC_RDONLY);
    std::vector<size_t> optCtlr;
    std::vector<double> value;
    boost::dynamic_bitset<> safeCtlr;
    boost::dynamic_bitset<> win;
    size_t cdims[2];
    reader.read_discrete_controller(win, safeCtlr, cdims, optCtlr, value);
    double xrlb[] = {-3, -3, -theta};
    double xrub[] = {3, 3, theta};
    rocs::grid rel_grid(3, eta_r, xrlb, xrub);
    rel_grid.gridding();



    /**
     * Launch patcher
     **/
    rocs::Patcher local;
    struct stat buffer;
    if(stat(graphfile.c_str(), &buffer) == 0) {
    	/* Read from a file */
    	std::cout << "\nReading winning graph...\n";
    	rocs::h5FileHandler graphRdr(graphfile, H5F_ACC_RDONLY);
    	tb = clock();
    	graphRdr.read_winning_graph(local);
    	te = clock();
    	std::cout << "Time of reading graph: " << (float)(te - tb)/CLOCKS_PER_SEC << '\n';
    } else {
    	std::cout << "Graph file doesn't exist.\n";
    	return 1;
    }
    int horizon = 3;  //set forward propagation horizon



    /**
     * Simulation
     **/
    const double drange[] = {3.0, 3.0};
    /* set ode solver */
    const double tsim = 0.1; //simulation time
    const double dt = 0.001; //integration step size for odeint
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;
    boost::dynamic_bitset<> u_safe(cdims[1]);

    /* Initial condition of the ego robot */
    rocs::Rn x{1, 1, M_PI/3.0};
    rocs::Rn u{0, 0};
    long long x_index = x_grid.val_to_id(x);
    std::vector<long long>::iterator p1;
    p1 = std::lower_bound(w_x0.begin(), w_x0.end(), x_index);
    if(*p1 != x_index) {
	std::cout << x[0] <<  " "  << x[1] << " " << x[2];
	std::cout << " is not in the winning set, change a new initial state." << std::endl;
	return 1;
    }
    /* Initial condition of 3 obstacles */
    rocs::Rn uomax{0.8, 0.8};
    rocs::Rn uomin{-0.8, -0.8};
    double v;
    
    /* State in local relative coordinate */
    rocs::Rn xr1{0,0,0};
    rocs::Rn xr2{0,0,0};
    rocs::Rn xr3{0,0,0};
    long long xr_id1, xr_id2, xr_id3;

    /* obstacle 1 */
    rocs::Rn xo1{3.0, 10.0, -1.5};
    rocs::Rn uo1{0.1, 0.0};
    auto obstacle_behavior1 = [&uo1, tsim](int j) {
				  double t = tsim*j;
				  if(t >= 0) {
				      uo1[0] = 0.2;
				      uo1[1] = 0.0;
				  }
			      };
    

    /* obstacle 2 */
    rocs::Rn xo2{5.3, 0.0, 2.0};
    rocs::Rn uo2{0.5, 0.0};
    auto obstacle_behavior2 = [&uo2, tsim](int j) {
				  double t = tsim*j;
				  if(t >= 0) {
				      if(t < 14) {
					  uo2[0] = 0.5;
					  uo2[1] = 0.0;
				      } else {
					  uo2[0] = 0.5;
					  uo2[1] = 0.6;
				      }
				  }
			      };

    /* obstacle 3 */
    rocs::Rn xo3{10.0, 6.4, M_PI};
    rocs::Rn uo3{0.15, 0.0};
    auto obstacle_behavior3 = [&uo3, tsim](int j) {
				  double t = tsim*j;
				  if(t >= 0) {
				      if(t < 31) {
					  uo3[0] = 0.2;
					  uo3[1] = 0.0;
				      } else if(t < 35) {
					  uo3[0] = 0.2;
					  uo3[1] = 0.3;
				      } else {
					  uo3[0] = 0.2;
					  uo3[1] = 0.0;
				      }
				  }
			      };


    /* Open files for writing results and logs */
    std::string simfile = "traj_closedloop_ca_4.txt";
    std::ofstream ctlrWtr(simfile);
    if(!ctlrWtr.is_open())
	ctlrWtr.open(simfile, std::ios::out);
    std::string simother1 = "traj_other_robot_4_1.txt";
    std::string simother2 = "traj_other_robot_4_2.txt";
    std::string simother3 = "traj_other_robot_4_3.txt";
    std::ofstream otherWtr1(simother1);
    std::ofstream otherWtr2(simother2);
    std::ofstream otherWtr3(simother3);
    if(!otherWtr1.is_open())
    	otherWtr1.open(simother1, std::ios::out);
    if(!otherWtr2.is_open())
    	otherWtr2.open(simother2, std::ios::out);
    if(!otherWtr3.is_open())
    	otherWtr3.open(simother3, std::ios::out);
    /********** Logging **********/
    std::ofstream logger;
    logger.open("logs_4.txt", std::ios::out);
    /********** Logging **********/


    /* Simulation loop */
    srand(time(NULL));
    float tpat = 0;
    int max_num_achieve_acc=5, max_num_iteration=1500; //3000000;
    rocs::UintSmall q;
    int i, j;
    double t0, tc; //record the time when an obstacle detected
    bool o1, o2, o3;
    o1 = false;
    o2 = false;
    o3 = false;
    std::vector<CTRL>::iterator p7;
    NODE_POST p5;
    std::cout << "\nLaunching simulation...\n";
    for(q = q0, i = j = 0; i<max_num_achieve_acc && j<max_num_iteration; ++j) {
	/* choose obstacle trajectory */
	obstacle_behavior1(j);
	obstacle_behavior2(j);
	obstacle_behavior3(j);

	/* Determine the control u by the product state (x_index, q) */
	p5 = nts_ctrlr[encode3[x_index]];
	p7 = std::lower_bound(ctrl.begin()+p5.pos, ctrl.begin()+p5.pos+p5.num_a,
			      q, [](const CTRL &item, const int val) {
				     return item.q < val;});
	if(p7->q != q) {
	    // std::cout << j << "th iteration: " << "Reach accepting state " << i << " times.\n";
	    // std::cout << "Automaton state: " << q << '\n';
	    // std::cout << "Position in ctrl list: " << p5.pos << ", number of actions: " << p5.num_a << '\n';
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

	tc = j*tsim; //record the current time

	/* Get the relative state xr */
	xr1[0] = xo1[0]-x[0];
	xr1[1] = xo1[1]-x[1];
	xr1[2] = xo1[2]-x[2];
	xr2[0] = xo2[0]-x[0];
	xr2[1] = xo2[1]-x[1];
	xr2[2] = xo2[2]-x[2];
	xr3[0] = xo3[0]-x[0];
	xr3[1] = xo3[1]-x[1];
	xr3[2] = xo3[2]-x[2];
	if(xr1[2] > M_PI) //convert angle into [-pi, pi]
	    xr1[2] -= 2*M_PI;
	if(xr1[2] < -M_PI)
	    xr1[2] += 2*M_PI;
	if(xr2[2] > M_PI) //convert angle into [-pi, pi]
	    xr2[2] -= 2*M_PI;
	if(xr2[2] < -M_PI)
	    xr2[2] += 2*M_PI;
	if(xr3[2] > M_PI) //convert angle into [-pi, pi]
	    xr3[2] -= 2*M_PI;
	if(xr3[2] < -M_PI)
	    xr3[2] += 2*M_PI;
	inertia_to_body(xr1, x); //convert from the inertial to body frame
	inertia_to_body(xr2, x);
	inertia_to_body(xr3, x);

	/* Deal with obstacles */
	if(std::fabs(xr1[0])*std::fabs(xr1[0]) +
	   std::fabs(xr1[1])*std::fabs(xr1[1]) < drange[0]*drange[1]) {
	    o1 = true;
	    std::cout << "Obstacle 1 detected.\n";
	    std::cout << "Relative coordinate: ";
	    if(rel_grid._bds.isout(xr1)) { // xr is the out-of-domain node
		std::cout << rel_grid._nv;
	    } else {
		xr_id1 = rel_grid.val_to_id(xr1); //uniform grid rel_grid
		std::cout << xr_id1;
	    }
	    std::cout << "([" << xr1[0] << ',' << xr1[1] << ',' << xr1[2] << "])\n";
	} else {
	    o1 = false;
	}//end if obstacle is within range

	if(std::fabs(xr2[0])*std::fabs(xr2[0]) +
	   std::fabs(xr2[1])*std::fabs(xr2[1]) < drange[0]*drange[1]) {
	    o2 = true;
	    std::cout << "Obstacle 2 detected.\n";
	    std::cout << "Relative coordinate: ";
	    if(rel_grid._bds.isout(xr2)) { // xr is the out-of-domain node
		std::cout << rel_grid._nv;
	    } else {
		xr_id2 = rel_grid.val_to_id(xr2); //uniform grid rel_grid
		std::cout << xr_id2;
	    }
	    std::cout << "([" << xr2[0] << ',' << xr2[1] << ',' << xr2[2] << "])\n";
	} else {
	    o2 = false;
	}//end if obstacle is within range

	if(std::fabs(xr3[0])*std::fabs(xr3[0]) +
	   std::fabs(xr3[1])*std::fabs(xr3[1]) < drange[0]*drange[1]) {
	    o3 = true;
	    std::cout << "Obstacle 3 detected.\n";
	    std::cout << "Relative coordinate: ";
	    if(rel_grid._bds.isout(xr3)) { // xr is the out-of-domain node
		std::cout << rel_grid._nv;
	    } else {
		xr_id3 = rel_grid.val_to_id(xr3); //uniform grid rel_grid
		std::cout << xr_id3;
	    }
	    std::cout << "([" << xr3[0] << ',' << xr3[1] << ',' << xr3[2] << "])\n";
	} else {
	    o3 = false;
	}//end if obstacle is within range

	
	/* Perform RH (receding horizon) */
	if(o1||o2||o3) {
	    /********** Logging **********/
	    logger << tc << '(' << j << "): "
	       << x_index << '(' << x[0] << ',' << x[1] << ',' << x[2] << "), ";
	    if(o1) {
		logger << "o1," << xr_id1 << '(' << xr1[0] << ',' << xr1[1] << ',' << xr1[2] << "), ";
		for(size_t k = 0; k < cdims[1]; ++k)
		    logger << safeCtlr[xr_id1*cdims[1]+k];
		logger << '\n';
	    }
	    if(o2) {
		logger << "o2," << xr_id2 << '(' << xr2[0] << ',' << xr2[1] << ',' << xr2[2]<< "), ";
		for(size_t k = 0; k < cdims[1]; ++k)
		    logger << safeCtlr[xr_id2*cdims[1]+k];
		logger << '\n';
	    }
	    if(o3) {
		logger << "o3," << xr_id3 << '(' << xr3[0] << ',' << xr3[1] << ',' << xr3[2] << ")";
		for(size_t k = 0; k < cdims[1]; ++k)
		    logger << safeCtlr[xr_id3*cdims[1]+k];
		logger << '\n';
	    }
	    std::cout << "o1:"<< o1 <<", o2:"<< o2 <<", o3:" << o3 <<'\n';
	    // std::cout << "Current state id of the robot:" << xp << ", key:" << key << '\n';
	    /********** Logging **********/
	    
	    size_t xp = x_index*nNodes+q; // (idx,idq)=idx*nNodes+idq (buchi.c)
	    size_t uid;
	    
	    /********** Logging **********/
	    size_t row, key = local._encode[xp];
	    row = local._idmap[key];
	    // std::cout << "id_n0n1: " << xp << ',' << " id_np: " << key << '\n';
	    for(size_t k = 0; k < u_grid._nv; ++k) {
		if(local._winfts._npost[row*local._na+k]) {//outdeg>0
		    logger << k << ' ';
		}
	    }
	    logger << ", ";
	    /********** Logging **********/
	    
	    /* Determine safe inputs for RH module */
	    std::cout << "Use RH to find safe control inputs.\n ";
	    u_safe.set();
	    int nSafeAct = 0;
	    rocs::Rn uu(u_dim);
	    double dtemp, u_dist=(uub[0]-ulb[0])*(uub[0]-ulb[0])+(uub[1]-ulb[1])*(uub[1]-ulb[1]);
	    for(size_t k = 0; k < cdims[1]; ++k) {
		if(o1)
		    u_safe[k] = u_safe[k] & safeCtlr[xr_id1*cdims[1]+k];
		if(o2)
		    u_safe[k] = u_safe[k] & safeCtlr[xr_id2*cdims[1]+k];
		if(o3)
		    u_safe[k] = u_safe[k] & safeCtlr[xr_id3*cdims[1]+k];
		
		if(u_safe[k]) {
		    /************* Logging **************/
		    // std::cout << k << ',';
		    // std::cout << "Outdeg: "
		    // 	      << local._winfts._npost[row*local._na+k]
		    // 	      << '\n';
		    /************* Logging **************/
		    if(local._winfts._npost[row*local._na+k]) {//outdeg>0
			++nSafeAct; //count the # of safe valid controls
			/* take the closest to the original one */
			u_grid.id_to_val(uu, k);
			dtemp = (uu[0]-u[0])*(uu[0]-u[0])+(uu[1]-u[1])*(uu[1]-u[1]);
			if(dtemp < u_dist) {
			    uid = k;
			    u_dist= dtemp;
			}
			/************* Logging **************/
			// std::cout << "Post nodes: ";
			// size_t pidstart=local._winfts._ptrpost[row*local._na+k];
			// for(size_t pid=pidstart; pid<pidstart+local._winfts._npost[row*local._na+k]; ++pid)
			//     std::cout << local._winfts._idpost[pid] << ' ';
			// std::cout << '\n';
			/************* Logging **************/
		    }
		    // /************* Logging **************/
		    // logger << k << ' ';
		    // /************* Logging **************/
		}
		/********** Logging **********/
		logger << u_safe[k];
		/********** Logging **********/
	    }
	    /********** Logging **********/
	    u_grid.id_to_val(uu, uid);
	    logger << '\n'
		   << p7->u << ' ' << u[0] << ',' << u[1] << ';'
		   << uid << ' ' << uu[0] << ',' << uu[1] << ' ' << u_dist << ';';
	    /********** Logging **********/

	    /* Get a safe control from patcher */
	    if(u_safe[p7->u]) {
		std::cout << "The static strategy is safe.\n";
		// u_grid.id_to_val(u, uid);
		logger << p7->u << ' ';
	    } else if(nSafeAct > 1) {
		tb = clock();
		int suc = local.solve_local_reachability(xp, horizon, u_safe);
		te = clock();
		tpat += (float)(te - tb)/CLOCKS_PER_SEC;
		std::cout << "Average time for RH: "
			  << tpat/(j+1) << ".\n";
		if(!suc) {//if fails in re-planning, use the closest u
		    std::cout << "Apply the safe control input closest to the original one.\n";
		    u_grid.id_to_val(u, uid);
		    logger << uid << ' ';
		} else {
		    std::cout << "Apply controls returned by local reachability.\n";
		    u_grid.id_to_val(u, local._ctlr);
		    logger << local._ctlr << ' ';
		}
	    } else if(nSafeAct > 0) {
		std::cout << "Apply the only one safe control.\n";
		u_grid.id_to_val(u, uid);
		logger << uid << ' ';
	    } else {
		std::cout << "No safe controls. Freeze.\n";
		u[0] = 0;
		u[1] = 0;
		logger << p7->u << ' ';
	    }
	    /********** Logging **********/
	    logger << u[0] << ',' << u[1] << "\n\n";
	    /********** Logging **********/
	}//end RH
	

	/* Write to txt file */
	// ctlrWtr << j << ':';
	for(int d = 0; d < x_dim; ++d) {
	    /* the controlled state */
	    ctlrWtr << x[d];
	    /* the state of the other robot */
	    otherWtr1 << xo1[d];
	    otherWtr2 << xo2[d];
	    otherWtr3 << xo3[d];
	    if(d<x_dim-1) {
		ctlrWtr << ',';
		otherWtr1 << ',';
		otherWtr2 << ',';
		otherWtr3 << ',';
	    }
	}
	ctlrWtr << ';';
	otherWtr1 << ';';
	otherWtr2 << ';';
	otherWtr3 << ';';
	for(int d = 0; d < u_dim; ++d) {
	    ctlrWtr << u[d];
	    otherWtr1 << uo1[d];
	    otherWtr2 << uo2[d];
	    otherWtr3 << uo3[d];
	    if(d<u_dim-1) {
		ctlrWtr << ',';
		otherWtr1 << ',';
		otherWtr2 << ',';
		otherWtr3 << ',';
	    }
	}
	ctlrWtr << '\n';
	otherWtr1 << '\n';
	otherWtr2 << '\n';
	otherWtr3 << '\n';

	/* Integrate ego robot trajectory */
	boost::numeric::odeint::integrate_const(rk45, car_dynamics(u), x, 0.0, tsim, dt);
	if(x[2] > M_PI) //convert angle into [-pi, pi]
	    x[2] -= 2*M_PI;
	if(x[2] < -M_PI)
	    x[2] += 2*M_PI;

	/* integrate obstacle trajectory after it appears */
	boost::numeric::odeint::integrate_const(rk45, car_dynamics(uo1), xo1, 0.0, tsim, dt);
	if(xo1[2] > M_PI) //convert angle into [-pi, pi]
	    xo1[2] -= 2*M_PI;
	if(xo1[2] < -M_PI)
	    xo1[2] += 2*M_PI;
	boost::numeric::odeint::integrate_const(rk45, car_dynamics(uo2), xo2, 0.0, tsim, dt);
	if(xo2[2] > M_PI) //convert angle into [-pi, pi]
	    xo2[2] -= 2*M_PI;
	if(xo2[2] < -M_PI)
	    xo2[2] += 2*M_PI;
	boost::numeric::odeint::integrate_const(rk45, car_dynamics(uo3), xo3, 0.0, tsim, dt);
	if(xo3[2] > M_PI) //convert angle into [-pi, pi]
	    xo3[2] -= 2*M_PI;
	if(xo3[2] < -M_PI)
	    xo3[2] += 2*M_PI;

	/* Update automaton state q */
	x_index = x_grid.val_to_id(x);
	q = q_prime[p5.label*nNodes + q];
	if(isacc[q])
	    std::cout << "******** achieve acc \'" << q << "\' " << ++i << " time ********" << std::endl;
    }

    ctlrWtr.close();
    otherWtr1.close();
    otherWtr2.close();
    otherWtr3.close();
    logger.close();

    return 0;
}
