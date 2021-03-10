/**
 *  sim_onlineca_replan.cpp
 *
 *  Simulate the LTL motion planning with collision avoidance
 *  with moving obstacles.
 *  Additional replanning if the robot is blocked.
 *
 *  Created by Yinan Li on Feb. 19, 2021.
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
     * ./simRe (dbafile ctlrfile cafile graphfile eta[0] eta_r[0])
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
     * Target and avoid region selection for replanning
     */
    rocs::UintSmall q = 0;
    const double ro = 0.8;
    const double rt = 0.35;
    double per = 0.8;
    std::vector<long long> avoid, target;
    std::vector< std::pair<size_t, int> > replan;
    boost::dynamic_bitset<> plhdr(cdims[1], true);
    auto local_region = [&x_grid, q, nNodes](const rocs::Rn &xo, const double ro,
					     const std::vector<long long> &encode) {
			    /* collect the grid points intersect with the box */
			    rocs::ivec box{rocs::interval(xo[0]-ro, xo[0]+ro),
					   rocs::interval(xo[1]-ro, xo[1]+ro),
					   x_grid._bds[2]};
			    // std::cout << "The box enclose the region: " << box << '\n';
			    std::vector<size_t> ids = x_grid.subset(box, 0, 0);
			    std::vector<long long> ix(ids.begin(), ids.end());

			    /* determine the grid cell intersect with the r-ball around xo */
			    double xl, xr, yr, yl, xsqr, ysqr;
			    double rx = x_grid._gw[0]/2.0;
			    double ry = x_grid._gw[1]/2.0;
			    std::vector<double> x(x_grid._dim);
			    for(auto &id:ix) {
				x_grid.id_to_val(x, id);
				// std::cout << x[0] << ' ' << x[1] << ' ' << x[2] <<'\n';
				xl = x[0]-xo[0] - rx;
				xr = x[0]-xo[0] + rx;
				yl = x[1]-xo[1] - ry;
				yr = x[1]-xo[1] + ry;
				xsqr = (xr*xr) > (xl*xl) ? (xl*xl) : (xr*xr);
				ysqr = (yr*yr) > (yl*yl) ? (yl*yl) : (yr*yr);
				if(xsqr+ysqr <= ro*ro) {
				    // std::cout << id <<',';
				    id = encode[id*nNodes+q]; //n0xn1->np
				} else
				    id = -1;
			    }

			    return std::move(ix);
			};

    auto local_target = [&local, &local_region]
	(const rocs::Rn &xr, const double ro, const double rt,
	 const rocs::Rn &x0, const std::vector<long long> &encode) {
			    if(xr[0]<0) {
				return std::vector<long long>();
			    }
			    double drsqr = xr[0]*xr[0]+xr[1]*xr[1];
			    double dr = std::sqrt(drsqr);
			    double dt = 2*ro;
			    double a = std::asin(ro/dr);
			    double b = std::atan2(xr[1], xr[0]);
			    double c = std::atan2(dr, dt);

			    double rho, theta;
			    if(c > a) {
				theta = c;
				rho = std::sqrt(drsqr+dt*dt);
			    } else {
				theta = a;
				rho = dr/std::cos(a);
			    }
			    // double rho1 = 2*std::sqrt(drsqr-ro*ro);
			    // double rho2 = dr/std::cos(a);
			    // double rho = rho1>rho2 ? rho1:rho2;

			    double xt1 = rho*std::cos(b+theta);
			    double yt1 = rho*std::sin(b+theta);
			    double xt2 = rho*std::cos(b-theta);
			    double yt2 = rho*std::sin(b-theta);
			    rocs::Rn z1{xt1, yt1, 0};
			    rocs::Rn z2{xt2, yt2, 0};
			    /********** Logging **********/
			    std::cout << dr << ',' << b << ',' << a << ',' << c << ','
				      << theta << ',' << rho << '\n';
			    std::cout << "z1=(" << xt1 << ' ' << yt1 << ' ' << 0;
			    std::cout << ", z2=(" << xt2 << ' ' << yt2 << ' ' << 0 << ")\n";
			    /********** Logging **********/

			    body_to_inertia(z1, x0);
			    body_to_inertia(z2, x0);
			    std::cout << "Converting to inertia: ";
			    std::cout << "z1=(";
			    for(int i = 0; i < 3; ++i) {
				z1[i] += x0[i];
				std::cout << z1[i] << ' ';
			    }
			    std::cout << "), z2=(";
			    for(int i = 0; i < 3; ++i) {
				z2[i] += x0[i];
				std::cout << z2[i] << ' ';
			    }
			    std::cout << ")\n";

			    /********** Logging **********/
			    std::ofstream logger;
			    logger.open("logs_newgoals.txt", std::ios::out);
			    for(auto &i:z1)
			    	logger << i << ' ';
			    logger << '\n';
			    for(auto &i:z2)
			    	logger << i << ' ';
			    logger.close();
			    /********** Logging **********/

			    std::vector<long long> target1 = local_region(z1, rt, encode);
			    std::vector<long long> target2 = local_region(z2, rt, encode);
			    int n1 = std::count_if(target1.begin(), target1.end(),
						   [&local](long long i){
						       return i>-1?local._idmap[i]:0;
						   });
			    int n2 = std::count_if(target2.begin(), target2.end(),
						   [&local](long long i){
						       return i>-1?local._idmap[i]:0;
						   });

			    std::cout << n1 << ',' << n2 << '\n';

			    if(n1 >= n2) // return the target set with more winning nodes
				return std::move(target1);
			    else
				return std::move(target2);
			};


    /**
     * Simulation
     **/
    const double drange[] = {3.0, 3.0};
    /* set ode solver */
    const double tsim = 0.1; //simulation time
    const double dt = 0.001; //integration step size for odeint
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;

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
    /* Initial condition of the other robot */
    rocs::Rn uomax{0.8, 0.8};
    rocs::Rn uomin{-0.8, -0.8};
    double v;
    rocs::Rn xo{0.0, 3.0, 0};
    rocs::Rn uo{0.0, 0.0};
    /* State in local relative coordinate */
    rocs::Rn xr{0,0,0};
    long long xr_id;
    boost::dynamic_bitset<> u_safe(cdims[1]);
    rocs::Rn ui{0.3, 0.0};

    /* Test case 1: a designed pattern for the obstacle */
    auto obstacle_behavior1 = [&uo, &xo, tsim](int j) {
				  double t = tsim*j;
				  if(t >= 0) {
				      if(t == 0) {
					  rocs::Rn xi{0.0, 3.0, 0};
					  xo.assign(xi.begin(), xi.end());
				      }
				      if(t < 1.5) {
					  uo[0] = 0.33;
					  uo[1] = 0.0;
				      } else if(t < 3.0) {
					  uo[0] = 0.8;
					  uo[1] = 0.7;
				      } else if(t < 10.8) {
					  uo[0] = 0.7;
					  uo[1] = 0;
				      } else if(t < 18.0) {
					  uo[0] = 0.7;
					  uo[1] = -0.6;
				      } else if(t < 21) {
					  uo[0] = 0.7;
					  uo[1] = 0.6;
				      } else {
					  uo[0] = 0.5;
					  uo[1] = 0.0;
				      }
				  }
			      };

    /* Test case 2: circular trajectory for the obstacle */
    auto obstacle_behavior2 = [&uo, &xo, tsim](int j) {
				  double t = tsim*j;
				  if(t >= 0) {
				      if(t == 0) {
					  rocs::Rn xi{4.5, 0.0, M_PI/2.0};
					  xo.assign(xi.begin(), xi.end());
				      }
				      if(t < 8.1) {
					  uo[0] = 0.65;
					  uo[1] = 0.0;
				      } else if(t < 9.0) {
					  uo[0] = 0.4;
					  uo[1] = -0.3;
				      } else {
					  uo[0] = 0.4;
					  uo[1] = 0.3;
				      }
				  }
			      };

    /* Test case 3: random v, w for the obstacle */
    double th = 0;
    auto obstacle_behavior3 = [&th,&uo, &xo, tsim, &uomin, &uomax](int j) {
				  double t = tsim*j;
				  double t1 = 15.9;
				  double t2 = 2.1;
				  if(t >= 0) {
				      if(t == 0) {
					  rocs::Rn xi{10.0, 6.4, M_PI};
					  xo.assign(xi.begin(), xi.end());
				      }
				      if(t < t1) {
					  uo[0] = 0.4;
					  uo[1] = 0.0;
				      } else if(t < t1+t2) {
					  uo[0] = 0.4;
					  uo[1] = 0.3;
				      } else {
					  uo[1] = 0.0;
					  if(th>=8*tsim) {
					      uo[0] = -uo[0];
					      th = 0;
					  }
					  th += tsim;
				      }
				      // if(t < 18.9) { // straight line, constant speed
				      // 	  uo[0] = 0.4;
				      // 	  uo[1] = 0.0;
				      // } else { // random v, w
				      // 	  uo[0] = (uomax[0]-uomin[0])*((double)rand()/RAND_MAX)
				      // 	      - (uomax[0]-uomin[0])/2;
				      // 	  uo[1] = (uomax[1]-uomin[1])*((double)rand()/RAND_MAX)
				      // 	      - (uomax[1]-uomin[1])/2;
				      // }
				  }
			      };


    /* Open files for writing results and logs */
    std::string simfile = "traj_closedloop_ca_3_t.txt";
    std::ofstream ctlrWtr(simfile);
    if(!ctlrWtr.is_open())
	ctlrWtr.open(simfile, std::ios::out);
    std::string simother = "traj_other_robot_3_t.txt";
    std::ofstream otherWtr(simother);
    if(!otherWtr.is_open())
    	otherWtr.open(simother, std::ios::out);
    /********** Logging **********/
    std::ofstream logger;
    logger.open("logs_3_t.txt", std::ios::out);
    /********** Logging **********/


    /* Simulation loop */
    srand(time(NULL));
    float tpat = 0;
    int max_num_achieve_acc=5, max_num_iteration=1500; //3000000;
    // rocs::UintSmall q;
    int i, j;
    double t0, tc; //record the time when an obstacle detected
    bool entry = false;
    rocs::Rn xb(xr); //record the position of the robot when an obstacle detected
    bool replanned = false;
    std::vector<CTRL>::iterator p7;
    NODE_POST p5;
    std::cout << "\nLaunching simulation...\n";
    for(q = q0, i = j = 0; i<max_num_achieve_acc && j<max_num_iteration; ++j) {
	/* choose obstacle trajectory */
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


	/* Get the relative state xr */
	xr[0] = xo[0]-x[0];
	xr[1] = xo[1]-x[1];
	xr[2] = xo[2]-x[2];
	if(xr[2] > M_PI) //convert angle into [-pi, pi]
	    xr[2] -= 2*M_PI;
	if(xr[2] < -M_PI)
	    xr[2] += 2*M_PI;
	inertia_to_body(xr, x); //convert from the inertial to body frame

	/* Deal with obstacles */
	if(std::fabs(xr[0])*std::fabs(xr[0]) +
	   std::fabs(xr[1])*std::fabs(xr[1]) < drange[0]*drange[1]) {
	// if(std::fabs(xr[0])<drange[0] && std::fabs(xr[1])<drange[1]) {
	    if(!entry) {
		entry = true;
		t0 = j*tsim; //record the first time an obstacle is detected
		xb = x;
		/********** Logging **********/
		logger << "t0=" << t0 << ",xb="
		       << '(' << xb[0] << ',' << xb[1] << ',' << xb[2] << ")\n";
		/********** Logging **********/
	    }
	    tc = j*tsim; //record the current time

	    std::cout << "An obstacle detected.\n";
	    std::cout << "Relative coordinate: ";
	    if(rel_grid._bds.isout(xr)) { // xr is the out-of-domain node
		std::cout << rel_grid._nv;
	    } else {
		xr_id = rel_grid.val_to_id(xr); //uniform grid rel_grid
		std::cout << xr_id;
	    }
	    std::cout << "([" << xr[0] << ',' << xr[1] << ',' << xr[2] << "])\n";

	    /********** Logging **********/
	    logger << tc << '(' << j << "):"
		   << x_index << '(' << x[0] << ',' << x[1] << ',' << x[2] << "),"
		   << xr_id << '(' << xr[0] << ',' << xr[1] << ',' << xr[2] << ")\n";
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

	    /* Determine if a replanning is needed to avoid stuck */
	    if(tc-t0>3) {
		if(std::fabs(x[0]-xb[0])<0.4 && std::fabs(x[1]-xb[1])<0.4) {
		    if(!replanned) {//replanning
			std::cout << "Replanning...\n";
			avoid = local_region(xo, ro, local._encode); // np-based id's
			target = local_target(xr, ro, rt, x, local._encode);
			if(!target.empty()) {
			    // /********** Logging **********/
			    // rocs::h5FileHandler wtr("logs_replan.h5", H5F_ACC_TRUNC);
			    // long long ix;
			    // std::vector<long long> wa, wt;
			    // for(auto id:avoid) {
			    // 	if(id>-1?local._idmap[id]:0) {
			    // 	    ix = local._decode[id]/nNodes;
			    // 	    wa.push_back(ix);
			    // 	}
			    // }
			    // wtr.write_array<long long>(wa, "avoid");
			    // for(auto id:target) {
			    // 	if(id>-1?local._idmap[id]:0) {
			    // 	    ix = local._decode[id]/nNodes;
			    // 	    wt.push_back(ix);
			    // 	}
			    // }
			    // wtr.write_array<long long>(wt, "target");
			    // wtr.write_array<double>(x, "x");
			    // wtr.write_array<double>(xo, "xo");

			    // // std::cin.get();
			    // /********** Logging **********/
			    tb = clock();
			    replan = local.solve_local_reachavoid(xp, per*target.size(), avoid, target, plhdr);
			    te = clock();
			    tpat += (float)(te - tb)/CLOCKS_PER_SEC;
			    std::cout << "Average time for Detouring: "
				      << tpat/(j+1) << ".\n";
			    if(replan[0].first == xp)
				replanned = true;
			    else
				std::cout << "Replanning failed.\n";
			} else {
			    std::cout << "Replanning failed: no target points.\n";
			}
		    }
		} else {
		    entry = false; //reset t0 if robot make progress in 3 secs
		}
	    }
	    /* Perform replan if replan is activated. Replan is activated when
	     * - replanning is successful, and
	     * - local target hasn't reached. */
	    std::vector<std::pair<size_t, int> >::iterator uiter;
	    if(replanned && std::find(target.begin(),target.end(),local._encode[xp])==target.end()) {
		/* If replan succeeds and the local target hasn't been reached */
		uiter = std::find_if(replan.begin(), replan.end(),
				     [&xp](std::pair<size_t, int> &item){return item.first==xp;});
		if(uiter != replan.end()) {//replan can be followed correctly
		    std::cout << "Apply the replanned reach control.\n";
		    uid = uiter->second;
		    u_grid.id_to_val(u, uid);
		    /********** Logging **********/
		    logger << uid << ' ';
		    /********** Logging **********/
		} else {
		    replanned = false;
		    std::cout << "Error in executing replan. Abort replan and follow RH.\n";
		}
	    } else {
		replanned = false;
	    }

	    /* Perform RH (receding horizon) module if replan is not activated. */
	    if(!replanned) {
		/* Determine safe inputs for RH module */
		std::cout << "Use RH to find safe control inputs: ";
		int nSafeAct = 0;
		rocs::Rn uu(u_dim);
		double dtemp, u_dist=(uub[0]-ulb[0])*(uub[0]-ulb[0])+(uub[1]-ulb[1])*(uub[1]-ulb[1]);
		for(size_t k = 0; k < cdims[1]; ++k) {
		    u_safe[k] = safeCtlr[xr_id*cdims[1]+k];
		    /************* Logging **************/
		    if(u_safe[k]) {
			// std::cout << k << ',';
			// std::cout << "Outdeg: "
			// 	      << local._winfts._npost[row*local._na+k]
			// 	      << '\n';
			if(local._winfts._npost[row*local._na+k]) {//outdeg>0
			    ++nSafeAct; //count the # of safe valid controls
			    /* take the closest to the original one */
			    u_grid.id_to_val(uu, k);
			    dtemp = (uu[0]-u[0])*(uu[0]-u[0])+(uu[1]-u[1])*(uu[1]-u[1]);
			    if(dtemp < u_dist) {
				uid = k;
				u_dist= dtemp;
			    }
			    // std::cout << "Post nodes: ";
			    // size_t pidstart=local._winfts._ptrpost[row*local._na+k];
			    // for(size_t pid=pidstart; pid<pidstart+local._winfts._npost[row*local._na+k]; ++pid)
			    //     std::cout << local._winfts._idpost[pid] << ' ';
			    // std::cout << '\n';
			}
			logger << k << ' ';
		    }
		    /************* Logging **************/
		}
		/********** Logging **********/
		u_grid.id_to_val(uu, uid);
		logger << '\n'
		       << p7->u << ' ' << u[0] << ',' << u[1] << ';'
		       << uid << ' ' << uu[0] << ',' << uu[1] << ' ' << u_dist << ';';
		/********** Logging **********/

		/* Get a safe control from patcher */
		if(u_safe[p7->u]) {
		// if(safeCtlr[xr_id*cdims[1]+uid]) {
		    std::cout << "The static strategy is safe.\n";
		    // u_grid.id_to_val(u, uid);
		    // logger << uid << ' ';
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
	    }//end RH
	    /********** Logging **********/
	    logger << u[0] << ',' << u[1] << "\n\n";
	    /********** Logging **********/

	} else {
	    entry = false;
	    replanned = false; //reset replanning flag
	}//end if obstacle is within range

	/* Write to txt file */
	// ctlrWtr << j << ':';
	for(int d = 0; d < x_dim; ++d) {
	    /* the controlled state */
	    ctlrWtr << x[d];
	    /* the state of the other robot */
	    otherWtr << xo[d];
	    if(d<x_dim-1) {
		ctlrWtr << ',';
		otherWtr << ',';
	    }
	}
	ctlrWtr << ';';
	otherWtr << ';';
	for(int d = 0; d < u_dim; ++d) {
	    ctlrWtr << u[d];
	    otherWtr << uo[d];
	    if(d<u_dim-1) {
		ctlrWtr << ',';
		otherWtr << ',';
	    }
	}
	ctlrWtr << '\n';
	otherWtr << '\n';

	/* Integrate ego robot trajectory */
	boost::numeric::odeint::integrate_const(rk45, car_dynamics(u), x, 0.0, tsim, dt);
	if(x[2] > M_PI) //convert angle into [-pi, pi]
	    x[2] -= 2*M_PI;
	if(x[2] < -M_PI)
	    x[2] += 2*M_PI;

	/* integrate obstacle trajectory after it appears */
	boost::numeric::odeint::integrate_const(rk45, car_dynamics(uo), xo, 0.0, tsim, dt);
	if(xo[2] > M_PI) //convert angle into [-pi, pi]
	    xo[2] -= 2*M_PI;
	if(xo[2] < -M_PI)
	    xo[2] += 2*M_PI;

	/* Update automaton state q */
	x_index = x_grid.val_to_id(x);
	q = q_prime[p5.label*nNodes + q];
	if(isacc[q])
	    std::cout << "******** achieve acc \'" << q << "\' " << ++i << " time ********" << std::endl;
    }

    ctlrWtr.close();
    otherWtr.close();
    logger.close();

    return 0;
}
