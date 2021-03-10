/**
 *  testStuckAvoidance.cpp
 *
 *  Test the local reachability to avoid stuck.
 *
 *  Created by Yinan Li on Feb. 12, 2021.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <utility>

#include "src/grid.h"
#include "src/definitions.h"
#include "src/abstraction.hpp"
#include "src/bsolver.hpp"
#include "src/patcher.h"
#include "src/hdf5io.h"
#include "src/DBAparser.h"

#include "odes.hpp"


int main()
{
    std::string specfile, ctlrfile, cafile, graphfile;
    specfile = "dba1.txt";
    ctlrfile = "controller_dba1_0.2-0.2-0.2.h5";
    graphfile = "gwin_s2.h5";
    const double eta[] = {0.2, 0.2, 0.2};
    
    // cafile = "controller_safety_abst-0.8-1-0.3-0.2.h5";
    cafile = "controller_safety_abst-0.8-0.9-0.1-0.1.h5";
    const double eta_r[] = {0.1, 0.1, 0.1};

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

    /* Functions avoiding stuck */
    rocs::UintSmall q = 0;
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
			    std::cout << dr << ',' << b << ',' << a << ',' << c << ','
				      << theta << ',' << rho << '\n';
			    std::cout << "z1=(" << xt1 << ' ' << yt1 << ' ' << 0;
			    std::cout << ", z2=(" << xt2 << ' ' << yt2 << ' ' << 0 << ")\n";
			    
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
			    
			    if(n1 > n2) // return the target set with more winning nodes
				return std::move(target1);
			    else
				return std::move(target2);
			};

    
    /* Test for a given x, xo, and xr */
    rocs::Rn x{4.02228, 5.45505, 3.0572};
    rocs::Rn xo{3.08789, 5.67336, -1.31934};
    rocs::Rn xr{0.949472,-0.138768,1.90665};
    const double ro = 0.8;
    const double rt = 0.35;
    std::vector<long long> avoid = local_region(xo, ro, local._encode); // np-based id's
    std::vector<long long> target = local_target(xr, ro, rt, x, local._encode);
    /* Write the corresponding grid points to files */
    rocs::h5FileHandler wtr("logs_replan.h5", H5F_ACC_TRUNC);
    long long ix;
    std::vector<long long> wa, wt;
    for(auto id:avoid) {
	if(id>-1?local._idmap[id]:0) {
	    ix = local._decode[id]/nNodes;
	    wa.push_back(ix);
	}
    }
    wtr.write_array<long long>(wa, "avoid");
    for(auto id:target) {
    	if(id>-1?local._idmap[id]:0) {
    	    ix = local._decode[id]/nNodes;
    	    wt.push_back(ix);
    	}
    }
    wtr.write_array<long long>(wt, "target");


    /* Local rechability: 
     * - forward propagation eliminating avoid area, and
     * - backward reachability to calculate control. */
    NODE_POST p5;
    long long x_index = x_grid.val_to_id(x);
    p5 = nts_ctrlr[encode3[x_index]];
    size_t xp = x_index*nNodes+q;
    std::cout << "Current node id: " << xp << '\n';
    boost::dynamic_bitset<> u_safe(local._na, true);
    std::cout << "Number of winning graph nodes: " << local._nwin << '\n';
    std::cout << "Number of controls: " << local._na << '\n';
    std::vector< std::pair<size_t, int> > replan =
    	local.solve_local_reachavoid(xp, 0.8*target.size(),
    				     avoid, target, u_safe);
    // int first = local.solve_local_reachavoid(xp, 0.8*target.size(),
    // 					     avoid, target, u_safe);

    std::ofstream logger;
    logger.open("logs_newctlr.txt", std::ios::out);
    if(replan[0].first == xp) {
	for(auto &item : replan)
	    logger << item.first << ',' << item.second << '\n';
    } else
	std::cout << "Replanning failed.\n";
    logger.close();
    

    return 0;
}
