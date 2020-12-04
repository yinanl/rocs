/**
 *  testAbst.cpp
 *
 *  Test abstraction by using a car kinematics model.
 *
 *  Created by Yinan Li on Nov. 14, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>

#include "src/grid.h"
#include "src/definitions.h"
#include "src/abstraction.hpp"
#include "src/hdf5io.h"


/* user defined dynamics */
struct car_ode {
    rocs::Rn u;
    car_ode (const rocs::Rn param): u (param) {}
    /**
     * ODE model
     * @param x system state: [x,y,theta], n=3
     * @param dxdt vector field
     * @param t time
     */
    void operator() (rocs::Rn &x, rocs::Rn &dxdt, double t) const
    {
	dxdt[0] = u[0]*std::cos(x[2]);
	dxdt[1] = u[0]*std::sin(x[2]);
	dxdt[2] = u[1];
    }
};

const double h = 0.3;  // sampling time
const double dt = 0.001; //integration step size for odeint

struct carde { // discrete-time model (difference equation)
    static const int n = 3;  // system dimension
    static const int m = 2;
    /**
     * Discrete-time dynamics
     * @param h sampling time
     * @param x system state: [x,y,theta], n=3
     * @param u control array (size of 2, velocity and steering angle)
     * @param nu the number of different control values
     */
    template<typename S>
    carde(S &dx, const S &x, rocs::Rn u) {
	if (std::fabs(u[0]) < 1e-6) { //v=0
	    dx[0] = x[0];
	    dx[1] = x[1];
	    dx[2] = x[2] + u[1] * h;
	} else if (std::fabs(u[1]) < 1e-6) { //w=0
	    dx[0] = x[0] + u[0]* cos(x[2])*h;
	    dx[1] = x[1] + u[0]* sin(x[2])*h;
	    dx[2] = x[2];
	} else { //v,w not 0
	    dx[0] = x[0] + u[0]/u[1]*2*sin(u[1]*h/2.)*cos(x[2]+u[1]*h/2.);
	    dx[1] = x[1] + u[0]/u[1]*2*sin(u[1]*h/2.)*sin(x[2]+u[1]*h/2.);
	    dx[2] = x[2] + u[1] * h;
	}
    }
    
}; // struct carde


struct twoagent {
    static const int n = 3;  // system dimension
    static const int nu = 2;  // control dimension
    rocs::ivec d{rocs::interval(-0.8, 0.8),
		 rocs::interval(-0.8, 0.8)};

    /* template constructor
     * @param[out] dx
     * @param[in] x = [xr, yr, psir]
     * @param u = [v, w]
     * @param d = [v', w']
     */
    template<typename S>
    twoagent(S *dx, const S *x, rocs::Rn u) {
	dx[0] = -u[0] + d[0]*cos(x[2]) + u[1]*x[1];
	dx[1] = d[0]*sin(x[2]) - u[1]*x[0];
	dx[2] = d[1] - u[1];
    }
};


int main()
{
    /* Config */
    clock_t tb, te;
    boost::numeric::odeint::runge_kutta_cash_karp54<rocs::Rn> rk45;

    /**
     * Case I
     */
    /* Set the state and control space */
    const int xdim = 3;
    const int udim = 2;
    
    double xlb[] = {-3, -3, -M_PI};
    double xub[] = {3, 3, M_PI};
    double eta[] = {0.2, 0.2, 0.2};
    
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /**
     * Define the two-agent system
     */
    double t = 0.3;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 5;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    rocs::CTCntlSys<twoagent> safety("collision-free", t,
    				     twoagent::n, twoagent::nu,
    				     delta, &controlparams);

    safety.init_workspace(xlb, xub);
    safety.init_inputset(mu, ulb, uub);
    safety.allocate_flows();

    rocs::abstraction< rocs::CTCntlSys<twoagent> > abst(&safety);
    abst.init_state(eta, xlb, xub);
    std::cout << "# of in-domain nodes: " << abst._x._nv << '\n';
    /**
     * Assign 1 to the target invariant set and 0 to others.
     * Mark 0 for any box intersect or inside the cylinder: x^2+y^2<=rmin^2, any phi.
     * The invariant set is the region outside of the cylinder.
     */
    auto inv_set = [&abst, &eta](size_t i) {
    		       const double rmin = 1.21;
    		       std::vector<double> x(abst._x._dim);
    		       abst._x.id_to_val(x, i);
    		       double xl = x[0] - eta[0]/2.;
    		       double xr = x[0] + eta[0]/2.;
    		       double yl = x[1] - eta[1]/2.;
    		       double yr = x[1] + eta[1]/2.;
    		       double xsqr = (xr*xr) > (xl*xl) ? (xl*xl) : (xr*xr);
    		       double ysqr = (yr*yr) > (yl*yl) ? (yl*yl) : (yr*yr);
    		       if(xsqr + ysqr < rmin*rmin)
    			   return 0;
    		       else
    			   return 1;
    		   };
    abst.assign_labels(inv_set);
    abst.assign_label_outofdomain(1); //out of domain is safe
    
    std::string transfile = "abstca_0.2-0.2-0.2.h5";
    struct stat buffer;
    if(stat(transfile.c_str(), &buffer) == 0) {
    	/* Read from a file */
    	std::cout << "Reading transitions...\n";
    	rocs::h5FileHandler transRdr(transfile, H5F_ACC_RDONLY);
    	tb = clock();
    	transRdr.read_transitions(abst._ts);
    	te = clock();
    } else {
    	std::cout << "No transition file found. Computing transitions...\n";
    	/* Robustness margins */
    	double e1[] = {0,0,0};
    	double e2[] = {0,0,0};
    	tb = clock();
    	abst.assign_transitions(e1, e2);
    	te = clock();
    	/* Write transitions to file */
    	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
    	transWtr.write_transitions(abst._ts);
    }
    float time = (float)(te - tb)/CLOCKS_PER_SEC;
    std::cout << "Time of reading/computing abstraction: " << time << '\n';
    std::cout << "# of all nodes: " << abst._ts._nx << '\n';
    std::cout << "# of actions: " << abst._ts._nu << '\n';
    std::cout << "# of transitions: " << abst._ts._ntrans << '\n';


    /**
     * Case II
     */
    // /* Set the state space */
    // const int xdim = 3;
    // const int udim = 2;
    // const double theta = 3.5;
    // double xlb[] = {0, 0, -theta};
    // double xub[] = {10, 10, theta};
    // double eta[] = {0.2, 0.2, 0.2};
    // /* Set the control values */
    // double ulb[] = {-1.0, -1.0};
    // double uub[] = {1.0, 1.0};
    // double mu[] = {0.3, 0.3};
    // /* Define the control system */
    // rocs::DTCntlSys<carde> car("DBA", h, carde::n, carde::m);
    // car.init_workspace(xlb, xub);
    // car.init_inputset(mu, ulb, uub);
    
    // rocs::abstraction< rocs::DTCntlSys<carde> > abst(&car);
    // abst.init_state(eta, xlb, xub);
    // std::cout << "# of in-domain nodes: " << abst._x._nv << '\n';
    // /* Assign the label of avoid area to -1 */
    // rocs::UintSmall nAvoid = 4;
    // double obs[4][4] = {
    // 	{1.6, 5.7, 4.0, 5.0},
    // 	{3.0, 5.0, 5.0, 8.0},
    // 	{4.3, 5.7, 1.8, 4.0},
    // 	{5.7, 8.5, 1.8, 2.5}
    // };
    // auto label_avoid = [&obs, &nAvoid, &abst, &eta](size_t i) {
    // 		     std::vector<double> x(abst._x._dim);
    // 		     abst._x.id_to_val(x, i);
    // 		     double c1= eta[0]/2.0+1e-10;
    // 		     double c2= eta[1]/2.0+1e-10;
    // 		     for(size_t i = 0; i < nAvoid; ++i) {
    // 			 if ((obs[i][0]-c1) <= x[0] && x[0] <= (obs[i][1]+c1) &&
    // 			     (obs[i][2]-c2) <= x[1] && x[1] <= (obs[i][3]+c2))
    // 			     return -1;
    // 		     }
    // 		     return 0;
    // 		 };
    // abst.assign_labels(label_avoid);
    // abst.assign_label_outofdomain(-1);
    // std::vector<size_t> obstacles;
    // for (size_t i = 0; i < abst._x._nv; ++i) {
    // 	if (abst._labels[i] < 0)
    // 	    obstacles.push_back(i);
    // }

    // /* Compute/Read abstraction */
    // float tabst;
    // std::string transfile = "abstfull_0.2-0.2-0.2.h5";
    // struct stat buffer;
    // if(stat(transfile.c_str(), &buffer) == 0) {
    // 	/* Read from a file */
    // 	std::cout << "Reading transitions...\n";
    // 	rocs::h5FileHandler transRdr(transfile, H5F_ACC_RDONLY);
    // 	tb = clock();
    // 	transRdr.read_transitions(abst._ts);
    // 	te = clock();
    // } else {
    // 	std::cout << "No transition file found. Computing transitions...\n";
    // 	/* Robustness margins */
    // 	double e1[] = {0,0,0};
    // 	double e2[] = {0,0,0};
    // 	tb = clock();
    // 	abst.assign_transitions(e1, e2);
    // 	te = clock();
    	
    // 	/* Write abstraction to file */
    // 	rocs::h5FileHandler transWtr(transfile, H5F_ACC_TRUNC);
    // 	transWtr.write_transitions(abst._ts);
    // }
    // tabst = (float)(te - tb)/CLOCKS_PER_SEC;
    // std::cout << "Time of reading/computing abstraction: " << tabst << '\n';
    // std::cout << "# of all nodes: " << abst._ts._nx << '\n';
    // std::cout << "# of actions: " << abst._ts._nu << '\n';
    // std::cout << "# of transitions: " << abst._ts._ntrans << '\n';


    /* Test */
    size_t na = abst._ts._nu;
    size_t nx = abst._ts._nx;
    size_t si, sk, k;
    bool suc = 0;
    int np;    
    
    
    /* Test post-pre consistency */
    std::cout << "Checking post-pre consistency...\n";
    for(size_t i = 0; i < nx; ++i) {
	for(size_t j = 0; j < na; ++j) {
	    si = abst._ts._ptrpost[i*na+j];
	    for(size_t p=si; p<si+abst._ts._npost[i*na+j]; ++p) {
		k = abst._ts._idpost[p];
		/* Test if the pre of post by j contains i */
		suc = 0;
		sk = abst._ts._ptrpre[k*na+j];
		// /********** logging **********/
		// if(i == 0 && j == 16 && k == 0) {
		//     std::cout << "The predecessors of " << k << " with " << j << ": ";
		// }
		// /********** logging **********/
		for(size_t pp=sk; pp<sk+abst._ts._npre[k*na+j]; ++pp) {
		    // /********** logging **********/
		    // if(i == 0 && j == 16 && k == 0) {
		    // 	std::cout << "idpre["<< pp << "]="
		    // 		  << abst._ts._idpre[pp] << '\n';
		    // }
		    // /********** logging **********/
		    if(abst._ts._idpre[pp] == i) {
			suc = 1;
			break;
		    }
		}
		if(i == 0 && j == 16 && k == 0) {
		    std::cout << '\n';
		}
		if(!suc) {//two cases: npre(k,j)=0 or no i in npre(k, j)
		    std::cout << "Post and pre transitions are inconsistent "
			      << i << "->(" << j << ")->" << k << '\n';
		    return -1;
		}
		    
	    }
	}
    }
    std::cout << "Every post transition has its corresponding pre transition.\n";

    for(size_t i = 0; i < nx; ++i) {
	for(size_t j = 0; j < na; ++j) {
	    si = abst._ts._ptrpre[i*na+j];
	    for(size_t p=si; p<si+abst._ts._npre[i*na+j]; ++p) {
		k = abst._ts._idpre[p];
		/* Test if the post of pre by j contains i */
		suc = 0;
		sk = abst._ts._ptrpost[k*na+j];
		for(size_t pp=sk; pp<sk+abst._ts._npost[k*na+j]; ++pp) {
		    if(abst._ts._idpost[pp] == i) {
			suc = 1;
			break;
		    }
		}
		if(!suc) {//two cases: npost(k,j)=0 or no i in npost(k, j)
		    std::cout << "Post and pre transitions are inconsistent "
			      << k << "->(" << j << ")->" << i << '\n';
		    return -1;
		}
		    
	    }
	}
    }
    std::cout << "Every pre transition has its corresponding post transition.\n";


    // /* Test reachable set computation */
    // std::cout << "Checking post transitions by rechable set computation...\n";
    // rocs::Rn x(xdim);
    // rocs::Rn u(udim);
    // rocs::Rn xpost(xdim);
    // rocs::ivec box(xdim);
    // std::vector<rocs::ivec> reachset(na, rocs::ivec(xdim));
    // // std::vector<rocs::Rn> corners(std::pow(2, xdim), rocs::Rn(xdim));
    // rocs::Rn corner(xdim);
    // int quo, rem;
    // rocs::ivec margin{rocs::interval(-rocs::EPSIVAL, rocs::EPSIVAL),
    // 		      rocs::interval(-rocs::EPSIVAL, rocs::EPSIVAL),
    // 		      rocs::interval(-rocs::EPSIVAL, rocs::EPSIVAL)};
    // rocs::ivec yt(xdim);
    // for(size_t i = 0; i < nx; ++i) {
    // 	// std::cout << "State x= " << '(' << x[0] << ',' << x[1] << ',' << x[2] << "):\n";
    // 	if(i < abst._x._nv) { //belongs to xgrid
    // 	    /* Compute the reachable set */
    // 	    abst._x.id_to_val(x, i); //x is the center of the box i
    // 	    for(int d = 0; d < xdim; ++d)
    // 		box.setval(d, rocs::interval(x[d]-eta[d]/2., x[d]+eta[d]/2.));
    // 	    car.get_reach_set(reachset, box);
	    
    // 	    /* Test valid control inputs */
    // 	    for(size_t j = 0; j < na; ++j) {
    // 		if(abst._ts._npost[i*na+j] > 0) {
    // 		    car._ugrid.id_to_val(u, j); //get control values
    // 		    /* Test if the reachable set covers ode solutions of all corners */
    // 		    for(int k = 0; k < std::pow(2, xdim); ++k) {
    // 			quo = k;
    // 			for(int d = 0; d < xdim; ++d) {
    // 			    if(quo % 2) {
    // 				corner[d] = x[d]+eta[d]/2.; //upper bound
    // 			    } else {
    // 				corner[d] = x[d]-eta[d]/2.; //lower bound
    // 			    }
    // 			    quo /= 2;
    // 			}
    // 			// std::cout << "Corner "
    // 			// 	  << '(' << corner[0] << ',' << corner[1] << ',' << corner[2] << ")\n";
    // 			boost::numeric::odeint::integrate_const(rk45, car_ode(u), corner, 0.0, h, dt);
    // 			yt = reachset[j] + margin;
    // 			if(!yt.isin(corner)) {
    // 			    std::cout << "The reachable set is incorrect with u="
    // 				      << '(' << u[0] << ',' << u[1] << "):"
    // 				      << '(' << corner[0] << ',' << corner[1] << ',' << corner[2] << ')'
    // 				      << " is not in " << yt << '\n'
    // 				      << "Test terminates.\n";
    // 			    return -1;
    // 			}
    // 		    }
    // 		    /* Test if all post nodes are in the reachable set (soundness) */
    // 		    si = abst._ts._ptrpost[i*na+j];
    // 		    for(size_t p = si; p<si+abst._ts._npost[i*na+j]; ++p) {
    // 			abst._x.id_to_val(xpost, abst._ts._idpost[p]); //xpost: post node center
    // 			for(int d = 0; d < xdim; ++d) //box: post interval centered at xpost
    // 			    box.setval(d, rocs::interval(xpost[d]-eta[d]/2., xpost[d]+eta[d]/2.));
    // 			if(reachset[j].isout(box)) { //box and reachset[j] should intersect
    // 			    std::cout << "Post transition for xid,uid=" << i << ',' << j
    // 				      << " is incorrect.\n"
    // 				      << "Test terminates.\n";
    // 			    return -1;
    // 			}
    // 		    }
    // 		}
    // 	    }//end for control values
    // 	} else { //out-of-domain node
    // 	    std::cout << "Checking the out-of-domain node...\n";
    // 	}
    // }

    return 0;
}
