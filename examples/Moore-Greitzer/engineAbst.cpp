/**
 *  engineAbst.cpp
 *
 *  Abstraction-based reach-stay-avoid control of an engine ode.
 *  
 *  Created by Yinan Li on July 3, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include <string>
#include <cmath>

#include "src/grid.h"
#include "src/definitions.h"
#include <boost/numeric/odeint.hpp>

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
struct mgode {

    static const int n = 3;  // system dimension
    static const int nu = 2;  // control dimension
    
    /* template constructor
     * @param[out] dx
     * @param[in] x
     * @param u
     */
    template<typename S>
    mgode(S *dx, const S *x, rocs::Rn u) {
	dx[0] = cx * (aH+H2*(x[0]-W)*(W2-(x[0]-W)*(x[0]-W)) - x[1]) + u[0];
	dx[1] = cy * (x[0] - x[2]*sqrt(x[1]));
	dx[2] = u[1];
    }
};


int main(int argc, char *argv[])
{
    /* set the state space */
    double xlb[3]{0.45, 0.6, 0.5};
    double xub[3]{0.55, 0.7, 0.8};
    
    /* set the control values */
    double Lmu = 0.01;
    double Lu = 0.05;
    double ulb[2]{-Lu, -Lmu};
    double uub[2]{Lu, Lmu};
    double mu[2]{Lu/5, 2*Lmu/5};

    /* set the sampling time and disturbance */
    double t = 0.1;
    double delta = 0.01;
    /* parameters for computing the flow */
    int kmax = 10;
    double tol = 0.01;
    double alpha = 0.5;
    double beta = 2;
    rocs::params controlparams(kmax, tol, alpha, beta);
    
    /* define the control system */
    rocs::CTCntlSys<mgode> engine("Moore-Greitzer", t, mgode::n, mgode::nu,
				  delta, &controlparams);
    engine.init_workspace(xlb, xub);
    engine.init_inputset(mu, ulb, uub);
    engine.allocate_flows();


    /* Abstraction */
    rocs::abstraction<rocs::CTCntlSys<mgode>> abst(&engine);
    const double eta[]{0.001, 0.001, 0.001};
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    double obs[3][2]{{0.520, 0.526}, {0.658, 0.664}, {0.5, 0.8}};
    auto avoid = [&obs, abst, eta](size_t& id) {
		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, id);
    		     double c1 = eta[0]/2.0+1e-10;
    		     double c2 = eta[1]/2.0+1e-10;
		     double c3 = eta[2]/3.0+1e-10;
		     if ((obs[0][0]-c1) <= x[0] && x[0] <= (obs[0][1]+c1) &&
			 (obs[1][0]-c2) <= x[1] && x[1] <= (obs[1][1]+c2) &&
			 (obs[2][0]-c3) <= x[2] && x[2] <= (obs[2][1]+c3))
			 return true;
		     else
			 return false;
		 };
    abst.assign_labels(avoid, -1);
    std::vector<size_t> obstacles;
    for (size_t i = 0; i < abst._x._nv; ++i) {
    	if (abst._labels[i] < 0)
    	    obstacles.push_back(i);
    }
    /* Robustness margins */
    double e1[] = {0,0};
    double e2[] = {0.001, 0.001};
    tb = clock();
    abst.assign_transitions(e1, e2);
    te = clock();
    std::cout << "Time of computing abstraction: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';
    std::cout << "# of transitions: " << abst._ts._ntrans << '\n';


    /*** Control synthesis ***/
    /* case I */
    // std::stringstream convert(argv[1]);
    // double e;
    // convert >> e;
    double e = 0.003;
    double glb1[]{0.5039-e, 0.6605-e, 0.5};
    double gub1[]{0.5039+e, 0.6605+e, 0.8};
    double olb1[]{0.520, 0.658, 0.5};
    double oub1[]{0.526, 0.664, 0.8};
    // /* case II */
    // // double e = 0.003;
    // double glb2[]{0.4519-e, 0.6513-e, 0.5};
    // double gub2[]{0.4519+e, 0.6513+e, 0.8};
    // double olb2[]{0.497, 0.650, 0.5};
    // double oub2[]{0.503, 0.656, 0.8};
    

    engine.release_flows();
    return 0;
}
