/**
 *  carAbst.cpp
 *
 *  The vehicle example from scots using abstraction-based control synthesis.
 *
 *  Created by Yinan Li on Feb. 18, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "src/abstraction.hpp"
#include "src/dsolver.h"
#include "src/matlabio.h"
#include "src/txtfileio.h"

#include "car.hpp"


int main(int argc, char *argv[])
{
    clock_t tb, te;

    /* set the state space */
    const double theta = 3.4;
    double xlb[] = {0, 0, -theta};
    double xub[] = {10, 10, theta};

    /* set the control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* Robustness margins */
    double e1[] = {0,0};
    double e2[] = {0.001, 0.001};

    /* define the control system */
    rocs::DTCntlSys<carde> car("reach goal", h, carde::n, carde::m);
    car.init_workspace(xlb, xub);
    car.init_inputset(mu, ulb, uub);

    /* Abstraction */
    const double eta[] = {0.2, 0.2, 0.2}; /* set precision */
    rocs::abstraction<rocs::DTCntlSys<carde>> abst(&car);
    abst.init_state(eta, xlb, xub);
    std::cout << "The number of abstraction states: " << abst._x._nv << '\n';
    /* Set avoid area by label -1 */
    double obs[15][4] = {
    			 { 1  , 1.2, 0  ,   9 },
    			 { 2.2, 2.4, 0  ,   5 },
    			 { 2.2, 2.4, 6  ,  10 },
    			 { 3.4, 3.6, 0  ,   9 },
    			 { 4.6, 4.8, 1  ,  10 },
    			 { 5.8, 6  , 0  ,   6 },
    			 { 5.8, 6  , 7  ,  10 },
    			 { 7  , 7.2, 1  ,  10 },
    			 { 8.2, 8.4, 0  ,  8.5},
    			 { 8.4, 9.3, 8.3,  8.5},
    			 { 9.3, 10 , 7.1,  7.3},
    			 { 8.4, 9.3, 5.9,  6.1},
    			 { 9.3, 10 , 4.7,  4.9},
    			 { 8.4, 9.3, 3.5,  3.7},
    			 { 9.3, 10 , 2.3,  2.5}
    };
    auto avoid = [&obs,abst,eta](size_t& id) {
    		     std::vector<double> x(abst._x._dim);
    		     abst._x.id_to_val(x, id);
    		     double c1= eta[0]/2.0+1e-10;
    		     double c2= eta[1]/2.0+1e-10;
    		     for(size_t i=0; i<15; i++) {
    			 if ((obs[i][0]-c1) <= x[0] && x[0] <= (obs[i][1]+c1) &&
    			     (obs[i][2]-c2) <= x[1] && x[1] <= (obs[i][3]+c2))
    			     return true;
    		     }
    		     return false;
    		 };

    abst.assign_labels(avoid, -1);
    std::vector<size_t> obstacles;
    for (size_t i = 0; i < abst._x._nv; ++i) {
    	if (abst._labels[i] < 0)
    	    obstacles.push_back(i);
    }
    tb = clock();
    abst.assign_transitions(e1, e2);
    te = clock();
    std::cout << "Time of computing abstraction: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';
    std::cout << "# of transitions: " << abst._ts._ntrans << '\n';

    // /* Write post and pre transitions of grid points (within a given area)
    //  * to txt files. */
    // double lb[] = {9.2, 4.0, -3.2};
    // double ub[] = {9.6, 4.4, 3.2};
    // std::cout << "Write post transitions to trans_post.txt \n";
    // rocs::txtWriter twtr("trans_post.txt");
    // if (twtr.open()) {
    // 	twtr.write_post_transitions(abst._ts, abst._x, lb, ub);
    // 	twtr.close();
    // 	std::cout << "Done.\n";
    // }

    // std::cout << "Write pre transitions to trans_pre.txt \n";
    // if (twtr.open("trans_pre.txt", std::ios::out)) {
    // 	twtr.write_pre_transitions(abst._ts, abst._x, lb, ub);
    // 	twtr.close();
    // 	std::cout << "Done.\n";
    // }


    /**
     * Reachability control
     */
    double glb[] = {9, 0, -theta};  // goal area
    double gub[] = {9.5, 0.5, theta};
    std::vector<size_t> target = abst.get_discrete_states(glb, gub);
    rocs::DSolver solver(&(abst._ts));
    tb = clock();
    solver.reachability(target);
    te = clock();
    std::cout << "Time of solving reachability control: "
    	      << (float)(te - tb)/CLOCKS_PER_SEC << '\n';


    /**
     * Display and save memoryless controllers.
     */
    rocs::matWriter wtr("data_abstbased.mat");
    wtr.open();
    wtr.write_uniform_grids(abst._x, "xgrid");
    wtr.write_real_array(obstacles, "avoidset");
    wtr.write_real_array(target, "goalset");
    // wtr.write_transitions(abst._ts, "trans_post", "post", "postptr", "trans_pre", "pre", "preptr");
    wtr.write_discrete_controller(solver, "leastctlr", "optctlr");
    wtr.close();


    return 0;
}
