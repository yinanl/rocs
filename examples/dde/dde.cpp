/**
 *  dcdc.cpp
 *
 *  A DCDC converter invariance control example.
 *
 *  Created by Yinan Li on May 10, 2018.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"




/* Discrete-time dynamics of DCDC converter */
struct dde {
    static const int n = 2;  // state dimension
    
    /**
     * Constructors: 
     * real-valued (arma::vec) and interval-valued (rocs::ivec)
     * @param[out] y the next state after the sampling time.
     * @param[in] x the current state.
     */
    template<typename S>
    dde(S &y, const S &x) {
	y[0] = x[0] - 0.5*x[1];
	y[1] = x[0];
    }
    
};

struct dde4 {
    static const int n = 4;  // state dimension
    
    /**
     * Constructors: 
     * real-valued (arma::vec) and interval-valued (rocs::ivec)
     * @param[out] y the next state after the sampling time.
     * @param[in] x the current state.
     */
    template<typename S>
    dde4(S &y, const S &x) {
	y[0] = 0.0809*x[0] - 0.0588*x[1] + 0.8257*x[2] - 0.1308*x[3];
	y[1] = 0.0588*x[0] + 0.0809*x[1] + 0.1308*x[2] + 0.8257*x[3];
	y[2] = x[0];
	y[3] = x[1];
    }
    // dde4(S &y, const S &x) {
    // 	y[0] = 0.1*x[0] + 0.1*x[2] - 0.2*x[3];
    // 	y[1] = 0.4*x[0] + 0.1*x[1] + 0.4*x[2] +0.5*x[3];
    // 	y[2] = x[0];
    // 	y[3] = x[1];
    // }
    
};

template<typename T>
T cset(const T &x) {
    const double c = 1.0;
    T y(8);
    y[0] = x[0] + x[1] - c;
    y[1] = x[0] - x[1] -c;
    y[2] = -x[0] + x[1] -c;
    y[3] = -x[0] - x[1] -c;
    y[4] = x[2] + x[3] - c;
    y[5] = x[2] - x[3] -c;
    y[6] = -x[2] + x[3] -c;
    y[7] = -x[2] - x[3] -c;

    return y;
}




int main()
{
    // /* case 1 */
    // /* set the state space */
    // double xlb[] = {-2, -2};
    // double xub[] = {2, 2};

    // /* define the control system */
    // rocs::DTSys<dde> ddeInv("dcdc", 1, dde::n);
    // ddeInv.init_workspace(xlb, xub);
    
    // /* set the specifications */
    // double glb[] = {-1, -1};
    // double gub[] = {1, 1};

    // /* solve the problem */
    // rocs::CSolver solver(&ddeInv, rocs::RELMAXG);
    // double eps[]{0.001, 0.001};
    // solver.init(rocs::GOAL, glb, gub);
    // solver.invariance_control(&ddeInv, eps);
    // solver.print_controller_info();


    /* case 2 */
    double xlb[] = {-2, -2, -2, -2};
    double xub[] = {2, 2, 2, 2};
    rocs::DTSys<dde4> ddeInv("dcdc", 1, dde4::n);
    ddeInv.init_workspace(xlb, xub);
    
    rocs::CSolver solver(&ddeInv, 0, rocs::ABSMAX);
    double eps[]{0.1, 0.1, 0.1, 0.1};
    solver.init(rocs::GOAL, &cset<rocs::ivec>, eps);
    solver.invariance_control(&ddeInv, eps);
    solver.print_controller_info();
    std::cout << solver._ctlr._root->_box.getdim() <<'\n';


    /* save the problem data and the solution */
    rocs::matWriter wtr("data_ddeInv.mat");
    wtr.open();
    wtr.write_problem_setting(ddeInv, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();

    return 0;
}
