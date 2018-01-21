/**
 *  temp.cpp
 *
 *  Temperature control example
 *
 *  Created by yinan li on August 14, 2017.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>

#include "src_itvl/csolver.h"
#include "temp.hpp"




int main()
{
    const double XD = 2;
    const double UD = 1;
    double tau = 10;  // 50 secs
    
    /* state space */
    /* x[0]= room temperature, x[1]=heater temperature */
    double xlb[] = {10, 12};
    double xub[] = {28, 30};

    rocs::input_type U {{1}, {2}, {3}, {4}};  // four modes

    /* target area */
    double glb[] = {18, 20};
    double gub[] = {20, 22};

    /* avoid area */
    double olb[] = {22, 26};
    double oub[] = {28, 30};

    /* functor of dynamics */
    TPC *ptrTPC = new TPC(U, tau);


    /* define an invariance control problem */
    rocs::CntlProb tpcRI("tpc", XD, UD, xlb, xub, ptrTPC);
    std::cout << tpcRI << '\n';

    /* use solver to design a controller */
    rocs::CSolver *solver = new rocs::CSolver(&tpcRI);

    solver->init(rocs::GOAL, glb, gub);
    solver->init(rocs::AVOID, olb, oub);
    
    // solver->reach_stay(0.1, ABSMAX, 0.1, ABSMAX);
    // solver->invariance_control(0.5, ABSMAX);
    solver->cobuchi(0.15, rocs::ABSMAX, 0.15, rocs::ABSMAX);

    solver->print_controller_info();
    // solver->print_controller();


    /* save controller to file */
    tpcRI.write2mat_settings("data_tpc_spec.mat");
    solver->write2mat_controller("data_tpc_cbox.mat");
    solver->serialize_controller("data_tpc_ctree.mat");


    delete solver;
    
    delete ptrTPC;

    return 1;
}
