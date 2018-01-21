/**
 *  dcdc.cpp
 *
 *  DCDC converter example (benchmark)
 *
 *  Created by yinan li on Feb. 23, 2017.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>
#include "src_itvl/problem.h"
#include "src_itvl/csolver.h"
#include "dcdc.hpp"


int main()
{
    const int XD = 2;
    const int UD = 1;
    
    /* state space */
    // double xlb[] = {-2, -1.5};
    // double xub[] = {2, 3};
    double xlb[] = {-2, 0.70};
    double xub[] = {2, 1.50};

    rocs::input_type U {{1}, {2}};  // two modes

    /* invariant area */
    double glb[] = {1.15, 1.09};
    double gub[] = {1.55, 1.17};

    /* functor of dynamics */
    DCDC *ptrDC = new DCDC(U, TS);


    /* define an invariance control problem */
    rocs::CntlProb dcdcInv("dcdc", XD, UD, xlb, xub, ptrDC);
    std::cout << dcdcInv << '\n';

    /* use solver to design a controller */
    rocs::CSolver *solver = new rocs::CSolver(&dcdcInv);

    solver->init(rocs::GOAL, glb, gub);
    
    solver->invariance_control(0.001, rocs::RELMAXG);

    solver->print_controller_info();


    /* save controller to file */
    dcdcInv.write2mat_settings("data_dcdc_spec.mat");
    solver->write2mat_controller("data_dcdc_cbox.mat");
    solver->serialize_controller("data_dcdc_ctree.mat");


    delete solver;
    
    delete ptrDC;

    return 1;
}
