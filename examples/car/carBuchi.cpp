/**
 *  DemoCarBuchi.cpp
 *
 *  Buchi objective of a 2d car \cite{Gunther2016}
 *  
 *  Created by Yinan Li on June 4, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <iostream>
#include "src_itvl/problem.h"
#include "src_itvl/csolver.h"
#include "car.hpp"




int main()
{
    const int XD = 3;
    const int UD = 2;
    
    /* state space */
    double xlb[] = {7.3, 0, -M_PI};
    double xub[] = {10, 2, M_PI};

    /* control space */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* functor of dynamics */
    double tau = 0.3;
    car *ptrCar = new car(UD, ulb, uub, mu, tau);

    /* create a control problem */
    rocs::CntlProb carReach("vehicle", XD, UD, xlb, xub, ptrCar);
    std::cout << carReach << '\n';

    /* specifications */
    double glb[] = {9, 0, -M_PI};  // goal area
    double gub[] = {9.5, 0.5, M_PI};
    
    double olb[] = {8, 0.3, -M_PI};  // obstacle
    double oub[] = {8.4, 1.2, M_PI};
    
    /* create a solver */
    rocs::CSolver *solver = new rocs::CSolver(&carReach);
    solver->init(rocs::GOAL, glb, gub);
    solver->init(rocs::AVOID, olb, oub);
    
    solver->buchi(rocs::ABSMAX, 0.2);

    solver->print_controller_info();
    
    /* save controller to file */
    carReach.write2mat_settings("data_car_spec3.mat");
    solver->write2mat_controller("data_car_cbox_spec3.mat");
    solver->serialize_controller("data_car_ctree_spec3.mat");

    std::cout << '\n';
    

    delete solver;

    delete ptrCar;

    return 1;
}
