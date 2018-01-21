/**
 *  DemoCar2.cpp \cite{Gunther2016}
 *
 *  A second specification: multiple obstacles.
 *
 *  Created by Yinan Li on April 23, 2017.
 *
 *  Hybrid Systems Group, UW.
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
    double xub[] = {10, 10, M_PI};

    /* control space */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* functor of dynamics */
    double tau = 0.3;
    car *ptrCar = new car(UD, ulb, uub, mu, tau);

    /* specifications */
    double glb[] = {9, 0, -M_PI};  // goal area
    double gub[] = {9.5, 0.5, M_PI};
    
    double olb1[] = {8.2, 0, -M_PI};  // obstacles
    double oub1[] = {8.4, 8.5, M_PI};
    
    double olb2[] = {8.4, 8.3, -M_PI};
    double oub2[] = {9.3, 8.5, M_PI};
    
    double olb3[] = {9.3, 7.1, -M_PI};
    double oub3[] = {10.0, 7.3, M_PI};
    
    double olb4[] = {8.4, 5.9, -M_PI};
    double oub4[] = {9.3, 6.1, M_PI};

    double olb5[] = {9.3, 4.7, -M_PI};
    double oub5[] = {10.0, 4.9, M_PI};
    
    double olb6[] = {8.4, 3.5, -M_PI};
    double oub6[] = {9.3, 3.7, M_PI};

    double olb7[] = {9.3, 2.3, -M_PI};
    double oub7[] = {10.0, 2.5, M_PI};

    
    /* create a control problem */
    rocs::CntlProb carReach("vehicle", XD, UD, xlb, xub, ptrCar);
    std::cout << carReach << '\n';
    
    /* create a solver */
    rocs::CSolver *solver = new rocs::CSolver(&carReach);
    solver->init(rocs::GOAL, glb, gub);
    solver->init(rocs::AVOID, olb1, oub1);
    solver->init(rocs::AVOID, olb2, oub2);
    solver->init(rocs::AVOID, olb3, oub3);
    solver->init(rocs::AVOID, olb4, oub4);
    solver->init(rocs::AVOID, olb5, oub5);
    solver->init(rocs::AVOID, olb6, oub6);
    solver->init(rocs::AVOID, olb7, oub7);

    solver->reachability_control(0.2, rocs::ABSMAX);

    solver->print_controller_info();
    std::cout << '\n';
    
    /* save the specification and controller */
    carReach.write2mat_settings("data_car_spec2.mat");
    solver->write2mat_controller("data_car_cbox_spec2.mat");
    solver->serialize_controller("data_car_ctree_spec2.mat");


    delete solver;
    delete ptrCar;

    return 1;
}
