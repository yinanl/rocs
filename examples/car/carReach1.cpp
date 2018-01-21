/**
 *  DemoCar1.cpp (Specification 1)
 *
 *  Reachability control of a 2d car \cite{Gunther2016}
 *  
 *  Created by Yinan Li on Feb. 21, 2017.
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
    // std::vector<std::vector<double>> U { {-1, -1}, {-1, 0}, {-1, 1},
    // 					 {0, -1}, {0, 0}, {0, 1},
    // 					 {1, -1}, {1, 0}, {1, 1} };
    
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};

    /* functor of dynamics */
    double tau = 0.3;
    // car *ptrCar = new car(U, tau);
    car *ptrCar = new car(UD, ulb, uub, mu, tau);

    /* create a control problem */
    rocs::CntlProb carReach("vehicle", XD, UD, xlb, xub, ptrCar);
    std::cout << carReach << '\n';

    /* specifications */
    double glb[] = {9, 0, -M_PI};  // goal area
    double gub[] = {9.5, 0.5, M_PI};
    
    double olb[] = {8, 0.3, -M_PI};  // obstacle
    double oub[] = {8.4, 1.2, M_PI};
    // double olb[] = {8.2, 0, -M_PI};
    // double oub[] = {8.4, 2, M_PI};
    
    /* create a solver */
    // CSolver solver2(&carReach);
    // solver2.init(0.001);
    rocs::CSolver *solver1 = new rocs::CSolver(&carReach);
    solver1->init(rocs::GOAL, glb, gub);
    solver1->init(rocs::AVOID, olb, oub);
    
    solver1->reachability_control(0.2, rocs::ABSMAX);

    solver1->print_controller_info();
    
    /* save controller to file */
    carReach.write2mat_settings("data_car_spec1.mat");
    solver1->write2mat_controller("data_car_cbox_spec1.mat");
    solver1->serialize_controller("data_car_ctree_spec1.mat");

    std::cout << '\n';
    

    delete solver1;

    delete ptrCar;

    return 1;
}
