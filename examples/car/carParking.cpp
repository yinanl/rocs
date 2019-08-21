/**
 *  carParking.cpp
 *
 *  Automatic parallel parking of a vehicle.
 *  
 *  Created by Yinan Li on May 12, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>

#include "src/system.hpp"
#include "src/csolver.h"

#include "src/matlabio.h"

#include "car.hpp"


/*** constraints ***/
const double d = 0.5; // marginal parking space

const double H = 1.0; // car width
const double L = 2.0; // car length
const double D = 0.5; // distance to the curb

/* rear and front car corners collide with the car body */
template<typename T>
T cst_v5(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}
template<typename T>
T cst_v6(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]-L) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]-L) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]-L) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]-L) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}
template<typename T>
T cst_v7(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]-(2*L+d)) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]-(2*L+d)) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]-(2*L+d)) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]-(2*L+d)) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}
template<typename T>
T cst_v8(const T &x) {
    T y(4);
    y[0] = sin(x[2])*(x[0]-(3*L+d)) - cos(x[2])*(x[1]-H) - H/2.0;
    y[1] = cos(x[2])*(x[0]-(3*L+d)) + sin(x[2])*(x[1]-H) - L/2.0;
    y[2] = -sin(x[2])*(x[0]-(3*L+d)) + cos(x[2])*(x[1]-H) - H/2.0;
    y[3] = -cos(x[2])*(x[0]-(3*L+d)) - sin(x[2])*(x[1]-H) - L/2.0;
    return y;
}

/* the car body collides with the rear car area: consider 4 corners */
template<typename T>
T cst_v1rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(H/2.0) - H;
    return y;
}
template<typename T>
T cst_v2rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(-H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(-H/2.0) - H;
    return y;
}
template<typename T>
T cst_v3rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(-H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(-H/2.0) - H;
    return y;
}
template<typename T>
T cst_v4rear(const T &x) {
    T y(2);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(H/2.0) - L;
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(H/2.0) - H;
    return y;
}

/* the car body collides with the front car area: consider 4 corners */
template<typename T>
T cst_v1front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(H/2.0) - (3.0*L+d);
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(L/2.0) + sin(x[2])*(H/2.0) + (2.0*L+d);
    return y;
}
template<typename T>
T cst_v2front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(L/2.0) - sin(x[2])*(-H/2.0) - (3*L+d);
    y[1] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(-H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(L/2.0) + sin(x[2])*(-H/2.0) + (2.0*L+d);
    return y;
}
template<typename T>
T cst_v3front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(-H/2.0) - (3*L+d);
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(-H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(-L/2.0) + sin(x[2])*(-H/2.0) + (2.0*L+d);
    return y;
}
template<typename T>
T cst_v4front(const T &x) {
    T y(3);
    y[0] = x[0] + cos(x[2])*(-L/2.0) - sin(x[2])*(H/2.0) - (3*L+d);
    y[1] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(H/2.0) - H;
    y[2] = -x[0] - cos(x[2])*(-L/2.0) + sin(x[2])*(H/2.0) + (2.0*L+d);
    return y;
}

/* the car body collides with the curb area */
template<typename T>
T cst_v1curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(H/2.0) + D;
    return y;
}
template<typename T>
T cst_v2curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(L/2.0) + cos(x[2])*(-H/2.0) + D;
    return y;
}
template<typename T>
T cst_v3curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(-H/2.0) + D;
    return y;
}
template<typename T>
T cst_v4curb(const T &x) {
    T y(1);
    y[0] = x[1] + sin(x[2])*(-L/2.0) + cos(x[2])*(H/2.0) + D;
    return y;
}


int main() {

    /* set the state space */
    double alb = -M_PI/2.5; double aub = M_PI/2.5;
    
    double xlb[] = {0, 0, alb};
    double xub[] = {8, 4, aub};

    /* set the control values */
    double ulb[] = {-1.0, -1.0};
    double uub[] = {1.0, 1.0};
    double mu[] = {0.3, 0.3};
    
    /* define the control system */
    rocs::DTCntlSys<carde> carParking("parallel parking", h, carde::n, carde::m);
    carParking.init_workspace(xlb, xub);
    carParking.init_inputset(mu, ulb, uub);

    /* load problem settings */
    double ceps = 0.02; // precision for collision area over-approximation
    double glb[] = {L+L/2.0+ceps, 0.5, -M_PI*3.0/180.0};
    double gub[] = {L+L/2.0+d-2*ceps, 0.6, M_PI*3.0/180.0};
    
    rocs::CSolver solver(&carParking);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    
    solver.init(rocs::AVOID, &cst_v5<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v6<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v7<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v8<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v1rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v1front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v1curb<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v2rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v2front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v2curb<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v3rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v3front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v3curb<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v4rear<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v4front<rocs::ivec>, ceps);
    solver.init(rocs::AVOID, &cst_v4curb<rocs::ivec>, ceps);
    solver.init_avoid_area();
    
    solver._ctlr.retract();
    solver.print_controller_info();

    // /* setting 1: wide space */
    // double glb[] = {3.25, 0.5, -M_PI*3.0/180.0};  // goal area
    // double gub[] = {3.75, 0.6, M_PI*3.0/180.0};
    
    // double olb1[] = {0.0, 0.0, alb};  // rear car: [0,2]x[0,1]
    // double oub1[] = {3.0, 1.7, aub};
    // double olb2[] = {5.0, 0.0, alb};  // front car: [6,8]x[0,1]
    // double oub2[] = {8.0, 1.7, aub};

    // /* setting 2: narrow space */
    // double glb[] = {3.02, 0.5, -M_PI*3.0/180.0};  // goal area
    // double gub[] = {3.22, 0.6, M_PI*3.0/180.0};
    
    // double olb1[] = {0.0, 0.0, alb};  // rear car: [0,2]x[0,1]
    // double oub1[] = {3.0, 1.7, aub};
    // double olb2[] = {3.3, 0.0, alb};  // front car: [4.3,6.3]x[0,1]
    // double oub2[] = {6.3, 1.7, aub};

    // rocs::CSolver solver(&carParking);
    // solver.init(rocs::GOAL, glb, gub);
    // solver.init(rocs::AVOID, olb1, oub1);
    // Solver.init(rocs::AVOID, olb2, oub2);

    // solver.cobuchi(&carParking, 0.1, rocs::RELMAXG, 0.2, rocs::RELMAXW, 0.04, true);
    // solver.print_controller_info();


    /* solve the problem */
    double seps = 0.02;
    solver.cobuchi(&carParking, seps, rocs::ABSMAX, 0.2, rocs::ABSMAX, seps, true);
    solver.print_controller_info();

    /* save data to file */
    rocs::matWriter wtr("data_carParking.mat");
    wtr.open();
    wtr.write_real_number(H, "H");
    wtr.write_real_number(L, "L");
    wtr.write_real_number(D, "D");
    wtr.write_real_number(d, "d");
    wtr.write_problem_setting(carParking, solver);
    wtr.write_sptree_controller(solver);
    wtr.close();
    
    
    return 0;
}
