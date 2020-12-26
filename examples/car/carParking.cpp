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
#include "src/hdf5io.h"
// #include "src/matlabio.h"

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

    rocs::CSolver solver(&carParking, 0, rocs::ABSMAX);
    solver.init(rocs::GOAL, glb, gub);
    solver.init_goal_area();
    const double e[]{ceps, ceps, ceps};
    solver.init(rocs::AVOID, &cst_v5<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v6<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v7<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v8<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v1rear<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v1front<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v1curb<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v2rear<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v2front<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v2curb<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v3rear<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v3front<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v3curb<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v4rear<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v4front<rocs::ivec>, e);
    solver.init(rocs::AVOID, &cst_v4curb<rocs::ivec>, e);
    solver.init_avoid_area();

    solver._ctlr.retract();
    solver.print_controller_info();

    /* solve the problem */
    const double ei[] = {0.02, 0.02, 0.02};
    const double er[] = {0.2, 0.2, 0.2};
    const double ermin[] = {0.02, 0.02, 0.02};
    solver.cobuchi(&carParking, ei, er, true, ermin);
    // solver.cobuchi(&carParking, seps, rocs::ABSMAX, 0.2, rocs::ABSMAX, seps, true);
    solver.print_controller_info();

    /* save data to file */
    // rocs::matWriter wtr("data_carParking.mat");
    // wtr.open();
    // wtr.write_real_number(H, "H");
    // wtr.write_real_number(L, "L");
    // wtr.write_real_number(D, "D");
    // wtr.write_real_number(d, "d");
    // wtr.write_problem_setting(carParking, solver);
    // wtr.write_sptree_controller(solver);
    // wtr.close();
    std::string datafile = "controller_carParking.h5";
    rocs::h5FileHandler ctlrWtr(datafile, H5F_ACC_TRUNC);
    ctlrWtr.write_number<double>(H, "H");
    ctlrWtr.write_number<double>(L, "L");
    ctlrWtr.write_number<double>(D, "D");
    ctlrWtr.write_number<double>(d, "d");
    ctlrWtr.write_problem_setting< rocs::DTCntlSys<carde> >(carParking);
    ctlrWtr.write_ivec_array(solver._goal, "G");
    ctlrWtr.write_ivec_array(solver._obs, "xobs");
    ctlrWtr.write_sptree_controller(solver);


    return 0;
}
