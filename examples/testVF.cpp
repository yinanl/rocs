/*
 *  testVF.cpp
 *
 *  test interval-valued vector field
 *
 *  Created by yinan li on Apl. 06, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ExampleDynamics
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>

#include <cmath>
#include <iostream>
#include "src_itvl/vectorfield.h"
#include "examples/dcdc/dcdc.hpp"
#include "examples/car/car.hpp"
#include "examples/ipdl/ipdl.hpp"
#include "examples/temp/temp.hpp"


typedef rocs::ivec (*fcst)(const rocs::ivec&);

/**
 * test case 1: Boost DCDC converter
 */
BOOST_AUTO_TEST_CASE(test_dcdc)
{
    rocs::ivec x(2);
    x[0] = rocs::interval(1.2, 1.25);
    x[1] = rocs::interval(1.11, 1.12);

    rocs::input_type U {{1}, {2}};  // two modes

    DCDC dcdc = DCDC(U, TS);
    std::vector<rocs::ivec> y = dcdc(x);

    std::cout << y[0] << '\n';
    std::cout << y[1] << '\n';

    BOOST_CHECK(true);
    
}


/*
 * test case 2: car
 */
BOOST_AUTO_TEST_CASE(test_car)
{
    rocs::ivec x(3);
    x[0] = rocs::interval(8.0, 8.2);
    x[1] = rocs::interval(0, 0.15);
    x[2] = rocs::interval(0.3927, 0.5890);

    std::vector<std::vector<double>> U { {0.6, -0.9}, {0.9, -0.9}, {0.9, -0.6}, {0.9, -0.3} };

    double tau = 0.3;
    
    car *ptrCar = new car(U, tau);

    std::vector<rocs::ivec> y = (*ptrCar)(x);

    for(int i = 0; i < y.size(); ++i )
	std::cout << y[i] << '\n';

    delete ptrCar;
    
    BOOST_CHECK(true);
    
}


/**
 * test case 3: inverted pendulum
 */
BOOST_AUTO_TEST_CASE(test_ipdl)
{
    /* state and input space */
    const int XD = 2;
    const int UD = 1;
    double xlb[XD] = {-0.25, -0.05};
    double xub[XD] = {0.25, 0.05};
    double ulb[UD] = {-0.1};
    double uub[UD] = {0.1};
    double mu[UD] = {0.02};
    double tau = 0.01;  // sampling time
    double dt = 0.002;  // ode integration step


    /* define control problem */
    ipdl *ptrIpdl = new ipdl(UD, ulb, uub, mu, tau, dt);
    
    /* test time map */
    rocs::ivec x1(2);
    x1.setval(0, rocs::interval(0.02,0.03));
    x1.setval(1, rocs::interval(-0.02,-0.01));
    std::cout << "Test state:" << x1 << '\n';

    std::vector<rocs::ivec> y1 = (*ptrIpdl)(x1);
    std::vector<double> yc;
    ptrIpdl->_ugrid.print_data();
    
    std::cout << "End-state interval:\n";
    for (int i = 0; i < y1.size(); ++i) {
	yc = y1[i].mid();
	std::cout << y1[i] << ", {"
		  << yc[0] << ',' << yc[1] << "}\n";
	
    }

    delete ptrIpdl;
    BOOST_CHECK(true);
}


/**
 * test case 5: room temperature control
 */
BOOST_AUTO_TEST_CASE(test_tpc)
{
    rocs::ivec x(2);
    x[0] = rocs::interval(24, 25);
    x[1] = rocs::interval(25, 27);

    double TS = 50;

    rocs::input_type U {{1}, {2}, {3}, {4}};  // two modes

    TPC tpc = TPC(U, TS);
    
    std::vector<rocs::ivec> y = tpc(x);

    std::cout << y[0] << '\n';
    std::cout << y[1] << '\n';
    std::cout << y[2] << '\n';
    std::cout << y[3] << '\n';

    BOOST_CHECK(true);
    
}
