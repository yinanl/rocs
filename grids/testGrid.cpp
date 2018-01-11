/**
 *  testGrid.cpp
 *
 *  Boost test of grid class.
 *
 *  Created by Yinan Li on Jan. 21, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Grid
#include <boost/test/unit_test.hpp>


#include "grid.h"


/*
 * test case 1: initialization of grid
 */
BOOST_AUTO_TEST_CASE(test_init_grid)
{
    const int XD = 3;
    
    /* define the state and input space by upper/lower bounds */
    double xlb[XD] = {0, 0, -M_PI};
    double xub[XD] = {10, 10, M_PI};
    double eta[XD] = {0.2, 0.2, 2*M_PI/35.0};

    grid xg(XD, eta, xlb, xub);
    xg.gridding();
    xg.print_info();
    xg.write2mat_data("grids_unicycle.mat", "xrv");
    
    // grid xg1;
    // xg1.gridding(XD, eta, xlb, xub);
    // xg1.print_data();

    double rp[XD] = {0.2, 0.2, 0.2};  // rate of eta
    std::vector< std::vector<double> > sub = xg.subgridding(xg._data[0], rp);
    for (size_t row = 0; row < sub.size(); ++row) {

	std::cout << row << " ";

	for (int col = 0; col < XD; ++col) {

	    std::cout << sub[row][col] << ' ';
	}

	std::cout << '\n';
    }
}


/*
 * test case 2: multidimensional grid test
 */
BOOST_AUTO_TEST_CASE(test_grid_subset_nd)
{
    const int XD = 3;
    double xlb[XD] = {1, -0.4, -5};
    double xub[XD] = {6, 1.2, -3};
    
    /* gridding 1 */
    double eta[XD] = {1, 0.1, 0.5};
    grid xg(XD, eta, xlb, xub);
    // xg.gridding();
    // xg.print_info();
    // xg.print_data();

    /* area 1: partially outside */
    ivec box(XD);
    box[0] = interval(2.6, 3.1);
    box[1] = interval(1.1, 1.23);
    box[2] = interval(-4, -3.9);

    std::vector<size_t> act1 = xg.subset(box, false, false);
    std::vector<size_t> ref(2);
    ref[0] = 296;
    ref[1] = 302;
    // ref[0] = 236;
    // ref[1] = 237;

    BOOST_CHECK(act1 == ref);
    BOOST_CHECK(xg.subset(box, false, true).empty());
    BOOST_CHECK(xg.subset(box, true, false).empty());


    /* gridding 2  */
    double eta2[XD] = {1.3, 0.14, 0.4};
    grid xg2(XD, eta2, xlb, xub);
    xg2.gridding();
    // xg2.print_info();
    // xg2.print_data();

    /* area 2: fully inside */
    ivec box2(XD);
    box2[0] = interval(2.6, 3.1);
    box2[1] = interval(1.1, 1.2);
    box2[2] = interval(-4, -3.9);
    std::vector<size_t> act2 = xg2.subset(box2, true, false);
    std::vector<size_t> ref2(1);
    ref2[0] = 189;
    BOOST_CHECK(act2 == ref2);
    BOOST_CHECK(xg2.subset(box2, true, true).empty());

    /* area 3: fully inside, no empty */
    ivec box3(XD);
    box3[0] = interval(2.9, 5.4);
    box3[1] = interval(-0.08, 0.42);
    box3[2] = interval(-3.6, -3.2);
    std::vector<size_t> act31 = xg2.subset(box3, true, false);
    std::vector<size_t> act32 = xg2.subset(box3, true, true);
    size_t d1[] = {201, 202, 203, 205, 206, 207, 209, 210, 211, 213, 214, 215, 217, 218, 219};
    size_t d2[] = {206, 210, 214};
    std::vector<size_t> ref31(d1, d1+15);
    std::vector<size_t> ref32(d2, d2+3);
    BOOST_CHECK(act31 == ref31);
    BOOST_CHECK(act32 == ref32);
}


/*
 * test case 3: one dimensional grid subset
 */
BOOST_AUTO_TEST_CASE(test_grid_subset_1d)
{
    /* define the state and input space by upper/lower bounds */
    double xlb[1] = {-1.4};
    double xub[1] = {-1.15};
    double eta[1] = {0.03};

    grid xg(1, eta, xlb, xub);
    xg.gridding();
    // xg.print_info();
    // xg.print_data();


    /* area 1: partially out of range */
    ivec box(1);
    box[0] = interval(-1.279, 0);

    std::vector<size_t> act1 = xg.subset(box, false, false);
    size_t d[] = {4, 5, 6, 7 ,8};
    std::vector<size_t> ref(d, d+5);
    BOOST_CHECK(act1 == ref);


    /* area 2: completely out of range */
    ivec box2(1);
    box2[0] = interval(-3, -2.7);
    BOOST_CHECK(xg.subset(box2, false, false).empty());

    
    /* area 3: on the grid boundary */
    ivec box3(1);
    box3[0] = interval(-1.32, -1.25);
    
    std::vector<size_t> ref31(2);
    ref31[0] = 3;
    ref31[1] = 4;
    std::vector<size_t> act31 = xg.subset(box3, true, true);
    BOOST_CHECK(act31 == ref31);
    
    std::vector<size_t> ref32(3);
    ref32[0] = 3;
    ref32[1] = 4;
    ref32[2] = 5;
    std::vector<size_t> act32 = xg.subset(box3, true, false);
    BOOST_CHECK(act32 == ref32);
}


/*
 * test case 4: test ipdl grids
 */
BOOST_AUTO_TEST_CASE(test_grid_ipdl)
{
    const int XD = 2;
    double xlb[XD] = {-0.1, -0.05};
    double xub[XD] = {0.1, 0.05};
    double eta[XD] = {0.001, 0.001};

    grid xg(XD, eta, xlb, xub);
    xg.gridding();
    xg.print_info();
}
