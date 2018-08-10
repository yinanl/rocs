/**
 *  testIvec.cpp
 *
 *  Test interval vector class
 *
 *  Created by Yinan Li on August 19, 2016.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE IntervalVectorClass


#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>
#include <cmath>
#include "src/interval_vector.h"


/*
  test case 1: initialization
*/
BOOST_AUTO_TEST_CASE(test_init_basics)
{
    rocs::ivec a; // test constructors
    rocs::ivec b(3);

    b.setval(0, rocs::interval(-1, 1)); // test setting elements in 2 ways
    b[1] = rocs::interval(-0.1, 0.01);
    b[2] = rocs::interval(0, 2);
    
    rocs::ivec c(b); //test copy constructor, equal to "ivec d = b;"
    rocs::ivec d(3);
    d = b; //test assignment
    
    rocs::ivec e(2);
    e[0] = rocs::interval(1, 0);
    e[1] = rocs::interval(rocs::NINF, rocs::PINF);

    /* test isempty() */
    BOOST_CHECK(a.isempty()); 
    BOOST_CHECK(!b.isempty());
    BOOST_CHECK(!c.isempty());
    BOOST_CHECK(!d.isempty());
    BOOST_CHECK(e.isempty());

    /* test setters and getters */
    BOOST_CHECK(3 == b.getdim());
    BOOST_CHECK(2 == b.maxwidth());
    BOOST_CHECK(0 == b.maxdim());
    std::cout << b.maxwidth() << '\n';
    
    double binfe[] = {-1, -0.1, 0};
    std::vector<double> infexp(binfe, binfe + sizeof(binfe)/sizeof(double));
    double bsupe[] = {1, 0.01, 2};
    std::vector<double> supexp(bsupe, bsupe + sizeof(bsupe)/sizeof(double));

    std::vector<double> rexp(3);
    rexp[0] = 1;
    rexp[1] = 0.055;
    rexp[2] = 1;

    BOOST_CHECK(infexp == b.getinf());
    BOOST_CHECK(supexp == b.getsup());
    BOOST_CHECK(rexp == b.radius());
    std::vector<double> brad = b.radius();
    for(int i = 0; i < 3; ++i)
	std::cout << brad[i] << ',';
    std::cout << '\n';

}

/*
  test case 2: printing
*/
BOOST_AUTO_TEST_CASE(test_printing)
{
    rocs::ivec a(3);

    a.setval(0, rocs::interval(-1, 1)); 
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 2);
    
    std::cout << a << std::endl;
}

/*
  test case 3: test accessing elements
*/
BOOST_AUTO_TEST_CASE(test_intersection)
{
    rocs::ivec a(3);

    a.setval(0, rocs::interval(-1, 1));
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 2);
    
    rocs::ivec x(3);
    rocs::ivec y(3);

    x[0] = rocs::interval(1.11, 2.3);
    x[1] = rocs::interval(-0.5, 0);
    x[2] = rocs::interval(0, 2);

    y[0] = rocs::interval(-1, 1);
    y[1] = rocs::interval(-0.1, -0.08);
    y[2] = rocs::interval(0.1, 1.15);

    rocs::ivec z = y;
    z[1] = rocs::interval(2, 1);

    BOOST_CHECK(a.isout(x));
    BOOST_CHECK(a.isin(y));
    BOOST_CHECK(a.isout(z));

    
    double arr1[] = {0, 0, 0}; // inside
    std::vector<double> s1(arr1, arr1 + sizeof(arr1)/sizeof(double));
    double arr2[] = {1, 1, 1}; // outside
    std::vector<double> s2(arr2, arr2 + sizeof(arr2)/sizeof(double));

    BOOST_CHECK(a.isin(s1));
    BOOST_CHECK(a.isout(s2));
}

/*
  test case 4: test bisections
*/
BOOST_AUTO_TEST_CASE(test_bisections)
{
    rocs::ivec a(3);

    a[0] = rocs::interval(-1, 1);
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 2);

    rocs::ivec left = lowerhalf(a, 2);
    rocs::ivec right = upperhalf(a, 2);

    rocs::ivec refr(3);
    rocs::ivec refl(3);
    refl[0] = rocs::interval(-1, 1);
    refl[1] = rocs::interval(-0.1, 0.01);
    refl[2] = rocs::interval(0, 1);
    refr[0] = rocs::interval(-1, 1);
    refr[1] = rocs::interval(-0.1, 0.01);
    refr[2] = rocs::interval(1, 2);

    BOOST_CHECK(left == refl);
    BOOST_CHECK(right == refr);
}

/*
  test case 4: test operator overloads
*/
BOOST_AUTO_TEST_CASE(test_addsub_overloads)
{
    rocs::ivec a(3);
    rocs::ivec b(3);
    rocs::ivec c(2);
    a[0] = rocs::interval(-1, 1);
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 2);

    b[0] = rocs::interval(0.3, 0.5);
    b[1] = rocs::interval(0.1, 2);
    b[2] = rocs::interval(-3, -1);

    c[0] = rocs::interval(4, 7);
    c[1] = rocs::interval(3, 4);

    double val = 0.04;
    rocs::ivec abadd(3);
    rocs::ivec avadd(3);
    rocs::ivec absub(3);
    rocs::ivec avsub(3);
    abadd[0] = rocs::interval(-0.7, 1.5);
    abadd[1] = rocs::interval(0, 2.01);
    abadd[2] = rocs::interval(-3, 1);
    absub[0] = rocs::interval(-1.5, 0.7);
    absub[1] = rocs::interval(-2.1, -0.09);
    absub[2] = rocs::interval(1, 5);
    avadd[0] = rocs::interval(-0.96, 1.04);
    avadd[1] = rocs::interval(-0.06, 0.05);
    avadd[2] = rocs::interval(0.04, 2.04);
    avsub[0] = rocs::interval(0.26, 0.46);
    avsub[1] = rocs::interval(0.06, 1.96);
    avsub[2] = rocs::interval(-3.04, -1.04);

    rocs::ivec r = a - c;
    BOOST_CHECK(a+b == abadd);
    BOOST_CHECK(a-b == absub);
    BOOST_CHECK(a+val == avadd);
    BOOST_CHECK(b-val == avsub);
    BOOST_CHECK(r.isempty());
}


/*
 * test case 5: test set operations
 */
BOOST_AUTO_TEST_CASE(test_set_operators)
{
    rocs::ivec a(3);
    rocs::ivec b(3);
    rocs::ivec c(2);
    rocs::ivec d(3);
    a[0] = rocs::interval(-1, 1);
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 2);

    b[0] = rocs::interval(0.3, 0.5);
    b[1] = rocs::interval(0.1, 2);
    b[2] = rocs::interval(-3, -1);

    c[0] = rocs::interval(4, 7);
    c[1] = rocs::interval(3, 4);

    d[0] = rocs::interval(0.4, 0.8);
    d[1] = rocs::interval(0.5, 1);
    d[2] = rocs::interval(-3.2, -2.9);

    rocs::ivec aorb(3);
    aorb[0] = rocs::interval(-1, 1);
    aorb[1] = rocs::interval(-0.1, 2);
    aorb[2] = rocs::interval(-3, 2);

    rocs::ivec bandd(3);
    bandd[0] = rocs::interval(0.4, 0.5);
    bandd[1] = rocs::interval(0.5, 1);
    bandd[2] = rocs::interval(-3, -2.9);
    
    // ivec r1 = intersect(a, b);
    // ivec r3 = intersect(a, c);
    
    BOOST_CHECK(intersect(a, b).isempty());
    BOOST_CHECK(hull(a, b) == aorb);
    BOOST_CHECK(intersect(b, d) == bandd);
    
}
