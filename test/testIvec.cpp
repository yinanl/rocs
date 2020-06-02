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
  test case 1: constructors and printing
*/
BOOST_AUTO_TEST_CASE(test_initialization)
{
    rocs::ivec a(3);  // the constructor ivec(int n)
    std::cout << a << std::endl;
    a.setval(0, rocs::interval(-1, 1)); 
    a[1] = rocs::interval(-0.1, 0.01);
    a[2] = rocs::interval(0, 2);
    std::cout << a << std::endl;
    rocs::ivec b{rocs::interval(-1, 1),
		 rocs::interval(-0.1, 0.01),
		 rocs::interval(0, 2)};  // the initialize_list constructor
    rocs::ivec c(b);  // copy constructor
    rocs::ivec d = b;  // assignment constructor
    rocs::ivec e;
    std::cout << e << std::endl;
    
    BOOST_CHECK_EQUAL(a, b);
    BOOST_CHECK_EQUAL(a, c);
    BOOST_CHECK_EQUAL(a, d);
    BOOST_CHECK(e.isempty());
}

/*
  test case 2: property access functions
*/
BOOST_AUTO_TEST_CASE(test_property_access_functions)
{
    rocs::ivec b{rocs::interval(-1, 1),
		 rocs::interval(-0.1, 0.01),
		 rocs::interval(0, 2)};
    rocs::ivec f{rocs::interval(1, 0),
		 rocs::interval(rocs::NINF, rocs::PINF)};
    
    /* test isempty() */
    BOOST_CHECK(!b.isempty());
    BOOST_CHECK(f.isempty());
    
    /* test setters and getters */
    BOOST_CHECK_EQUAL(3, b.getdim());
    BOOST_CHECK_CLOSE(2., b.maxwidth(), 1e-9);
    BOOST_CHECK_EQUAL(0, b.maxdim());
    
    std::vector<double> infexp{-1., -0.1, 0.};
    std::vector<double> supexp{1., 0.01, 2.};
    std::vector<double> rexp{1., 0.055, 1.};
    
    // BOOST_CHECK_EQUAL(infexp, b.getinf());
    // BOOST_CHECK_EQUAL(supexp, b.getsup());
    // BOOST_CHECK_EQUAL(rexp, b.radius());
    BOOST_CHECK(infexp == b.getinf());
    BOOST_CHECK(supexp == b.getsup());
    BOOST_CHECK(rexp == b.radius());
    std::vector<double> brad = b.radius();
    for(int i = 0; i < 3; ++i)
	std::cout << brad[i] << ',';
    std::cout << '\n';

}

/*
  test case 3: test set operations
*/
BOOST_AUTO_TEST_CASE(test_set_operations)
{
    rocs::ivec a{rocs::interval(-1, 1),
		 rocs::interval(-0.1, 0.01),
		 rocs::interval(0, 2)};
    rocs::ivec x{rocs::interval(1.11, 2.3),
		 rocs::interval(-0.5, 0),
		 rocs::interval(0, 2)};
    rocs::ivec y{rocs::interval(-1, 1),
		 rocs::interval(-0.1, -0.08),
		 rocs::interval(0.1, 1.15)};

    rocs::ivec z = y;
    z[1] = rocs::interval(2, 1);
    BOOST_CHECK(a.isout(x));
    BOOST_CHECK(a.isin(y));
    BOOST_CHECK(a.isout(z));

    std::vector<double> s1{0, 0, 0};
    std::vector<double> s2{1, 1, 1};
    BOOST_CHECK(a.isin(s1));
    BOOST_CHECK(a.isout(s2));

    
    rocs::ivec b{rocs::interval(0.3, 0.5),
		 rocs::interval(0.1, 2),
		 rocs::interval(-3, -1)};
    rocs::ivec c{rocs::interval(4, 7),
		 rocs::interval(3, 4)};
    rocs::ivec d{rocs::interval(0.4, 0.8),
		 rocs::interval(0.5, 1),
		 rocs::interval(-3.2, -2.9)};
    rocs::ivec aorb{rocs::interval(-1, 1),
		    rocs::interval(-0.1, 2),
		    rocs::interval(-3, 2)};
    rocs::ivec bandd{rocs::interval(0.4, 0.5),
		     rocs::interval(0.5, 1),
		     rocs::interval(-3, -2.9)};
    
    BOOST_CHECK(intersect(a, b).isempty());
    BOOST_CHECK(hull(a, b) == aorb);
    BOOST_CHECK(intersect(b, d) == bandd);
}

/*
  test case 4: test bisections
*/
BOOST_AUTO_TEST_CASE(test_bisections)
{
    rocs::ivec a{rocs::interval(-1, 1),
		 rocs::interval(-0.1, 0.01),
		 rocs::interval(0, 2)};

    rocs::ivec left = lowerhalf(a, 2);
    rocs::ivec right = upperhalf(a, 2);

    rocs::ivec refl{rocs::interval(-1, 1),
		    rocs::interval(-0.1, 0.01),
		    rocs::interval(0, 1)};
    rocs::ivec refr{rocs::interval(-1, 1),
		    rocs::interval(-0.1, 0.01),
		    rocs::interval(1, 2)};
    
    BOOST_CHECK_EQUAL(left, refl);
    BOOST_CHECK_EQUAL(right, refr);
}

/*
  test case 4: test operator overloads
*/
BOOST_AUTO_TEST_CASE(test_addsub_overloads)
{
    rocs::ivec a{rocs::interval(-1, 1),
		 rocs::interval(-0.1, 0.01),
		 rocs::interval(0, 2)};
    rocs::ivec b{rocs::interval(0.3, 0.5),
		 rocs::interval(0.1, 2),
		 rocs::interval(-3, -1)};
    rocs::ivec c{rocs::interval(4, 7),
		 rocs::interval(3, 4)};
    rocs::ivec d = {rocs::interval(-1, 1),
		    rocs::interval(-0.1, 0.01),
		    rocs::interval(0, 2)};

    double val = 0.04;
    rocs::ivec abadd{rocs::interval(-0.7, 1.5),
		     rocs::interval(0., 2.01),
		     rocs::interval(-3., 1.)};
    rocs::ivec absub{rocs::interval(-1.5, 0.7),
		     rocs::interval(-2.1, -0.09),
		     rocs::interval(1., 5.)};
    rocs::ivec avadd{rocs::interval(-0.96, 1.04),
		     rocs::interval(-0.06, 0.05),
		     rocs::interval(0.04, 2.04)};
    rocs::ivec avsub{rocs::interval(0.26, 0.46),
		     rocs::interval(0.06, 1.96),
		     rocs::interval(-3.04, -1.04)};
    
    rocs::ivec r = a - c;
    BOOST_CHECK_EQUAL(a+b, abadd);
    BOOST_CHECK_EQUAL(a-b, absub);
    BOOST_CHECK_EQUAL(a+val, avadd);
    BOOST_CHECK_EQUAL(b-val, avsub);
    BOOST_CHECK(r.isempty());
    a += b;
    BOOST_CHECK_EQUAL(a, abadd);
    d -= b;
    BOOST_CHECK_EQUAL(d, absub);
}
