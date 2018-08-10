/**
 *  Testival.cpp
 *
 *  Test interval class
 *
 *  Created by Yinan Li on July 20, 2016.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE IntervalClass


#include <boost/test/unit_test.hpp>
//#include <boost/test/unit_test_log.hpp>
#include <cmath>
#include <iomanip>
#include "src/interval.h"


/*
  test case 1: initialization
*/
BOOST_AUTO_TEST_CASE(test_init_basics)
{
    rocs::interval a(2,1);
    rocs::interval b(0.2, 1.5);
    rocs::interval c(b);
    rocs::interval d = b;
    rocs::interval e(rocs::PINF, rocs::PINF);

    /* test isempty() */
    BOOST_CHECK(a.isempty()); 
    BOOST_CHECK(!b.isempty());
    BOOST_CHECK(!c.isempty());
    BOOST_CHECK(!e.isempty());

    /* test equivalence function */
    BOOST_CHECK(c == b);
    BOOST_CHECK(d == b);
    BOOST_CHECK(a != b);

    /* test print function */
//    std::cout << b<< std::endl;
}

/*
  test case 2: + -
*/
BOOST_AUTO_TEST_CASE(test_plus_minus)
{
    using namespace rocs;
    
    interval e(2,1); // empty interval
    
    interval a(0.2, 1.5);
    interval b(-0.3, -0.0001);
    interval c = interval(-2.33, 1.09);
    
    /* test plus */
    BOOST_CHECK(interval(0.4, 3) == a + a);
    BOOST_CHECK(interval(-0.1, 1.4999) == b + a);
    BOOST_CHECK(interval(-0.1, 1.4999) == a + b);
    BOOST_CHECK(interval(NAN, NAN) == a + e);
    BOOST_CHECK(interval(7.2, 8.5) == a + 7);
    BOOST_CHECK(interval(-102.33, -98.91) == -100 + c);
    BOOST_CHECK(interval(NINF, PINF) == c + interval(NINF, PINF)); // support [-oo,+oo]
    
    /* test minus */
    BOOST_CHECK(interval(0.2001, 1.8) == a - b);
    BOOST_CHECK(interval(-1.8, -0.2001) == b - a);
    BOOST_CHECK(interval(-1.39, 2.3299) == b - c);
    BOOST_CHECK(interval(NAN, NAN) == a - e); // support NAN
    BOOST_CHECK(interval(-11.95, -11.6501) == b - 11.65);
    BOOST_CHECK(interval(-104.13, -100.71) == -103.04 - c);
    BOOST_CHECK(interval(NINF, PINF) == c - interval(NINF, PINF)); // support [-oo, +oo]

    //std::cout << PINF - PINF << "\n";
    //std::cout << (PINF == PINF);
    
}

/*
  test case 3: * /
*/
BOOST_AUTO_TEST_CASE(test_mul_div)
{
    using namespace rocs;
    
    interval e(2,1); // empty interval
    
    interval a(0.2, 1.5);
    interval b(-0.3, -0.0001);
    interval c(-2.33, 1.09);
    interval d(0, 1);
    interval g(-1.15, 0);
    
    interval p1(-1, PINF);
    interval p2(0,PINF);
    interval p3(3, PINF);
    
    interval np(NINF, PINF);

    interval n1(NINF, 5);
    interval n2(NINF,0);
    interval n3(NINF, -0.02);

    /* mul */
    BOOST_CHECK(interval(0.04, 2.25) == a * a);
    BOOST_CHECK(interval(-2.33*1.5, 1.09*1.5) == a * c);
    BOOST_CHECK(interval(-0.45, -0.00002) == a * b);
    BOOST_CHECK(interval(-2.33*1.5, 1.09*1.5) == c * a);
    BOOST_CHECK(interval(-2.5397, 5.4289) == c * c);
    BOOST_CHECK(interval(-0.327, 0.699) == c * b);
    BOOST_CHECK(interval(-0.45, -0.00002) == b * a);
    BOOST_CHECK(interval(-0.327, 0.699) == b * c);
    BOOST_CHECK(interval(0.00000001, 0.09) == b * b);
    BOOST_CHECK(interval(-0.3, 0) == d * b); //contain 0
    BOOST_CHECK(interval(-2.33, 1.09) == c * d);
    BOOST_CHECK(interval(-1.15*1.5, 0) == g * a);
    BOOST_CHECK(interval(0, 1.15*0.3) == g * b);
    BOOST_CHECK(interval(-1.15*1.09, 2.33*1.15) == c * g);
    BOOST_CHECK(interval(-1.15*1.09, 2.33*1.15) == g * c);
    BOOST_CHECK(interval(NAN, NAN) == e * a); //empty
    BOOST_CHECK(interval(-1.5, PINF) == p1 * a); //infinity
    BOOST_CHECK(interval(0, PINF) == p2 * a);
    BOOST_CHECK(interval(NINF, -0.0003) == b * p3);
    BOOST_CHECK(interval(NINF, 7.5) == a * n1);
    BOOST_CHECK(interval(0, PINF) == n2 * b);
    BOOST_CHECK(interval(NINF, PINF) == c * n3);
    BOOST_CHECK(interval(NINF, PINF) == p1 * n1);
    BOOST_CHECK(interval(NINF, PINF) == n1 * n2);
    BOOST_CHECK(interval(NINF, PINF) == n1 * n3);
    BOOST_CHECK(interval(0, PINF) == n2 * n3);
    BOOST_CHECK(interval(NINF, PINF) == np * np);
    BOOST_CHECK(interval(0, 0) == 0 * a);
    BOOST_CHECK(interval(0, 0) == a * 0);
    BOOST_CHECK(interval(0, 0) == 0 * p1);
    BOOST_CHECK(interval(0, 0) == 0 * np);
    BOOST_CHECK(interval(-1, PINF) == 1 * p1);
    BOOST_CHECK(interval(0.00005, 0.15) == -0.5 * b);
    BOOST_CHECK(interval(0, PINF) == -0.5 * n2);


    /* div */
    BOOST_CHECK(interval(0.2/1.5, 1.5/0.2) == a / a);
    BOOST_CHECK(interval(-15000, -0.2/0.3) == a / b);
    BOOST_CHECK(interval(-0.3/0.2, -0.0001/1.5) == b / a);
    BOOST_CHECK(interval(NINF, PINF) == a / c);
    BOOST_CHECK(interval(-2.33/0.2, 1.09/0.2) == c / a);
    BOOST_CHECK(interval(-10900, 23300) == c / b);
    BOOST_CHECK(interval(0, 5) == d / a);
    BOOST_CHECK(interval(0.2, PINF) == a / d);
    BOOST_CHECK(interval(-5.75, 0) == g / a);
    BOOST_CHECK(interval(NINF, -0.2/1.15) == a / g);
    BOOST_CHECK(interval(0, 11500) == g / b);
    BOOST_CHECK(interval(NINF, PINF) == d / g);
    BOOST_CHECK(interval(NAN, NAN) == a / interval(0, 0));
    BOOST_CHECK(interval(0, 0) == 0 / a);
    BOOST_CHECK(interval(0, 0) == 0 / b);
    BOOST_CHECK(interval(NAN, NAN) == 0 / c); // num = 0, 0\in den
    BOOST_CHECK(interval(NAN, NAN) == 0 / g);
    BOOST_CHECK(interval(NAN, NAN) == 0 / d);
    BOOST_CHECK(interval(1/1.5, 5) == 1 / a); // num > 0
    BOOST_CHECK(interval(-10000, -1/0.3) == 1 / b);
    BOOST_CHECK(interval(NINF, PINF) == 1 / c);
    BOOST_CHECK(interval(1, PINF) == 1 / d);
    BOOST_CHECK(interval(NINF, -1/1.15) == 1 / g);
    BOOST_CHECK(interval(-5, -1/1.5) == -1 / a);
    BOOST_CHECK(interval(1/0.3, 10000) == -1 / b);
    BOOST_CHECK(interval(NINF, PINF) == -1 / c);
    BOOST_CHECK(interval(NINF, -1) == -1 / d);
    BOOST_CHECK(interval(1/1.15, PINF) == -1 / g);
}

/*
  test case 4: sin(), cos(), tan(), atan()
*/
BOOST_AUTO_TEST_CASE(test_trigonometric_fcns)
{
    using namespace rocs;
    
    /* test sin() [cos() is implemented by sin()] */
    BOOST_CHECK(interval(-1, sin(-1.1)) == sin(interval(-1.6, -1.1))); //wid = 0.5
    BOOST_CHECK(interval(sin(-1), sin(-0.5)) == sin(interval(-1, -0.5)));
    BOOST_CHECK(interval(sin(-0.2), sin(0.3)) == sin(interval(-0.2, 0.3))); //0
    std::cout << sin(interval(-0.2, 0.3)) << '\n';
    BOOST_CHECK(interval(-1, sin(5)) == sin(interval(4.5, 5))); //3pi/2
    
    BOOST_CHECK(interval(sin(-1), sin(0.57)) == sin(interval(-1, 0.57))); //wid = 1.57
    BOOST_CHECK(interval(0, 1) == sin(interval(0, 1.571))); //pi/2
    BOOST_CHECK(interval(-1, sin(3.14)) == sin(interval(3.14, 4.72))); //3pi/2
    
    BOOST_CHECK(interval(-1, sin(-3.2)) == sin(interval(-3.2, -3.2+PIIVAL))); //wid = pi
    std::cout << sin(interval(-3.2, -3.2+PIIVAL)) << '\n';
    BOOST_CHECK(interval(-1, 0) == sin(interval(-PIIVAL, 0)));
    BOOST_CHECK(interval(-1, 1) == sin(interval(-PIIVAL/2, PIIVAL/2)));
    BOOST_CHECK(interval(sin(-1), 1) == sin(interval(-1, -1+PIIVAL)));

    BOOST_CHECK(interval(sin(-1), 1) == sin(interval(-1, 3.4))); //wid = 4.3
    BOOST_CHECK(interval(-1, 1) == sin(interval(1.5, 5.8)));
    BOOST_CHECK(interval(-1, sin(7.4)) == sin(interval(3.1, 7.4)));
    BOOST_CHECK(interval(-1, 1) == sin(interval(7, 11.3)));

    BOOST_CHECK(interval(-1, 1) == sin(interval(0, PI2IVAL))); //wid = 2pi
    BOOST_CHECK(interval(-1, 1) == sin(interval(PI2IVAL, 2*PI2IVAL)));
    BOOST_CHECK(interval(-1, 1) == sin(interval(-0.5, -0.5+PI2IVAL)));

    BOOST_CHECK(interval(-1, 1) == sin(interval(PIIVAL, 5*PIIVAL))); //wid = 4pi
    BOOST_CHECK(interval(-1, 1) == sin(interval(-6, 1))); //wid = 7


    /* tan() */
    BOOST_CHECK(interval(NINF, PINF) == tan(interval(-1.6, 2.54))); //wid = 4.14 > pi
    // wid = 3.14 < pi, pi/2 \in [inf, sup]
    BOOST_CHECK(interval(NINF, PINF) == tan(interval(-5.14, -2)));
    BOOST_CHECK(interval(tan(4.3), tan(4.5)) == tan(interval(4.3, 4.5)));  // k+=1 branch
    std::cout << tan(interval(4.3, 4.5)) << '\n';
    BOOST_CHECK(interval(NINF, PINF) == tan(interval(4.3, 4.8)));  // k+=1 branch
    BOOST_CHECK(interval(tan(-2), tan(-1.7)) == tan(interval(-2, -1.7)));  // k+=1 branch
    std::cout << tan(interval(-2, -1.7)) << '\n';
}


/*
  test case 5: sqr(), sqrt(), power()
*/
BOOST_AUTO_TEST_CASE(test_power_root)
{
    using namespace rocs;
    
    BOOST_CHECK(interval(0.01, 25) == sqr(interval(-5, -0.1))); //sqr
    BOOST_CHECK(interval(0.01, 25) == sqr(interval(0.1, 5)));
    BOOST_CHECK(interval(0, 25) == sqr(interval(-5, 0.1)));
    BOOST_CHECK(interval(0, 0.25) == sqr(interval(0, 0.5)));
    BOOST_CHECK(interval(0, 0.25) == sqr(interval(-0.5, 0)));

    BOOST_CHECK(pow(interval(0, 0.5), 2) == sqr(interval(0, 0.5))); //pow
    BOOST_CHECK(pow(interval(0, 0.5), 2) == sqr(interval(0, 0.5)));
    BOOST_CHECK(interval(-1, 8) == pow(interval(-1, 2), 3));
    BOOST_CHECK(interval(0, 1) == pow(interval(-1, 1), 6));
    BOOST_CHECK(interval(-1, 2) == pow(interval(-1, 2), 1));

    BOOST_CHECK(interval(1, 2) == sqrt(interval(1, 4)));//sqrt
    BOOST_CHECK(interval(0, 2) == sqrt(interval(-1, 4)));
    BOOST_CHECK(interval(0, 0) == sqrt(interval(-1, 0)));
}


/*
  test case 6: intersect, hull
*/
BOOST_AUTO_TEST_CASE(test_inter_hull)
{
    using namespace rocs;
    
    interval a(-0.7, 0.2);
    interval b(0, 0.5);
    interval c(-0.5, 0);
    interval d(-3, -0.6);
    interval e(NINF, -2);
    interval f(0.3, PINF);
    interval g(NINF, PINF);

    
    BOOST_CHECK(interval(NAN, NAN) == intersect(b, d)); //intersect
    BOOST_CHECK(c == intersect(a, c));
    BOOST_CHECK(interval(0, 0.2) == intersect(a, b));
    BOOST_CHECK(interval(0, 0) == intersect(b, c));
    BOOST_CHECK(interval(-3, -2) == intersect(e, d));
    BOOST_CHECK(interval(NAN, NAN) == intersect(e, a));
    BOOST_CHECK(interval(0.3, 0.5) == intersect(b, f));
    BOOST_CHECK(a == intersect(a, g));

    BOOST_CHECK(interval(-3, 0.5) == hull(b, d)); //hull
    BOOST_CHECK(interval(-0.7, 0.5) == hull(a, b));
    BOOST_CHECK(interval(-0.5, 0.5) == hull(b, c));
    BOOST_CHECK(interval(NINF, -0.6) == hull(e, d));
    BOOST_CHECK(interval(NINF, 0.2) == hull(e, a));
    BOOST_CHECK(a == hull(a, c));
    BOOST_CHECK(g == hull(a, g));
}


/*
  test case 7: bisections
*/
BOOST_AUTO_TEST_CASE(test_bisection)
{
    rocs::interval a(-0.7, 0.2);

    rocs::interval left = lowerhalf(a);
    rocs::interval right = upperhalf(a);

    BOOST_CHECK(rocs::interval(-0.7, -0.25) == left);
    BOOST_CHECK(rocs::interval(-0.25, 0.2) == right);
}
