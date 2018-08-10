/**
 *  adutils.h
 *
 *  Using FADBAD++ for auto differentiation:
 *  Template specialization for interval class.
 *
 *  Created by Yinan Li on Mar. 23, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _adutils_h
#define _adutils_h

#include "interval.h"

#include "FADBAD++/fadbad.h"


/* template specialization for auto diff */
template <> struct Op<rocs::interval>
{
  typedef rocs::interval Base;
  static Base myInteger(const int i) { return Base(i); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }
  static rocs::interval myPos(const rocs::interval& x) { return +x; }
  static rocs::interval myNeg(const rocs::interval& x) { return -x; }
  static rocs::interval& myCadd(rocs::interval& x, const rocs::interval& y) { return x+=y; }
  static rocs::interval& myCsub(rocs::interval& x, const rocs::interval& y) { return x-=y; }
  static rocs::interval& myCmul(rocs::interval& x, const rocs::interval& y) { return x*=y; }
  static rocs::interval& myCdiv(rocs::interval& x, const rocs::interval& y) { return x/=y; }
  static rocs::interval myInv(const rocs::interval& x) { return myOne()/x; }
  static rocs::interval mySqr(const rocs::interval& x) { return x*x; }
  static rocs::interval myPow(const rocs::interval& x, const int n) { return pow(x,n); }
  static rocs::interval myPow(const rocs::interval& x, const rocs::interval& y) { return pow(x,y); }
  static rocs::interval mySqrt(const rocs::interval& x) { return sqrt(x); }
  static rocs::interval myLog(const rocs::interval& x) { return log(x); }
  static rocs::interval myExp(const rocs::interval& x) { return exp(x); }
  static rocs::interval mySin(const rocs::interval& x) { return sin(x); }
  static rocs::interval myCos(const rocs::interval& x) { return cos(x); }
  static rocs::interval myTan(const rocs::interval& x) { return tan(x); }
  static rocs::interval myAsin(const rocs::interval& x) { return asin(x); }
  static rocs::interval myAcos(const rocs::interval& x) { return acos(x); }
  static rocs::interval myAtan(const rocs::interval& x) { return atan(x); }
  static bool myEq(const rocs::interval& x, const rocs::interval& y) { return x==y; }
  static bool myNe(const rocs::interval& x, const rocs::interval& y) { return x!=y; }
  static bool myLt(const rocs::interval& x, const rocs::interval& y) { return x<y; }
  static bool myLe(const rocs::interval& x, const rocs::interval& y) { return x<=y; }
  static bool myGt(const rocs::interval& x, const rocs::interval& y) { return x>y; }
  static bool myGe(const rocs::interval& x, const rocs::interval& y) { return x>=y; }
};


#endif
