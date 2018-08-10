/**
 *  interval.h
 *
 *  An interval class.
 *
 *  Created by Yinan Li on May 24, 2016.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _interval_h_
#define _interval_h_

#include <iostream>

#include "validated.h"


namespace rocs {
  
  class interval
  {
  private:
    double m_inf; /**< the lower bound */
    double m_sup; /**< the upper bound */
        
  public:

  interval(): m_inf(0), m_sup(0){};
  interval(double l, double u): m_inf(l), m_sup(u){};
  interval(double r): m_inf(r), m_sup(r){};

    /**
     * A copy constructor.
     * @param a another interval.
     */
    interval(const interval &a){ //copy
      m_inf= a.m_inf;
      m_sup= a.m_sup;
    }

    /**
     * A copy assignment.
     */
    interval& operator=(const interval&);

  
    double getinf() const { return m_inf; }
    double getsup() const { return m_sup; }
    void setinf(const double val) { m_inf= val; }
    void setsup(const double val) { m_sup= val; }

    bool isempty() const;
    double width() const { return roundup(m_sup - m_inf); }
    double mid() const { return roundup(m_inf + (m_sup - m_inf)/2.0); }

    /**
     * Interval overlap checks: [x],[y] intervals; a (real value)
     * [x].isout([y]): [x],[y] are disjoint; 
     * [x].isout(a): a is outside of [x];
     * [x].isin([y]): [y] is inside [x];
     * [x].isin(a): a is inside [x].
     */
    bool isout(const interval &) const; // if a is disjoint with it
    bool isout(const double ) const;
    bool isin(const interval &) const; // if a is contained in it
    bool isin(const double ) const;


    /**
     * Negative
     */
    interval operator-() const;
    /**
     * Positive
     */
    interval operator+() const;
    
    interval& operator += (const double x); // to be compatible with Profil (to return a reference)
    interval& operator -= (const double x);
    interval& operator *= (const double x);
    interval& operator /= (const double x);
    interval& operator += (const interval &x);
    interval& operator -= (const interval &x);
    interval& operator *= (const interval &x);
    interval& operator /= (const interval &x);

  
    friend interval operator+(const interval&, const interval&);
    friend interval operator+(const double, const interval&);
    friend interval operator+(const interval&, const double);
    friend interval operator*(const interval&, const interval&);
    friend interval operator*(const double, const interval&);
    friend interval operator*(const interval&, const double);
    friend interval operator-(const interval&, const interval&);
    friend interval operator-(const double, const interval&);
    friend interval operator-(const interval&, const double);
    friend interval operator/(const interval&, const interval&);
    friend interval operator/(const double, const interval&);
    friend interval operator/(const interval&, const double);
  
    friend bool operator==(const interval&, const interval&);
    friend bool operator!=(const interval&, const interval&);
    friend bool operator<(const interval&x, const interval& y);
    friend bool operator<=(const interval&x, const interval& y);
    friend bool operator>(const interval&x, const interval& y);
    friend bool operator>=(const interval&x, const interval& y);

    friend interval sin(const interval&); // y = sin(x)
    friend interval cos(const interval&); // y = cos(x)
    friend interval tan(const interval&); // y = tan(x)
    friend interval atan(const interval&); // y = atan(x)
    friend interval asin(const interval&); // y = asin(x)
    friend interval acos(const interval&); // y = acos(x)
    friend interval exp(const interval&); // y = e^x
    friend interval log(const interval&); // y = log_e(x)
    friend interval log2(const interval&); // y = log_2(x)
    friend interval abs(const interval&); // y=|x|
    friend interval sqrt(const interval&); // y = sqrt(x)
    friend interval sqr(const interval&); // y = x^2
    friend interval pow(const interval&, const int); // y = x^n
    friend interval pow(const interval &x, const interval &y); // r=x^y
    friend interval pow(const interval &x, const double y);
    friend interval nthroot(const interval&, int); // y = x^(1/n)
  

    /**
     * Intersections and unions.
     */
    friend interval intersect(const interval&, const interval&);
    friend interval hull(const interval&, const interval&);
  

    /**
     * A bisection.
     */
    friend interval lowerhalf(const interval&);
    friend interval upperhalf(const interval&);
  
  
    /**
     * I/O
     */
    friend std::ostream& operator<<(std::ostream&, const interval&);
  };

 
  inline interval invhull(double x, double y) {
    double inf, sup;
    if (x<=y) {
      inf = x;     sup = y;
    } else {
      inf = y;     sup = x;
    }
    return interval(inf, sup);
  }
 

} // namespace rocs

#endif
