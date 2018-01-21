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

#include "config.h"


namespace rocs {
 
class interval
{
 private:
  double m_inf; /**< the lower bound */
  double m_sup; /**< the upper bound */
        
 public:

  /**
   * Constructor: no arguments.
   */
  interval(): m_inf(0), m_sup(0){};

  /**
   * Constructor
   * @param l lower bound.
   * @param u upper bound.
   */
  interval(double l, double u): m_inf(l), m_sup(u){};

  /**
   * Copy constructor
   * @param a another interval.
   */
  interval(const interval &a){ //copy
    m_inf= a.m_inf;
    m_sup= a.m_sup;
  }

  
  /* getters and setters */
  double getinf() const { return m_inf; }
  double getsup() const { return m_sup; }
  void setinf(const double val) { m_inf= val; }
  void setsup(const double val) { m_sup= val; }


  /* basic properties */
  bool isempty() const;
  double width() const { return m_sup - m_inf; }
  double mid() const { return m_inf + (m_sup - m_inf)/2; }
        
        
  /**
   * Copy assignment 
   */
  interval& operator=(const interval&);

  
  /* interval overlap checks */
  bool isout(const interval &) const; // if a is disjoint with it
  bool isout(const double ) const;
  bool isin(const interval &) const; // if a is contained in it
  bool isin(const double ) const;


  /* unary operators */
  interval operator-() const;


  /* real arithmetic operations overload: +,-,*,/ */
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


  /* elementary functions */
  friend interval sin(const interval&); // y = sin(x)
  friend interval cos(const interval&); // y = cos(x)
  friend interval tan(const interval&); // y = tan(x)
  friend interval atan(const interval&); // y = atan(x)
  
  friend interval exp(const interval&); // y = e^x
  friend interval log(const interval&); // y = log_e(x)
  friend interval log2(const interval&); // y = log_2(x)
  
  friend interval sqrt(const interval&); // y = sqrt(x)
  friend interval sqr(const interval&); // y = x^2
  friend interval power(const interval&, int); // y = x^n
  friend interval nthroot(const interval&, int); // y = x^(1/n)


  /* intersections and unions */
  friend interval intersect(const interval&, const interval&);
  friend interval hull(const interval&, const interval&);
  

  /* bisections */
  friend interval lowerhalf(const interval&);
  friend interval upperhalf(const interval&);
  
  
  /* I/O */
  friend std::ostream& operator<<(std::ostream&, const interval&);

  
};

} // namespace rocs

#endif
