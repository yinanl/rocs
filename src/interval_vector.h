/**
 *  interval_vector.h
 *  
 *  An interval vector class with a given dimension size.
 *
 *  Created by Yinan Li on Aug. 08, 2016.
 *  
 *  Hybrid Systems Group, University of Waterloo.
 */




#ifndef _interval_vector_h_
#define _interval_vector_h_

#include <vector>

#include <armadillo>
#include "interval.h"


namespace rocs {

class ivec
{
 private:
  interval* _itvls; /**< A pointer to an interval. */
  int _dim; /**< The dimension of the interval vector. */

 public:

  /**
   * Constructors.
   */
  ivec(): _itvls(NULL), _dim(0) {};
  ivec(int n) { _dim = n; _itvls = new interval[n];}
  /**
   * A copy constructor.
   * y = ivec(x).
   */
  ivec(const ivec &ix);

  
  /**
   * A destructor.
   */
  ~ivec() { delete[] _itvls;}

  /**
   * A copy assignment.
   * y = x.
   */
  ivec& operator=(const ivec&);

  /**
   * A [] operator accessing elements.
   */
  interval& operator[](const int i) const { return _itvls[i];}


  
  /**
   * Get the lower bound of an interval vector.
   * (a1.inf,a2.inf,...,an.inf]), vec = [a1] x [a2] x... [an]
   */
  std::vector<double> getinf() const;
  /**
   * Get the upper bound of an interval vector.
   * (a1.sup, a2.sup,..., an.sup), vec = [a1] x [a2] x... [an]
   */
  std::vector<double> getsup() const;
  int getdim() const { return _dim;}
  void setval(int axis, const interval &x) { _itvls[axis] = x;}

  
  /**
   * An empty test.
   * x = \empty if any dimension is empty.
   */
  bool isempty() const;
  /**
   * Get the width of each dimension.
   * (w[a1],w[a2],...w[an]), vec = [a1] x [a2] x... [an]
   */
  std::vector<double> width() const;
  /**
   * Get the radius of each dimension.
   * (w[a1]/2,w[a2]/2,...w[an]/2), vec = [a1] x [a2] x... [an]
   */
  std::vector<double> radius() const;
  /**
   * Get the center point of an interval vector.
   * (mid[a1],mid[a2],...mid[an]), vec = [a1] x [a2] x... [an]
   */
  std::vector<double> mid() const;
  /**
   * Get the maximum width of all dimensions.
   * max([a1],[a2],...[an]), vec = [a1] x [a2] x... [an]
   */
  double maxwidth() const;
  /**
   * Determine the dimension which has the maximum width.
   */
  int maxdim() const;


  /**
   * A test to see if an interval vector is outside the other.
   * x.isout(y), 1 if x, y are disjoint
   */
  bool isout(const ivec&) const;
  bool isout(const std::vector<double>& ) const;
  /**
   * A test to see if an interval vector is inside the other.
   * x.isin(y), 1 if y is contained in x
   */
  bool isin(const ivec&) const;
  bool isin(const std::vector<double>& ) const;

  
  /**
   * Intersection of two interval vectors.
   */
  friend ivec intersect(const ivec&, const ivec&);
  /**
   * Interval hull of two interval vectors.
   */
  friend ivec hull(const ivec&, const ivec&);
  

  /**
   * Operation overloads between interval vectors.
   */
  friend bool operator==(const ivec&, const ivec&);
  friend bool operator!=(const ivec&, const ivec&);

  friend ivec operator-(const ivec&, const ivec&);
  friend ivec operator-(const double, const ivec&);
  friend ivec operator-(const ivec&, const double);
  friend ivec operator-(const std::vector<double>&, const ivec&);
  friend ivec operator-(const ivec&, const std::vector<double>&);

  friend ivec operator+(const ivec&, const ivec&);
  friend ivec operator+(const double, const ivec&);
  friend ivec operator+(const ivec&, const double);
  friend ivec operator+(const std::vector<double>&, const ivec&);
  friend ivec operator+(const ivec&, const std::vector<double>&);

  friend ivec operator*(const double, const ivec&);
  friend ivec operator*(const ivec&, const double);

  
  /**
   * Bisect x along the axis dimension.
   * lowerhalf, upperhalf.
   */
  friend ivec lowerhalf(const ivec&, const int);
  friend ivec upperhalf(const ivec&, const int);

  
  /**
   * A linear affine operation on an interval vector.
   * y = A[x] + b (b=0 -> y = A[x])
   * y = A * xc + |A|*|xr|, x = xc + [xr]
   */
  friend ivec linmap(const arma::mat&, const arma::vec&, const ivec&);

  
  /**
   * I/O 
   */
  friend std::ostream& operator<<(std::ostream&, const ivec&);

};


} // namespace rocs

#endif
