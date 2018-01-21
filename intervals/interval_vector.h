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
  interval* _itvls;
  int _dim;

 public:

  /* constructors */
  ivec(): _itvls(NULL), _dim(0) {};
  ivec(int n) { _dim = n; _itvls = new interval[n];}
  ivec(const ivec &ix); // copy constructor

  
  /* destructor */
  ~ivec() { delete[] _itvls;}

  ivec& operator=(const ivec&); //assignment
  interval& operator[](const int i) const { return _itvls[i];} //access element


  /* getters and setters */
  std::vector<double> getinf() const;
  std::vector<double> getsup() const;
  int getdim() const { return _dim;}
  void setval(int axis, const interval &x) { _itvls[axis] = x;}


  /* basic properties */
  bool isempty() const;

  std::vector<double> width() const;
  std::vector<double> radius() const;
  std::vector<double> mid() const;

  double maxwidth() const;
  int maxdim() const;


  /* intersection tests */
  bool isout(const ivec&) const;
  bool isout(const std::vector<double>& ) const;
  bool isin(const ivec&) const;
  bool isin(const std::vector<double>& ) const;


  /* set operations */
  friend ivec intersect(const ivec&, const ivec&);
  friend ivec hull(const ivec&, const ivec&);
  

  /* operation overload */
  friend bool operator==(const ivec&, const ivec&);
  friend bool operator!=(const ivec&, const ivec&);

  friend ivec operator-(const ivec&, const ivec&);
  friend ivec operator-(const double, const ivec&);
  friend ivec operator-(const ivec&, const double);

  friend ivec operator+(const ivec&, const ivec&);
  friend ivec operator+(const double, const ivec&);
  friend ivec operator+(const ivec&, const double);

  friend ivec operator*(const double, const ivec&);
  friend ivec operator*(const ivec&, const double);

  /* bisections */
  friend ivec lowerhalf(const ivec&, const int);
  friend ivec upperhalf(const ivec&, const int);


  /* linear operation */
  friend ivec linmap(const arma::mat&, const arma::vec&, const ivec&);
  
  /* I/O */
  friend std::ostream& operator<<(std::ostream&, const ivec&);

};


} // namespace rocs

#endif
