/**
 *  vectorfield.h
 *
 *  Abstract base class for vector field.
 *
 *  Created by Yinan Li on Feb. 18, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _vectorfield_h
#define _vectorfield_h


#include <vector>

#include "../intervals/interval_vector.h"

#include "../grids/grid.h"




typedef std::vector<double> state_type;
typedef std::vector< std::vector<double> > input_type;

/**
 * Base class defining dynamics.
 */
class VFunctor {

 public:
  
  input_type  _uset; /**< a set of controls */
  grid _ugrid;  /**< a grid of controls */
  size_t _unum; /**< number of control values */
  double _tau; /**< sampling time */
  double _dt; /**< ode integrate time step */

  /* constructors */
  VFunctor(){}
  VFunctor(double T) : _tau(T) {}
  VFunctor(double T, double dt) : _tau(T), _dt(dt) {}  /* with integration step */
  VFunctor(size_t un, double T) : _unum(un), _tau(T) {}
  
  VFunctor(input_type U, double T) : _uset(U), _unum(U.size()), _tau(T) {}
  VFunctor(input_type U, double T, double dt) : _uset(U), _unum(U.size()), _tau(T), _dt(dt) {}
  
  VFunctor(grid ug, double T) : _ugrid(ug), _unum(ug._nv), _tau(T) {}
  VFunctor(grid ug, double T, double dt) : _ugrid(ug), _unum(ug._nv), _tau(T), _dt(dt) {}
  
  VFunctor(const int dim, const double lb[], const double ub[],
	   const double mu[], double T) : _ugrid(dim, mu, lb, ub), _tau(T)
  {
    _ugrid.gridding();
    
    _unum = _ugrid._nv;
  }
  VFunctor(const int dim, const double lb[], const double ub[],
	   const double mu[], double T, double dt)
    : _ugrid(dim, mu, lb, ub), _tau(T), _dt(dt)
  {
    _ugrid.gridding();
    
    _unum = _ugrid._nv;
  }

  /**
   * Set a grid of controls to the functor.
   * @param u the grid of controls.
   */
  void setugrid(grid &u) { _ugrid = u;}
  
  /**
   * Interface of interval reachable set computation w.r.t. dynamics.
   *
   * Control inputs u as constant parameters.
   * @param x an interval.
   * @return a list of intervals y determined by u, \f$y=f(x,u)\f$.
   */
  virtual std::vector<ivec> operator()(const ivec &x) = 0;

  /**
   * Default function of real valued solution computation.
   *
   * Control inputs u as constant parameters.
   * @param x a real-valued state.
   * @return a list of solutions y determined by u, \f$y=f(x,u)\f$.
   */
  virtual std::vector<state_type> operator()(const state_type &x) {
    
    if (!_uset.empty()) {
      std::vector<state_type> y(_uset.size(), x);
      return y;
      
    } else if (!_ugrid._data.empty()) {
      std::vector<state_type> y(_ugrid._nv, x);
      return y;

    } else
      assert(false);
    
  }


  virtual ~VFunctor() {}

};


typedef VFunctor* PtrVF;

#endif
