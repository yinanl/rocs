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
#include "intervals/interval_vector.h"
#include "grids/grid.h"


namespace rocs {

typedef std::vector<double> state_type;
typedef std::vector< std::vector<double> > input_type;

/**
 * Base class defining dynamics.
 */
class VFunctor {

 public:
  
  grid _ugrid;  /**< a grid of controls */
  double _tau; /**< sampling time */
  double _dt; /**< ode integrate time step */

  
  /**
   * Constructors:
   *
   * Initialize the sampling time and number of input values.
   * @param T the sampling time.
   * @param dt the integration step.
   * @param un the number of input values.
   */
  VFunctor(){}
  VFunctor(double T, double dt = 0) : _tau(T), _dt(dt) {}
  
  /**
   * Constructor:
   *
   * Change input_type U to the grid data type.
   * @param U the inputs in the form of input_type.
   * @param T the sampling time.
   */
  VFunctor(input_type U, double T, double dt = 0) : _tau(T), _dt(dt) {
    
    _ugrid._nv = U.size();
    _ugrid._dim = U[0].size();
    _ugrid._data = U;
    
  }
  /**
   * Constructor:
   * 
   * Initialize from input grid.
   * @param ug input grid.
   */
  VFunctor(grid ug, double T, double dt = 0) : _ugrid(ug), _tau(T), _dt(dt) {}
  /**
   * Constructor:
   *
   * Initialize the input grid during functor construction.
   * @param dim dimension of inputs.
   * @param lb lower bound of input values.
   * @param ub upper bound of input values.
   * @param mu grid width.
   * @param T sampling time.
   */
  VFunctor(const int dim, const double lb[], const double ub[],
	   const double mu[], double T, double dt = 0)
    : _ugrid(dim, mu, lb, ub), _tau(T), _dt(dt)
  {
    _ugrid.gridding();
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
    
    if (!_ugrid._data.empty()) {
      std::vector<state_type> y(_ugrid._nv, x);
      return y;

    } else {
      std::cout << "vfunctor() failed: no input data.\n";
      assert(false);
    }
    
  }


  virtual ~VFunctor() {}

};


typedef VFunctor* PtrVF;


} // namespace rocs

#endif
