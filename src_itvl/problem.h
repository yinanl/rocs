/**
 * problem.h
 *
 * A control problem class.
 *
 * Created by Yinan Li on Jan. 04, 2017.
 *
 * Hybrid Systems Group, University of Waterloo.
 */



#ifndef _problems_h_
#define _problems_h_


#include "vectorfield.h"

#include "interval_paver.h"




/* typedef std::vector<ivec> (*fcn)(const ivec &x, */
/* 				 const std::vector< std::vector<double> > &u); */

/**
 * Control problem class.
 *
 * Single goal area, multiple avoiding areas.
 */
class CntlProb
{
public:

  const char *_name;  		/**< problem name (can't change) */
  
  int _xdim;			/**< state dimension */
  int _udim;			/**< input dimension */
  ivec _workspace;		/**< state space */
  PtrVF _vf;			/**< pointer to the dynamics functor */
  /* int _unum; */
  /* std::vector< std::vector<double> > _uspace;  // _uspace[_unum][_udim] */
    
  /* fcn _vf; */
  
  ivec _goal;			/**< the goal area (an interval) */
  std::vector<ivec> _obs;	/**< avoiding areas (an array of intervals) */

    
  /**
   * Construct from an interval vector.
   * @param name control problem name.
   * @param xd dimension of state space.
   * @param ud dimension of input space.
   * @param ws state space (an interval).
   * @param f pointer to the dynamics functor.
   */
  CntlProb(const char *name, const int xd, const int ud,
	   const ivec &ws, PtrVF f)
    : _name(name), _xdim(xd), _udim(ud), _workspace(ws), _vf(f) {}

  /**
   * Construct from arrays of upper and lower bounds.
   * @param name control problem name.
   * @param xd dimension of state space.
   * @param ud dimension of input space.
   * @param lb[] an array of lower bound.
   * @param ub[] an array of upper bound.
   * @param f pointer to the dynamics functor.
   */
  CntlProb(const char *name, const int xd, const int ud,
	   double lb[], double ub[], PtrVF f)
    : _name(name), _xdim(xd), _udim(ud), _vf(f)  // by lower/upper bounds
  {
    ivec xs(_xdim);

    /* assign workspace and constraint */
    for (int i = 0; i < _xdim; ++i) {
      
      xs[i] = interval(lb[i], ub[i]);
    }
    _workspace = xs;
  }

  /**
   * Write case settings to a .mat file.
   * - X: state space.
   * - U: input space.
   * - ts: sampling time.
   * - G: goal interval.
   * - obs: avoid intervals.
   * @param filename
   */
  void write2mat_settings(const char* filename);

  

  /**
   * I/O interface
   */
  friend std::ostream& operator<<(std::ostream&, const CntlProb&);
  
};


typedef CntlProb* PtrCntlProb;




#endif
