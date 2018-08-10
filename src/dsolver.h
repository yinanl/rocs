/**
 *  dsolver.h
 *
 *  A control synthesis solver class over a finite state system.
 *
 *  Created by Yinan Li on Sept. 04, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _dsolver_h_
#define _dsolver_h_


#include <climits>
#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "config.h"
#include "transition.h"


namespace rocs {
  
/**
 * A control synthesis solver class over a (finite) transition system.
 */
class DSolver {
  
 public:
  fts *_ts;  /**< the pointer to a transitionsystem class */
  
  boost::dynamic_bitset<> _winset;  /**< [N] states in the winning set are marked 1 */
  
  std::vector<size_t> _optctlr;  /**< [N] the optimal input for a state */
  std::vector<double> _value;  /**< [N] the value of the optimal input */
  boost::dynamic_bitset<> _leastctlr;  /**< [N*M] all feasible inputs for each state */
  

  /**
   * Constructor: initialize the sizes of member arrays.
   */
 DSolver(fts *ts) : _ts(ts), _winset(_ts->_nx), _optctlr(_ts->_nx),
    _value(_ts->_nx, PINF), _leastctlr(_ts->_nx * _ts->_nu) {}

  /**
   * Solve a reachability game (minimax goal) using djikstra's algorithm.
   *
   * The complexity is of O(c*_ntrans). There are edges might be checked more than once, 
   * because the maximal cost is taken for all \f$k\in (i,j,k), 
   * where i(start state) and j(input) are given.\f$
   * @param[in] target indices of target area.
   * @param[in] avoid indices of avoiding area (DEFAULT is empty).
   * @return an optimal controller in _optctlr member, 
   * and a least restrictive controller in _leastctlr.
   */
  void reachability(std::vector<size_t> &target,
		    std::vector<size_t> avoid = std::vector<size_t>());
  
  /**
   * Solve an invariance game.
   * 
   * @param[in] target indices of target area.
   * @return a least restrictive controller in _leastctlr.
   */
  void invariance(std::vector<size_t> &target);
};


} // namespace rocs

#endif
