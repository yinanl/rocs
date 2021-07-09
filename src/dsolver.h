/**
 *  A control synthesis solver class over a finite state system.
 *  Functions are limited.
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
#include <string>

#include "config.h"
#include "transition.hpp"


namespace rocs {
  
    /**
     * \brief A naive abstraction-based engine solver.
     */
    class DSolver {
  
    public:
	fts *_ts;  /**< The pointer to a transitionsystem class */

	size_t _nw;  /**< The number of winning states */
	boost::dynamic_bitset<> _win;  /**< [N] States in the winning set are marked 1 */
  
	std::vector<size_t> _optctlr;  /**< [N] The optimal input for a state */
	std::vector<double> _value;  /**< [N] The value of the optimal input */
	boost::dynamic_bitset<> _leastctlr;  /**< [N*M] All feasible inputs for each state */
  

	/**
	 * \brief A constructor
	 *
	 * Initialize the sizes of member arrays.
	 * @param[in] fts The pointer to a finite transition system
	 */
	DSolver(fts *ts) : _ts(ts), _nw(0), _win(_ts->_nx, false), _optctlr(_ts->_nx, 0),
			   _value(_ts->_nx, PINF), _leastctlr(_ts->_nx * _ts->_nu, false) {}

	/**
	 * \brief Solve a reachability game (minimax goal) using djikstra's algorithm.
	 *
	 * The complexity is of O(c*_ntrans). There are edges might be checked more than once, 
	 * because the maximal cost is taken for all \f$k\in (i,j,k)\f$, 
	 * where i(start state) and j(input) are given.
	 * @param[in] target Indices of target area
	 * @return an optimal controller in _optctlr member, and a least restrictive controller in _leastctlr.
	 */
	void reachability(std::vector<size_t> &target);
	
	/**
	 * \brief Solve an invariance game.
	 * 
	 * @param[in] target Indices of target area
	 * @return a least restrictive controller in _leastctlr.
	 */
	void invariance(std::vector<size_t> &target);

	// /**
	//  * \brief Solve a cobuchi game by the \f$\mu\f$-calculus formula:
	//  *
	//  * \f$\nu Y.\nu X.[\text{Pre}(Y)\cup(B\cap \text{Pre}(X))]\f$
	//  * @param[in] target Indices of target area.
	//  * @param[in] avoid Indices of avoiding area (DEFAULT is empty).
	//  * @return results are recorded in DSolver class members.
	//  */
	// void reach_stay(std::vector<size_t> &target, 
	// 		std::vector<size_t> avoid = std::vector<size_t>());

    };


} // namespace rocs

#endif
