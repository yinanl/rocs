/**
 *  Local control synthesis based on a product system of an abstraction and an automaton.
 *
 *  Created by Yinan Li on Oct. 13, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _patcher_h
#define _patcher_h

#include <utility>
#include <cassert>
#include <string>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

extern "C" {
#include "buchi.h"
}

#include "definitions.h"
#include "grid.h"
#include "transition.hpp"


namespace rocs {

    // typedef struct BiDGraph {
    // 	NODE_PRODUCT *nodes;
    // 	EDGE *in_edge;
    // 	int *outdeg;
    // 	BOOL *acc;
    // 	size_t *a_pos, *postid;
    // 	int *a_index;
    // } BiDGraph;
    
    /**
     * \brief A class that holds a winning graph of a control problem.
     *
     * A patcher holds:
     * - an fts: the entire solved product transition system (with only winning set), 
     * - a local controller
     */
    class Patcher
    {
    public:
	fts _winfts;  /**< A transition system structure for the winning graph */
	std::vector<long long> _idmap;  /**< An array that maps a nts_product node id to the index in the arrays in _winfts */
	std::vector<long long> _reachstep;  /**< An array of length of # of winning nodes recording the # of steps to reach accepting states */

	std::vector<long long> _encode;  /**< Map from (n0xn1)-based node index to np-based */
	std::vector<long long> _decode;  /**< Map from np-based node index to (n0xn1)-based */
	int _ctlr = -1;  /**< The recomputed safe local control index (>=0). If patching fails, _ctlr is the default value -1 */

	size_t _na;  /**< The number of DBA nodes */
	size_t _nwin;  /**< The number of winning product states */
	

	/**
	 * \brief The default constructor.
	 *
	 * Only use the default constructor. 
	 * All member variables will be assigned later.
	 *
	 * Disable copy/copy assignment and move/move assignment
	 */
	Patcher(){}
	Patcher(Patcher&&) = delete;
	Patcher(const Patcher&) = delete;
	Patcher& operator=(Patcher&&) = delete;
	Patcher& operator=(const Patcher&) = delete;

	/**
	 * \brief Initialize the values in
	 * _winfts, _idmap, _reachstep, _decode, _encode, _na, _nwin
	 */
	void initialize_winning_graph(HEAD &);
	

	/**
	 * \brief Solve for the _ctlr by backward reachability on the subgraph.
	 *
	 * A subgraph (fts) from current node on the _winfts by forward propagation will be constructed.
	 * @param[in] v The index of (x,q) pair
	 * @param[in] h The horizon (or steps) of forward propagation
	 * @param[in] safecnt The address of the bool vector indicating safe controls
	 */
	int solve_local_reachability(size_t v, int h,
				     boost::dynamic_bitset<> &safecnt);
	/**
	 * \brief Similar to solve_local_reachability(), but with avoid set.
	 *
	 * @param[in] v The index of (x,q) pair
	 * @param[in] N The minimun number of target nodes to be visited in forward propagation
	 * @param[in] avoid The set of avoid nodes
	 * @param[in] target The set of target nodes
	 * @param[in] safecnt The address of the bool vector indicating safe controls
	 */
	std::vector< std::pair<size_t, int> >
	solve_local_reachavoid(size_t v, size_t N,
			       std::vector<long long> &avoid,
			       std::vector<long long> &target,
			       boost::dynamic_bitset<> &safecnt);

	/**
	 * \brief Perform Replanning: determine a local area within a radius around a point.
	 *
	 * @param[in] xo The given \f$R^3\f$ point
	 * @param[in] ro The radius
	 * @param[in] x_grid The grid for the state space
	 * @param[in] q the The current automaton state
	 * @param[in] nNodes The total number of automaton states
	 * @return a vector of permisible local state.
	 */
	std::vector<long long> replan_region(const rocs::Rn &xo, const double ro,
						 rocs::grid &x_grid,
						 const rocs::UintSmall &q,
						 const rocs::UintSmall &nNodes);
	
    }; //end class Patcher
    
} //end namespace

#endif
