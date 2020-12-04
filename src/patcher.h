/**
 *  patcher.h
 *
 *  Local control synthesis based on a product system of an abstraction and an automaton.
 *
 *  Created by Yinan Li on Oct. 13, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _patcher_h
#define _patcher_h

#include <cassert>
#include <string>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

extern "C" {
#include "buchi.h"
}
#include "transition.hpp"
// #include "hdf5io.h"


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
     * A patcher holds:
     * 1)an fts: the entire solved product transition system (with only winning set), 
     * 2)a local controller
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
	 * Initialize the values in
	 * _winfts, _idmap, _reachstep, _decode, _encode, _na, _nwin
	 */
	void initialize_winning_graph(HEAD &);
	

	/**
	 * Construct a subgraph (fts) from current node on the _winfts by forward propagation.
	 * Solve for the _ctlr by backward reachability on the subgraph.
	 * @param v index of (x,q) pair
	 * @param h horizon (or steps) of forward propagation
	 * @param safecnt the address of the bool vector indicating safe controls
	 */
	int solve_local_reachability(size_t v, int h,
				      boost::dynamic_bitset<> &safecnt);
	
	
    }; //end class Patcher
    
} //end namespace

#endif
