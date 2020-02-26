/**
 *  transition.hpp
 *
 *  A transition system class.
 *
 *  Created by Yinan Li on Sept. 05, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _transition_h_
#define _transition_h_

#include <vector>


namespace rocs {
  
class fts {

 public:
    size_t _nx;  /**< number of discrete states */
    size_t _nu;  /**< number of control input */
    size_t _ntrans;  /**< number of transitions */
  
    std::vector<size_t> _idpre;  /**< [_ntrans] a list of pre's */
    std::vector<size_t> _npre;  /**< [Nx x Nu] a list of the number of the pre's 
				   for each state and input */
    std::vector<size_t> _ptrpre;  /**< [Nx x Nu] the index in _npre of the first pre 
				     for state i under input j: _ptrpre[i*Nu+j] */
    std::vector<size_t> _idpost;  /**< [_ntrans] a list of post's */
    std::vector<size_t> _npost;  /**< [Nx x Nu] a list of the number of the post's 
				    for each state and input [n x m] */
    std::vector<size_t> _ptrpost;  /**< [Nx x Nu] the index in _npost of the first post 
				      for state i under input j: _ptrpost[i*Nu+j] */
    std::vector<double> _cost;  /**< [Nx x Nu] a list of costs for every state i under input j */

    /**
     * Constructor: assign all numbers 0.
     */
    fts() : _nx(0), _nu(0), _ntrans(0) {}
    /**
     * Constructor: assign numbers of states and inputs, number of transitions is 0.
     * @param n number of states
     * @param m number of inputs
     */
    fts(size_t n, size_t m) : _nx(n), _nu(m), _ntrans(0),
			      _npre(n*m,0), _ptrpre(n*m,0),
			      _npost(n*m,0), _ptrpost(n*m,0),
			      _cost(n*m, 0) {}

    /**
     * Initialization: assign numbers of states and inputs, number of transitions is 0.
     * @param n number of states
     * @param m number of inputs
     */
    void init(size_t n, size_t m) {
	_nx = n;
	_nu = m;
	_ntrans = 0;
	_npost.resize(n*m, 0);
	_npre.resize(n*m, 0);
	_ptrpost.resize(n*m, 0);
	_ptrpre.resize(n*m, 0);
	_cost.resize(n*m, 0);
    }


    /**
     * Get post states.
     * @param ix a state id.
     * @param ju a input id.
     * @return a list of integers of post state ids from ix under ju.
     */
    std::vector<size_t> get_post(const size_t ix, const size_t ju) {

	// std::vector<size_t> post(_npost[ix*_nu+ju]);
	std::vector<size_t>::iterator ptr = _idpost.begin() + _ptrpost[ix*_nu+ju];
	std::vector<size_t> post(ptr, ptr + _npost[ix*_nu+ju]);

	return post;
    }


    /**
     * Get the costs corresponding to the post states.
     * @param ix a state id.
     * @param ju a input id.
     * @return a list of costs of transitions from ix under ju.
     */
    double get_cost(const size_t ix, const size_t ju) {
	return _cost[ix*_nu+ju];
    }
    

    /**
     * Get predecessor states.
     * @param ix a state id.
     * @param ju a input id.
     * @return a list of integers of preceding state ids leading to kx under ju.
     */
    std::vector<size_t> get_pre(const size_t ix, const size_t ju) {

	std::vector<size_t>::iterator ptr = _idpre.begin() + _ptrpre[ix*_nu+ju];
	std::vector<size_t> pre(ptr, ptr + _npre[ix*_nu+ju]);

	return pre;
    }
  
};


} // namespace rocs

#endif
