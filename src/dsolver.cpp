/**
 *  dsolver.cpp
 *
 *  A control synthesis solver class.
 *
 *  Created by Yinan Li on Sept. 04, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include <queue>
#include <iostream>

#include "dsolver.h"


namespace rocs {

    void DSolver::reachability(std::vector<size_t> &target) {
    	size_t n = _ts->_nx;
    	size_t m = _ts->_nu;
    	std::vector<int> outdeg(_ts->_npost);
    	std::vector<double> pjmax_cost(n*m,0);
	
    	/* initialize _win, _val, and queue */
    	std::queue<size_t> que;
	size_t t;
    	for (size_t i = 0; i < target.size(); ++i) {
	    t = target[i];
    	    _win[t] = true;
	    _leastctlr.set(t*m, m, true);
    	    _value[t] = 0;
    	    que.push(t);
    	}

	/* continuously compute predecessors */
	_nw = target.size();
    	size_t s, p; /* p->s */
    	while (!que.empty()) {
    	    s = que.front();
    	    que.pop();
    	    for (size_t j = 0; j < m; ++j) {
    		/* loop all predecessors of s under input j */
    		for (int l = 0; l < _ts->_npre[s*m+j]; ++l) {
    		    p = _ts->_idpre[_ts->_ptrpre[s*m+j]+l]; /* p(j)->s */
    		    /* evaluate the edge p(j)->s */
		    if (!_win[p]) {
		    	--outdeg[p*m+j];
		    	if (!outdeg[p*m+j]) {
		    	    _leastctlr[p*m+j] = true;
		    	    _win[p] = true;
		    	    que.push(p);
			    ++_nw;
		    	}
		    }
		    // --outdeg[p*m+j];
    		    // if (!outdeg[p*m+j]) {
    		    // 	/* if the out degree of p (under input j) is 0, then p is winning */
    		    // 	if (!_win[p]) {
    		    // 	    // std::cout << p << " is winning.\n";
    		    // 	    _win[p] = true; /* once marked true, can be marked false */
    		    // 	    que.push(p);
    		    // 	}
    		    // 	_leastctlr[p*m+j] = true; /* mark all possible control inputs */
    		    // 	/* compute the optimal cost */
    		    // 	if (_value[p] > _value[s]+_ts->_cost[p*m+j]) {
    		    // 	    _value[p] = _value[s]+_ts->_cost[p*m+j];
    		    // 	    _optctlr[p] = j;
    		    // 	}
    		    // }
    		} /* end loop predecessors of s under input j */
    	    } /* end loop control input */
    	} /* end while loop */
	
    } /* end function reachability */


    void DSolver::invariance(std::vector<size_t> &target) {
	size_t n = _ts->_nx;
	size_t m = _ts->_nu;
	boost::dynamic_bitset<> out(n*m, false);
	bool isout = true;
	
	/* initialize _win:  */
	size_t t;
	for (size_t i = 0; i < target.size(); ++i) {
	    t = target[i];
	    _win[t] = true;
	    _value[t] = 0;
	    _leastctlr.set(t*m, m, true);
	}

	size_t s;
	size_t nt = target.size();
	size_t oldnt = nt + 1;
	while (oldnt > nt) {
	    oldnt = nt;
	    for (size_t i = 0; i < n; ++i) {
		/* consider the state marked winning */
		if (_win[i]) {
		    isout = true;
		    for (size_t j = 0; j < m; ++j) {
			if (_ts->_npost[i*m+j]) {/* if j is a valid control input for i */
			    for(int l = 0; l < _ts->_npost[i*m+j]; ++l) {
				s = _ts->_idpost[_ts->_ptrpost[i*m+j]+l]; /* i,j->s */
				if (!_win[s]) {
				    out[i*m+j] = true;
				    _leastctlr[i*m+j] = false;
				    break;
				}
			    }
			} else {
			    _leastctlr[i*m+j] = false;
			    out[i*m+j] = true;
			}
			isout = isout & out[i*m+j];
		    } /* end loop inputs */
		    if (isout) {
			_win[i] = false;
			--nt;
		    }
		}
	    } /* end loop states */
	} /* end while */

	_nw = nt;
    }


} // namespace rocs
