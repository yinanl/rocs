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

    void DSolver::reachability(std::vector<size_t> &target,
			       std::vector<size_t> avoid) {
    
	/* initialize _winset to target */
	std::queue<size_t> que;
	for (int iter = 0; iter < target.size(); ++iter) {
	    _winset[target[iter]] = 1;
	    _value[target[iter]] = 0;
	    que.push(target[iter]);
	}
	/* initialize avoid map */
	boost::dynamic_bitset<> avoidmap(_ts->_nx, 0);
	for (int iter = 0; iter < avoid.size(); ++iter) {
	    avoidmap[avoid[iter]] = 1;
	}

	size_t s;
	while (!que.empty()) {
	    s = que.front();
	    que.pop();

	    for (int j = 0; j < _ts->_nu; ++j) {
	    
		std::vector<size_t> prejs = _ts->get_pre(s, j); // p,j->s
		for (int i = 0; i < prejs.size(); ++i) {
		    double maxpj = 0;  // temporary variable for max(p,j), p,j are fixed.
		
		    /* check if pre should be avoided */
		    if (avoidmap[prejs[i]])
			continue;
		    /* check all posts from p=prejs[i] under control j */
		    bool out = 0;
		    std::vector<size_t> postpj = _ts->get_post(prejs[i], j);
		    std::vector<double> costpj = _ts->get_cost(prejs[i], j);
		    
		    /* postpj must have at least one post prejs[i] */
		    for (int k = 0; k < postpj.size(); ++k) {
			if (!_winset[postpj[k]]) {
			    out = 1;
			    break;
			}
			/* max(p,j) = max_k{_cost(p,j,q)+_value[q]}, q=postpj[k] */
			maxpj = maxpj<_value[postpj[k]]+costpj[k] ? _value[postpj[k]]+costpj[k] : maxpj;
		    
		    } // end check posts of (i,j)
		
		    if (!out) {
			if (!_winset[prejs[i]]){
			    _winset[prejs[i]] = 1;
			    que.push(prejs[i]);
			}
		    
			_leastctlr[prejs[i]*_ts->_nu + j] = 1; // least restrictive control
		    
			if (_value[prejs[i]] > maxpj) { // optimal control: min_j{max(p,j)}
			    _value[prejs[i]] = maxpj;
			    _optctlr[prejs[i]] = j;
			}
		    }
		
		} // end pre's loop
	    } // end inputs loop
	
	} //end while
    
    }


    void DSolver::invariance(std::vector<size_t> &target) {

	/* initialize _winset to target */
	std::queue<size_t> que;
	for (int i = 0; i < target.size(); ++i) {
	    _winset[target[i]] = 1;
	    _value[target[i]] = 0;
	    que.push(target[i]);
	}
    
	size_t s;
	size_t nq = 0;
	bool newloop = 1;
    
	boost::dynamic_bitset<> winmap(_ts->_nx);
	boost::dynamic_bitset<> zeros(_ts->_nx, 0);
	while (nq != que.size() && !que.empty()) {
	    if (newloop) {
		winmap = winmap & zeros;  // set winmap to 0
		nq = que.size();
		newloop = 0;
	    }
	    s = que.front();
	    que.pop();
	
	    /* for a state s: input loop as outer loop */
	    for (int j = 0; j < _ts->_nu; ++j) {
		std::vector<size_t> prejs = _ts->get_pre(s, j); // p,j->s
		/* loop all pre's of state s under input j */
		for (int i = 0; i < prejs.size(); ++i) {
		    /* first make sure that presj[i] is in the winset */
		    if (!_winset[prejs[i]])
			continue;

		    /* check all posts from i under control j */
		    bool out = 1;
		    bool outflag = 0;
		    std::vector<size_t> postpj = _ts->get_post(prejs[i], j);
		    for (int k = 0; k < postpj.size(); ++k) {
			if (!_winset[postpj[k]]) {
			    outflag = 1;
			}
		    } // end check posts of (i,j)
		    if (!outflag)
			out = 0;
		
		    if (!out) {
			_leastctlr[prejs[i]*_ts->_nu + j] = 1;
			winmap[prejs[i]] = 1;
		    }  // assign least restrictive control
		
		} // end pre's loop
	    } // end inputs loop

	    if (que.empty()) {
		swap(winmap, _winset);
		newloop = 1;
		// std::cout << '\n' << "Winning set: ";
		for (int id = 0; id < target.size(); ++id) {
		    if (_winset[target[id]]) {
			que.push(target[id]);
		    }
		}
	    
	    }  // prepare for the next loop when que is empty.
	
	}  //end while
    
    }



} // namespace rocs
