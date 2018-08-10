/**
 *  abstraction.h
 *
 *  An abstraction class.
 *
 *  Created by Yinan Li on Aug. 30, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _abstraction_h_
#define _abstraction_h_

#include <armadillo>
#include <boost/dynamic_bitset.hpp>

#include "transition.h"
#include "system.hpp"


namespace rocs {

    /**
     * Boolean function defining a set.
     */
    typedef bool (*gset)(const ivec &x);
    
    /**
     * Weighting function (of edges) of transitions.
     */
    typedef double (*WGT)(std::vector<double> &x0,
			  std::vector<double> &x1,
			  std::vector<double> &u0);
    
    
    template<typename S>
    class abstraction {
    public:
	grid _x;  /**< a grid of states */

	S *_ptrsys; /**< pointer to the system object */
	WGT _wf;  /**< weighting callback function */
  
	fts _ts;  /**< the transition system to be constructed */
  
	/**
	 * Constructors: default and the other two with vector field functor.
	 *
	 * Assign user defined system dynamics and weighting function to an abstraction.
	 * @param sys the pointer to a system object.
	 */
	abstraction(S *sys):_ptrsys(sys) {}

	/**
	 * State initialization
	 *
	 * - Initialize state grid by bound and grid size.
	 * - Assign the subgridding parameter.
	 * @param n dimension
	 * @param eta[n] an array of grid size.
	 * @param xlb[n] an array of lower bounds.
	 * @param xub[n] an array of upper bounds.
	 */
	void init_state(const int n, const double eta[],
			const double xlb[], const double xub[]) {
    
	    /* _rp.assign(rp, rp+n); */
	    _x.init(n, eta, xlb, xub);
	    _x.gridding();
	}

	/**
	 * Initialize the size of _npre and _npost (n x m).
	 */
	bool init_transitions() {
	    if (_x._nv > 0 && _ptrsys->_ugrid._nv > 0) {
		_ts.init(_x._nv, _ptrsys->_ugrid._nv);
	    } else {
		std::cout << "Transition initialization failed: gridding problem.\n";
		return false;
	    }
	    return true;
	}

	/**
	 * Get discrete state indicies of a given area (the area is fully covered).
	 * @param xlb lower bound of the area.
	 * @param xub upper bound of the area.
	 * @return a list of indicies.
	 */
	std::vector<size_t> get_discrete_states(const double xlb[], const double xub[]) {
	    ivec states(_x._dim);
	    for (int i = 0; i < _x._dim; ++i) {
		states.setval(i, interval(xlb[i], xub[i]));
	    }

	    return _x.subset(states, false, false);  // allow area out of domain, and collect grids intersect the area
	}

	/**
	 * Get discrete state indicies of a given area:
	 * \f$\{x|g(x)=\text{True}\}\f$.
	 *
	 * @param g the function determining whether x is inside the set.
	 * @return a list of indicies.
	 */
	std::vector<size_t> get_discrete_states(gset g);
  

	/**
	 * Assign transitions with robustness margins e1,e2:
	 * @param e1 the robustness margin at the initial point.
	 * @param e2 the robustness margin at the end point.
	 * @return whether construction is successful.
	 */
	bool assign_transitions(const double e1[], const double e2[]);
  
	/**
	 * Assign transitions considering robustness margins:
	 * see assign_transitions()
	 *
	 * @param rp[] the pointer to an array of relative subgridding size
	 * @return whether construction is successful.
	 */
	bool assign_transitions_subgridding(const double rp[]);
    };
    

    template<typename S>
    std::vector<size_t> abstraction<S>::get_discrete_states(gset g) {
	std::vector<size_t> r;
	if (_x._nv>0 && !_x._data.empty()) {
	    ivec x(_x._dim);
	    for (int i = 0; i < _x._nv; ++i) {
		/* set interval x */
		for (int j = 0; j < _x._dim; ++j) {
		    x.setval(j, interval(_x._data[i][j]-_x._gw[j]/2, _x._data[i][j]+_x._gw[j]/2));
		}
		/* test */
		if (g(x)) {
		    r.push_back(i);
		}
	    }
		
	} else {
	    std::cout << "get_discrete_states: a grid of state space hasn't been iniitlized.\n";
	}
	return r;
    }
    

    template<typename S>
    bool abstraction<S>::assign_transitions(const double e1[], const double e2[]) {
	if (!init_transitions())
	    return false;

	int n = _x._dim;
	size_t nx = _x._nv;
	size_t nu = _ptrsys->_ugrid._nv;

	_ts._ntrans = 0;
	// arma::vec r(_x._gw);  // grid width

	ivec ie2(n);
	for (int j = 0; j < n; ++j)
	    ie2.setval(j, interval(-e2[j],e2[j]));
	
	ivec y0(n);
	double xc;
	/* loop state grids */
	for (size_t row = 0; row < nx; ++row) {

	    /* initialize y0 */
	    for (int j = 0; j < n; ++j) {
		xc = _x._data[row][j];
		y0.setval(j, interval(xc-_x._gw[j]/2.0-e1[j], xc+_x._gw[j]/2.0+e1[j]));
	    }
	    
	    /* compute reachable set */
	    std::vector<ivec> yt = _ptrsys->get_reach_set(y0);

	    if (!yt.empty()) {
		for (size_t col = 0; col < nu; ++col) { /* loop inputs */
		    yt[col] = yt[col] + ie2;
		    std::vector<size_t> posts = _x.subset(yt[col], true, false);
		    if (!posts.empty()) { /* assign post transition arrays */
			_ts._npost[row*nu + col] = posts.size();
			_ts._ptrpost[row*nu + col] = _ts._ntrans;
			_ts._ntrans += posts.size();
		
			for (std::vector<size_t>::iterator it = posts.begin(); it != posts.end(); ++it) {
		    
			    _ts._idpost.push_back(*it);
    		    
			    if (_wf) {
				_ts._cost.push_back((*_wf)(_x._data[row], _x._data[*it], _ptrsys->_ugrid._data[col]));
			    }
			    else {  // if no weighting function defined, assign all to 0
				_ts._cost.push_back(0);
			    }
		    
			    /* record the number of pres for state (*it) */
			    _ts._npre[(*it) * nu + col] ++;
			}
		    }  // end transition assignment
		}  // end input loop
	    } // end yt empty check
	    
	}  // end state loop
	
	std::cout << "# of transitions: " << _ts._ntrans << '\n';
	assert(_ts._idpost.size() == _ts._ntrans);

	/* determine pre's by post's: loop _npost and _idpost */
	_ts._idpre.resize(_ts._ntrans);  // initialize the size of pre's
	/* assign _ptrpre */
	size_t sum = 0;
	for (size_t row = 0; row < nx; ++row) {
	    for (size_t col = 0; col < nu; ++col) {
		_ts._ptrpre[row*nu + col] = sum;
		sum += _ts._npre[row*nu + col];
	    }
	}
	/* assign _idpre */
	std::vector<size_t> precount(_ptrsys->_ugrid._nv*_x._nv, 0);
	size_t idtspre;
	for (size_t row = 0; row < nx; ++row) {
	    for (size_t col = 0; col < nu; ++col) {
		for (int ip = 0; ip < _ts._npost[row*nu + col]; ++ip) {
		
		    idtspre = _ts._idpost[_ts._ptrpost[row*nu+col] + ip]*nu + col;
		    _ts._idpre[_ts._ptrpre[idtspre] + precount[idtspre]] = row;
		    precount[idtspre] ++;
		}
	    }
	}
	assert(_ts._idpre.size() == _ts._ntrans);
	return true;
    }


    template<typename S>
    bool abstraction<S>::assign_transitions_subgridding(const double rp[]) {
	if (!init_transitions())
	    return false;
    
	int n = _x._dim;
	size_t nx = _x._nv;
	size_t nu = _ptrsys->_ugrid._nv;
    
	_ts._ntrans = 0;
	boost::dynamic_bitset<> mark(nx); // all 0's by default
	boost::dynamic_bitset<> zeros(nx);
    
	/* compute the number of sub grid points */
	size_t subnv = 1;
	std::vector<double> subgw(n);
	std::vector<size_t> number(n);
	for (int k = 0; k < n; ++k) {
	    number[k] = ceil(1.0 / rp[k]);
	    subgw[k] = _x._gw[k] / number[k];
	    subnv *= number[k];
	}
    
	/* transition computation by interval subgridding: loop states */
	std::vector<double> xmin(n);
	ivec v(n);
	std::vector<size_t> subposts;
	std::vector<std::vector<double>> sub(subnv, std::vector<double> (n));
	std::vector<size_t>::iterator iter;
	int np = 0;
    
	for (size_t row = 0; row < nx; ++row) {
       	    /* compute reachable set by interval subgridding */
	    if (subnv == 1) {  // no subgridding
		sub[0] = _x._data[row];
	    }
	    else {  // subnv > 1
		for (int k = 0; k < n; ++k) {
		    xmin[k] = _x._data[row].at(k) - _x._gw[k]/2. + subgw[k]/2.;
		}
		_x.griddingHelper(sub, xmin, subgw, number, subnv);
	    }
	
	    /* compute post intervals for each subgrid w.r.t. all inputs */
	    std::vector< std::vector<ivec> > ys(subnv, std::vector<ivec> (nu));  // ys[vi][ui]
	    for (int vi = 0; vi < subnv; ++vi) {
	    
		for (int k = 0; k < n; ++k) {
		    v.setval(k, interval(sub[vi][k]-subgw[k]/2, sub[vi][k]+subgw[k]/2));
		} // assign v
		/* get a list of post intervals w.r.t. different inputs */
		ys[vi] = _ptrsys->get_reach_set(v);
	    } // end for loop (subgrid)
	
	
	    for (size_t col = 0; col < nu; ++col) {
		std::set<size_t> posts;
	    
		/* loop subgrids: collect all unique posts */
		for (int vi = 0; vi < subnv; ++vi) {		
		    subposts = _x.subset(ys[vi][col], true, false);
		    if (subposts.empty()) {  // out of domain
			np = 0;
			posts.clear();
			break;  // jump out of the subgrid loop
		    } else {
			posts.insert(subposts.begin(), subposts.end());
		    } // end if
		} // end collecting posts

		/* assign posts */
		if (!posts.empty()) {
		    _ts._npost[row*nu + col] = posts.size();
		    _ts._ptrpost[row*nu + col] = _ts._ntrans;
		    _ts._ntrans += posts.size();
		
		    for (std::set<size_t>::iterator it = posts.begin(); it != posts.end(); ++it) {
			_ts._idpost.push_back(*it);
		    
			if (_wf) {
			    _ts._cost.push_back((*_wf)(_x._data[row], _x._data[*it], _ptrsys->_ugrid._data[col]));
			}
			else {  // if no weighting function defined, assign all to 0
			    _ts._cost.push_back(0);
			}
		    
			/* record the number of pres for state (*it) */
			_ts._npre[(*it) * nu + col] ++;
		    }
		}  // end transition assignment
	    }  // end input loop
	
	}  // end for loop states

	std::cout << "# of transitions: " << _ts._ntrans << '\n';
	std::cout << "length of _idpost: " << _ts._idpost.size() << '\n';
	assert(_ts._idpost.size() == _ts._ntrans);
    
	/* determine pre's by post's: loop _ts._npost and _idpost */
	_ts._idpre.resize(_ts._ntrans);  // initialize the size of pre's

	/* assign _ptrpre */
	size_t sum = 0;
	for (size_t row = 0; row < nx; ++row) {
	    for (size_t col = 0; col < nu; ++col) {
		_ts._ptrpre[row*nu + col] = sum;
		sum += _ts._npre[row*nu + col];
	    }
	}
	/* assign _idpre */
	std::vector<size_t> precount(_ptrsys->_ugrid._nv*_x._nv, 0);
	size_t idtspre;
	for (size_t row = 0; row < nx; ++row) {
	    for (size_t col = 0; col < nu; ++col) {
		for (int ip = 0; ip < _ts._npost[row*nu + col]; ++ip) {
		
		    idtspre = _ts._idpost[_ts._ptrpost[row*nu+col] + ip]*nu + col;
		    _ts._idpre[_ts._ptrpre[idtspre] + precount[idtspre]] = row;
		    precount[idtspre] ++;
		}
	    }
	}
    
	std::cout << "length of _idpre: " << _ts._idpre.size() << '\n';
	assert(_ts._idpre.size() == _ts._ntrans);
    
	return true;
    }


} // namespace rocs


#endif
