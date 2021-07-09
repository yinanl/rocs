/**
 *  The definition of an abstraction class.
 *
 *  Created by Yinan Li on Aug. 30, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _abstraction_h_
#define _abstraction_h_

#include <armadillo>
#include <boost/dynamic_bitset.hpp>
#include <functional>
#include <set>

#include "transition.hpp"
#include "system.hpp"


namespace rocs {

    // /**
    //  * Boolean function defining a set.
    //  */
    // typedef bool (*gset)(const ivec &x);
    
    /**
     * \brief Weighting function (of edges) of transitions.
     */
    typedef double (*WGT)(std::vector<double> &x0,
    			  std::vector<double> &x1,
    			  std::vector<double> &u0);
    
    /**
     * \brief An abstraction class of a dynamical system.
     *
     * This is a template class where the typename %S will be replaced by the actual dynamical system type.
     */
    template<typename S>
    class abstraction {
	
    public:
	grid _x;  /**< A grid of states */
	S *_ptrsys; /**< A pointer to the system object */
	std::vector<int> _labels; /**< A vector of labels for the state grid _x */
	WGT _wf;  /**< Weighting callback function */
	fts _ts;  /**< The finite transition system */
  
	/**
	 * \brief A constructor.
	 *
	 * Assign user defined system dynamics to an abstraction.
	 *
	 * @param[in] sys The pointer to a system object
	 */
	abstraction(S *sys) : _x(sys->_xdim), _ptrsys(sys), _wf(NULL) {}
	/**
	 * \brief A constructor.
	 *
	 * Assign user defined system dynamics as well as a weighting function to an abstraction.
	 *
	 * @param[in] sys The pointer to a system object
	 * @param[in] wf A weighting function
	 */
	abstraction(S *sys, WGT wf) : _x(sys->_xdim), _ptrsys(sys), _wf(wf) {}

	/**
	 * \brief Initialize the states in the abstraction.
	 *
	 * This includes 
	 * - assigning a grid of states by the given bounds and grid size.
	 * - assigning labels
	 *
	 * @param[in] eta An array of grid size
	 * @param[in] xlb An array of lower bounds
	 * @param[in] xub An array of upper bounds
	 */
	void init_state(const double eta[], const double xlb[], const double xub[]) {
	    _x.gridding(eta, xlb, xub);
	    _labels.resize(_x._nv+1, 0); //the extra one is an out-of-domain node
	}

	/**
	 * \brief Initialize the finite transition system.
	 *
	 * Assign the number of predecessing (fts._npre) and successing (fts._npost) states.
	 */
	bool init_transitions() {
	    if (_x._nv > 0 && _ptrsys->_ugrid._nv > 0) {
		_ts.init(_x._nv+1, _ptrsys->_ugrid._nv);
	    } else {
		std::cout << "Transition initialization failed: gridding problem.\n";
		return false;
	    }
	    return true;
	}
	
	/**
	 * \brief Get state indicies of a given hyper-rectangle area.
	 *
	 * The area is given by its lower and upper bounds.
	 * It is allowed to be out of domain, but only the inside domain part is considered.
	 *
	 * @param xlb The lower bound of the area
	 * @param xub The upper bound of the area
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
	 * @param[in] g The function that determines if x is inside the set
	 * @return a list of indicies.
	 */
	std::vector<size_t> get_discrete_states(std::function<bool(const ivec &)> g) {
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

	/**
	 * \brief A template function that assigns labels to the state grid.
	 *
	 * @param[in] labeling A function that returns the label of a grid id
	 */
	template<typename F>
	void assign_labels(F labeling) {
	    for (size_t i = 0; i < _x._nv; ++i) {
		if(_labels[i] > -1)
		    _labels[i] = labeling(i);
	    }
	}

	/**
	 * \brief Assign the label to the out-of-domain part.
	 *
	 * @param[in] label The label given to the out-of-domain part
	 */
	void assign_label_outofdomain(int label) {
	    _labels[_x._nv] = label;
	}
	

	/**
	 * \brief Assign transitions with robustness margins e1,e2.
	 *
	 * This function constructs an fts according to the given system.
	 * 
	 * @param[in] e1 the robustness margin at the initial point.
	 * @param[in] e2 the robustness margin at the end point.
	 * @return whether construction is successful.
	 */
	bool assign_transitions(const double e1[], const double e2[]) {
	    /* Initialize _ts (the transition system) */
	    if (!init_transitions())
		return false;
	    int n = _x._dim;
	    size_t nx = _x._nv;
	    size_t nu = _ptrsys->_ugrid._nv;
	    // std::cout << "dim=" << n << ',' << "#states=" << nx << ",#inputs=" << nu <<'\n';

	    ivec ie2(n);
	    for (int j = 0; j < n; ++j)
	    	ie2.setval(j, interval(-e2[j],e2[j]));

	    std::vector<double> xmin(n), xmax(n);
	    for (int k = 0; k < n; ++k) {
	    	xmin[k] = _x._valmin[k]-_x._gw[k]/2.0;
	    	xmax[k] = xmin[k] + _x._gw[k]*_x._size[k];
	    }
	    ivec xbds(n);
	    for (int k = 0; k < n; ++k)
		xbds.setval(k, interval(xmin[k], xmax[k]));

	    /* Out-of-domain indicator:
	     * 0:inside, 1:fully outside, 2:partially outside
	     */
	    std::vector<int> out_of_domain(nx*nu,0);
	
	    ivec y0(n);
	    std::vector<ivec> yt(nu, ivec(n));
	    std::vector<double> ytl(n), ytu(n);
	    size_t r, postid;
	    double idu;
	    std::vector<size_t> il(n), iu(n);
	    std::vector<size_t> stpost(nx*nu, 0), rngpost(nx*nu*n, 0);
	    
	    /* loop state grids */
	    for (size_t row = 0; row < nx; ++row) {
		/* the avoid points (labeled by -1) has no outgoing edges */
		if (_labels[row] == -1)
		    continue;
		/* compute reachable set */
		for (int k = 0; k < n; ++k) {
		    y0.setval(k, interval(_x._data[row][k]-_x._gw[k]/2.0,
		    			  _x._data[row][k]+_x._gw[k]/2.0));
		    // y0.setval(k, interval(_x._data[row][k]-_x._gw[k]/2.0-e1[k],
		    // 			  _x._data[row][k]+_x._gw[k]/2.0+e1[k]));
		}		
		_ptrsys->get_reach_set(yt, y0);
		if (!yt.empty()) {
		    /* loop inputs */
		    for (size_t col = 0; col < nu; ++col) {
			// yt[col] += ie2;
			// /* Skip the yt[col] that is out of the domain _bds.
			//  * No need to delete the transitions to avoid area, since they are sinks.
			//  * */
			// if (!xbds.isin(yt[col]))
			//     continue;
			// /** Using this line will probably result in failure in synthesis **/
			// // if (!_x._bds.isin(yt[col]))
			// //     continue;
			
			if(xbds.isout(yt[col])) { //fully out of domain
			    out_of_domain[row*nu+col] = 1;
			    _ts._npost[row*nu+col] = 1;
			    _ts._ptrpost[row*nu+col] = _ts._ntrans;
			    ++_ts._ntrans;
			    continue;
			}
			if(!xbds.isin(yt[col])) {//partially inside domain
			    out_of_domain[row*nu+col] = 2;
			    /* Take intersection of xbds and yt[col] */
			    for(int k = 0; k < n; ++k) {
				ytl[k] = yt[col][k].getinf()<xmin[k] ? xmin[k] : yt[col][k].getinf();
				ytu[k] = yt[col][k].getsup()>xmax[k] ? xmax[k] : yt[col][k].getsup();
			    }
			} else { //fully inside domain
			    yt[col].getinf(ytl);
			    yt[col].getsup(ytu);
			}
			/* compute the indices of the bounds of the yt[col] */
			_ts._npost[row*nu + col] = 1;
			for (int k = 0; k < n; ++k) {
			    il[k] = static_cast<size_t> ((ytl[k]-_x._valmin[k])/_x._gw[k]+0.5);
			    idu = (ytu[k]-_x._valmin[k])/_x._gw[k] + 0.5;
			    iu[k] = static_cast<size_t> (idu);
			    if (std::fabs(idu - iu[k]) < EPSIVAL)
			    	iu[k] = iu[k] - 1;
			    // if (std::fabs(fmod(ytu[k]-_x._valmin[k], _x._gw[k]) - _x._gw[k]/2.) < EPSIVAL) /* fmod is slow */
			    // 	iu[k] = static_cast<size_t> ((ytu[k]-_x._valmin[k])/_x._gw[k]);
			    // else
			    // 	iu[k] = static_cast<size_t> ((ytu[k]-_x._valmin[k]+_x._gw[k]/2.0)/_x._gw[k]);
			    stpost[row*nu+col] += _x._base[k] * il[k]; /* compute the smallest post ID of all the post grid points */
			    rngpost[n*(row*nu+col)+k] = iu[k] - il[k] + 1;
			    _ts._npost[row*nu+col] *= rngpost[n*(row*nu+col)+k];
			}
			_ts._ptrpost[row*nu+col] = _ts._ntrans;
			_ts._ntrans += _ts._npost[row*nu+col];

			/* Add an out-of-domain transition */
			if(out_of_domain[row*nu+col]) {
			    ++_ts._npost[row*nu+col];
			    ++_ts._ntrans;
			}
			// /********** logging **********/
			// if(row == 0 && col == 16) {
			//     std::cout << "y0=" << y0 << '\n'
			// 	      << "yt=" << yt[col] << '\n'
			// 	      << "out_of_domain=" << out_of_domain[row*nu+col] << '\n';
			//     if(out_of_domain[row*nu+col]>1) {
			// 	std::cout << "Intersection with domain: ";
			// 	for(int k = 0; k < n; ++k)
			// 	    std::cout << '[' << ytl[k] << ',' << ytu[k] << "] ";
			// 	std::cout << '\n';
			//     }
			//     std::cout << "starting address of post transitions: "
			// 	      << _ts._ptrpost[row*nu+col] << '\n';
			//     std::cout << "# of post transitions: "
			// 	      << _ts._npost[row*nu+col] << '\n';
			// }
			// /********** logging **********/
		    } // end input loop
		} else {
		    yt.resize(nu, ivec(n));
		}// end yt empty check
	    }  //end loop state grids
	    /* Assign out-of-domain posts */
	    if(_labels[nx] >= 0) {//assign a self-loop if the out-of-domain node is not avoided
		for(size_t col = 0; col < nu; ++col) {
		    _ts._npost[nx*nu+col] = 1;
		    _ts._ptrpost[nx*nu+col] = _ts._ntrans;
		    ++_ts._ntrans;
		}
	    }
	    
	    /* assign _idpost, _cost and _npre */
	    _ts._idpost.resize(_ts._ntrans);
	    double w;
	    // for (size_t row = 0; row < nx; ++row) {
	    // 	for (size_t col = 0; col < nu; ++col) {
	    // 	    /* assign post grid point IDs to _idpost */
	    // 	    if (_ts._npost[row*nu+col] == 0)
	    // 		continue;
	    // 	    for (int l = 0; l < _ts._npost[row*nu+col]; ++l) {
	    // 		postid = stpost[row*nu+col];
	    // 		r = l;
	    // 		for (int k = 0; k < n; ++k) {
	    // 		    postid += (r % rngpost[n*(row*nu+col)+k])*_x._base[k];
	    // 		    r = r / rngpost[n*(row*nu+col)+k];
	    // 		}
	    // 		_ts._idpost[_ts._ptrpost[row*nu+col]+l] = postid;
	    // 		/* assign cost (worst case): maximum from all posts */
	    // 		if (_wf) {
	    // 		    w = (*_wf)(_x._data[row], _x._data[postid], _ptrsys->_ugrid._data[col]);
	    // 		    _ts._cost[row*nu+col] = _ts._cost[row*nu+col]<w ? w : _ts._cost[row*nu+col];
	    // 		}
			
	    // 		_ts._npre[postid*nu+col] ++;
	    // 	    }
	    // 	} /* end for col */
	    // } /* end for row */
	    int num_post;
	    size_t ptrout;
	    for (size_t row = 0; row < nx; ++row) {
		for (size_t col = 0; col < nu; ++col) {
		    // /********** logging **********/
		    // if(row == 0 && col == 16) {
		    // 	std::cout << "post nodes: ";
		    // }
		    // /********** logging **********/
		    /* assign post grid point IDs to _idpost */
		    if(out_of_domain[row*nu+col] != 1) {//intersect with domain
			if(out_of_domain[row*nu+col]) {
			    num_post = _ts._npost[row*nu+col] - 1;
			} else {
			    num_post = _ts._npost[row*nu+col];
			}
			for (int l = 0; l < num_post; ++l) {
			    postid = stpost[row*nu+col];
			    r = l;
			    for (int k = 0; k < n; ++k) {
				postid += (r % rngpost[n*(row*nu+col)+k])*_x._base[k];
				r = r / rngpost[n*(row*nu+col)+k];
			    }
			    _ts._idpost[_ts._ptrpost[row*nu+col]+l] = postid;
			    // /********** logging **********/
			    // if(row == 0 && col == 16) {
			    // 	std::cout << "ptrpost[" << _ts._ptrpost[row*nu+col]+l << "]="
			    // 		  << postid << '\n';
			    // }
			    // /********** logging **********/
			    /* assign cost (worst case): maximum from all posts */
			    if (_wf) {
				w = (*_wf)(_x._data[row], _x._data[postid], _ptrsys->_ugrid._data[col]);
				_ts._cost[row*nu+col] = _ts._cost[row*nu+col]<w ? w : _ts._cost[row*nu+col];
			    }
			
			    ++_ts._npre[postid*nu+col];
			    // /********** logging **********/
			    // if(postid == 0 && col == 16)
			    // 	std::cout << "current # of predecessor of node x=0 under u=16: "
			    // 		  << _ts._npre[postid*nu+col] << '\n';
			    // /********** logging **********/
			}
		    }
		    if(out_of_domain[row*nu+col]) {//fully or partially out of domain
			ptrout = _ts._ptrpost[row*nu+col]+_ts._npost[row*nu+col]-1;
			_ts._idpost[ptrout] = nx; //post node is xout
			++_ts._npre[nx*nu+col];
			// /********** logging **********/
			// if(row == 0 && col == 16) {
			//     std::cout << "ptrpost[" << ptrout << "]="
			// 	      << nx << '\n';
			// }
			// /********** logging **********/
		    }
		} /* end for col */
	    } /* end for row */
	    if(_labels[nx] >= 0) {//assign a self-loop if the out-of-domain node is not avoided
		for(size_t col = 0; col < nu; ++col) {
		    _ts._idpost[_ts._ptrpost[nx*nu+col]] = nx; //xout is a sink
		    ++_ts._npre[nx*nu+col];
		}
	    }
	
	    /* Determine pre's by post's: loop _npost and _idpost */
	    _ts._idpre.resize(_ts._ntrans);  // initialize the size of pre's
	    /* assign _ptrpre by _npre */
	    size_t sum = 0;
	    for (size_t row = 0; row < nx+1; ++row) {
		for (size_t col = 0; col < nu; ++col) {
		    _ts._ptrpre[row*nu+col] = sum;
		    sum += _ts._npre[row*nu+col];
		}
	    }
	    assert(sum == _ts._ntrans);
	    
	    /* assign _idpre */
	    std::vector<size_t> precount(nu*(nx+1), 0);
	    size_t idtspre;
	    for (size_t row = 0; row < nx+1; ++row) {
		for (size_t col = 0; col < nu; ++col) {
		    // /********** logging **********/
		    // if (row == 0 && col == 16) {
		    // 	std::cout << "Assign pre transitions:\n";
		    // }
		    // /********** logging **********/
		    for (int ip = 0; ip < _ts._npost[row*nu+col]; ++ip) {
			postid = _ts._idpost[_ts._ptrpost[row*nu+col]+ip];
			idtspre = postid * nu + col;
			_ts._idpre[_ts._ptrpre[idtspre]+precount[idtspre]++] = row;
			// precount[idtspre]++;
			// /********** logging **********/
			// if (row == 0 && col == 16) {
			//     std::cout << postid << ": idpre[" << _ts._ptrpre[idtspre]
			// 	      << '+' << precount[idtspre]-1 << "]="
			// 	      << row << '\n';
			// }
			// if(_ts._ptrpre[idtspre]+precount[idtspre]-1 == 46) {
			//     std::cout << "idpre[46] is filled at: "
			// 	      << row << "->(" << col << ")->" << postid
			// 	      << '\n';
			// }
			// /********** logging **********/
		    }
		}
	    }
	    return true;
	}
	
  
	/**
	 * \brief Assign transitions by subgridding.
	 * see assign_transitions()
	 *
	 * @param[in] rp[] the pointer to an array of relative subgridding size
	 * @return whether construction is successful.
	 */
	bool assign_transitions_subgridding(const double rp[]) {
	    /*********** logging ***********/
	    // std::fstream logfile;
	    // logfile.open("y.log", std::ios::out | std::ios::ate);
	    // std::fstream logfile2;
	    // logfile2.open("post946.log", std::ios::out | std::ios::ate);
	    /*********** logging ***********/
	    if (!init_transitions())
		return false;
	    int n = _x._dim;
	    size_t nx = _x._nv;
	    size_t nu = _ptrsys->_ugrid._nv;
	  
	    /* compute the number of sub grid points */
	    size_t subnv = 1;
	    std::vector<double> subgw(n);
	    std::vector<size_t> number(n);
	    for (int k = 0; k < n; ++k) {
		number[k] = ceil(1.0 / rp[k]);
		subgw[k] = _x._gw[k] / number[k];
		subnv *= number[k];
	    }
	    // std::cout << "number of subgrids: " << subnv << '\n';
    
	    /* transition computation by interval subgridding: loop states */
	    std::vector<double> xmin(n);
	    // std::vector<double> xc(n);
	    ivec v(n);
	    std::vector<size_t> subposts;
	    std::vector<std::vector<double>> sub(subnv, std::vector<double> (n));
	    std::vector<size_t>::iterator iter;
	    int np = 0;
	    /* loop state grids */
	    for (size_t row = 0; row < nx; ++row) {
		if (_labels[row] == -1) /* skip the avoid grid points (labeled by -1) */
		    continue;
		/* compute reachable set by interval subgridding */
		if (subnv == 1) {  // no subgridding
		    sub[0] = _x._data[row];
		} else {  // subnv > 1
		    for (int k = 0; k < n; ++k) {
			// xmin[k] = xc[k] - _gw[k]/2. + _rp[k]*_gw[k]/2.;
			xmin[k] = _x._data[row][k] - _x._gw[k]/2. + subgw[k]/2.;
		    }
		    _x.griddingHelper(sub, xmin, subgw, number, subnv);
		}
		std::vector< std::vector<ivec> > ys(subnv, std::vector<ivec> (nu));  // ys[vi][ui]
		for (int vi = 0; vi < subnv; ++vi) {
		    for (int k = 0; k < n; ++k) {
			v.setval(k, interval(sub[vi][k]-subgw[k]/2, sub[vi][k]+subgw[k]/2));
		    } // assign v
		    /* get a list of post intervals w.r.t. different inputs */
		    _ptrsys->get_reach_set(ys[vi], v);
		} // end for loop (subgrid)
		
		/* loop inputs */
		for (size_t col = 0; col < nu; ++col) {
		    std::set<size_t> posts;
		    /*********** logging ***********/
		    // logfile << col << ":\n";
		    /*********** logging ***********/
	    
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
			/*********** logging ***********/
			// logfile << ys[vi][col] << "(";
			// for (int i = 0; i < subposts.size(); ++i)
			//     logfile << subposts[i] << ',';
			// logfile << ")\n";
			/*********** logging ***********/
		    } // end collecting posts

		    /* assign posts */
		    if (!posts.empty()) {
			/*********** logging ***********/
			// logfile2 << col << '(' << posts.size() << "): ";
			/*********** logging ***********/
			_ts._npost[row*nu + col] = posts.size();
			_ts._ptrpost[row*nu + col] = _ts._ntrans;
			_ts._ntrans += posts.size();
		
			for (std::set<size_t>::iterator it = posts.begin(); it != posts.end(); ++it) {
			    _ts._idpost.push_back(*it);
		    
			    /* record the number of pres for state (*it) */
			    _ts._npre[(*it) * nu + col] ++;
			    /*********** logging ***********/
			    // logfile2 << *it << ", ";
			    /*********** logging ***********/
			}
			// logfile2 << '\n';
		    }  // end transition assignment
		}  // end input loop
		// /*********** logging ***********/
		// logfile << '\n';
		// /*********** logging ***********/
	
	    }  // end for loop states

	    std::cout << "# of transitions: " << _ts._ntrans << '\n';
	    std::cout << "length of _idpost: " << _ts._idpost.size() << '\n';
	    assert(_ts._idpost.size() == _ts._ntrans);

	    /* assign _cost */
	    double w;
	    size_t postid;
	    for (size_t row = 0; row < nx; ++row) {
		for (size_t col = 0; col < nu; ++col) {
		    /* assign post grid point IDs to _idpost */
		    if (_ts._npost[row*nu+col] == 0)
			continue;
		    for (int l = 0; l < _ts._npost[row*nu+col]; ++l) {
			/* assign cost (worst case): maximum from all posts */
			postid = _ts._ptrpost[row*nu+col]+l;
			if (_wf) {
			    w = (*_wf)(_x._data[row], _x._data[_ts._idpost[postid]], _ptrsys->_ugrid._data[col]);
			    _ts._cost[row*nu+col] = _ts._cost[row*nu+col]<w ? w : _ts._cost[row*nu+col];
			}
		    }
		} /* end for col */
	    } /* end for row */
    
	    /* determine pre's by post's: loop _ts._npost and _idpost */
	    _ts._idpre.resize(_ts._ntrans);  // initialize the size of pre's
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

	    // logfile.close();
	    // logfile2.close();
    
	    return true;
	}
	
    }; /* the abstraction class */

} // namespace rocs


#endif
