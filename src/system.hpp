/**
 * controlsystem.h
 *
 * A control system class.
 *
 * Created by Yinan Li on April 27, 2018.
 *
 * Hybrid Systems Group, University of Waterloo.
 */



#ifndef _system_h_
#define _system_h_

#include <stdlib.h>
#include "grid.h"
#include "interval.h"
#include "interval_vector.h"
#include "flow_taylor.hpp"


namespace rocs {

    /**
     * System class: <T,X,U>.
     * Must contain control to be used by CSolver, and U can be a singleton.
     */
    class System {
    public:
	const char *_name;  /**< problem name (can't change) */

	double _tau;
	
	int _xdim;	/**< state dimension */
	ivec _workspace;  /**< state space */

	int _udim;	/**< input dimension */
	grid _ugrid;  /**< a grid of controls */
	
	/**
	 * Constructor.
	 * @param name control problem name.
	 * @param xd dimension of state space.
	 * @param ws state space (an interval vector).
	 */
	System(const char *name, const double t, const int xd, const int ud = 1):
	    _name(name),_tau(t),_xdim(xd),_workspace(xd),_udim(ud),_ugrid(ud) {} //treat as no control if don't specify it.


	/**
	 * Initialize the workspace (X):
	 * @param ws an interval vector.
	 * @param lb[] an array of lower bound.
	 * @param ub[] an array of upper bound.
	 */
	void init_workspace(const ivec &ws) {
	    assert(ws.getdim() == _workspace.getdim());
	    _workspace = ws;
	}
	void init_workspace(const double lb[], const double ub[]) {
	    for (int i = 0; i < _xdim; ++i)
		_workspace[i] = interval(lb[i], ub[i]);
	}

	/**
	 * Initialize the input set {U}:
	 * @param mu[] an array of grid width.
	 * @param lb[] an array of lower bound.
	 * @param ub[] an array of upper bound.
	 */
	void init_inputset(const double mu[], const double lb[], const double ub[]) {
	    _ugrid.init(mu, lb, ub);
	    _ugrid.gridding();
	}
	/**
	 * Initialize the input set {U}:
	 * @param U an array of selected input values (double).
	 */
	void init_inputset(vecRn &U) {
	    _ugrid._nv = U.size();
	    _ugrid._dim = U[0].size();
	    _ugrid._data = U;
	}

	/**
	 * I/O interface
	 */
	friend std::ostream& operator<<(std::ostream&, const System&);
  
    };



    
    /**
     * Discrete-time system class with real-valued controls.
     */
    template<typename F>
    class DTCntlSys : public System {
    public:
	DTCntlSys(const char *name, const double t, const int xd, const int ud):
	    System(name, t, xd, ud) {}
	
	bool get_reach_set(std::vector<ivec> &xt, const ivec &x0) {
	    for (size_t i = 0; i < _ugrid._nv; ++i)
		F(xt[i], x0, _ugrid._data[i]);
	    return true;
	}
    
    };

    /**
     * Discrete-time switched system class with finite number of modes.
     */
    template<typename F>
    class DTSwSys : public System {
    public:
	DTSwSys(const char *name, const double t, const int xd, const int m):
	    System(name, t, xd) {_ugrid._nv = m; _ugrid._dim = 1;}
	
	bool get_reach_set(std::vector<ivec> &xt, const ivec &x0) {
	    for (size_t i = 0; i < _ugrid._nv; ++i)
		F(xt[i], x0, i+1);
	    return true;
	}
    
    };

    /**
     * Discrete-time system class without controls.
     */
    template<typename F>
    class DTSys : public System {
    public:
	DTSys(const char *name, const double t, const int xd):
	    System(name, t, xd) {_ugrid._nv = 1; _ugrid._dim = 1;}
	
	bool get_reach_set(std::vector<ivec> &xt, const ivec &x0) {      
	    F(xt[0], x0);
	    return true;
	}
    
    };
    


    /**
     * Continuous-time system class.
     */
    class ContinuousTime : public System {
    public:
	double _delta;
	params* _ptrparams;

	ContinuousTime(const char *name, const double t, const int xd, const int ud,
		       const double d, params *p):
	    System(name, t, xd, ud), _delta(d), _ptrparams(p) {}
	
	ContinuousTime(const char *name, const double t, const int xd,
		       const double d, params *p):
	    System(name, t, xd), _delta(d), _ptrparams(p) {}
    };
    

    /**
     * Continuous-time real-valued control system class.
     */
    template<typename F>
    class CTCntlSys : public ContinuousTime {
    public:
	std::vector< flowTaylor<F,interval,Rn>* > _flows;
	
	CTCntlSys(const char *name, const double t, const int xd, const int ud,
		  const double d, params *p):
	    ContinuousTime(name, t, xd, ud, d, p) {}

	
	void allocate_flows() {
	    assert(_ugrid._nv > 0);
	    _flows.resize(_ugrid._nv);
	    for (size_t i = 0; i < _flows.size(); ++i)
		_flows[i] = new flowTaylor<F, interval, Rn> (_ugrid._data[i], _ptrparams, _tau, _delta);
	}

	void release_flows() {
	    if (!_flows.empty()) {
		for (size_t i = 0; i < _flows.size(); ++i) {
		    delete _flows[i];
		}
		_flows.clear();
	    }
	}
    
	bool get_reach_set(std::vector<ivec> &xt, const ivec &x0) {
	    for (size_t i = 0; i < _ugrid._nv; ++i) {
		// _flows[i]->reachset_robust(xt[i], x0, _ptrparams->eps);
		if (!_flows[i]->reachset_robust(xt[i], x0, _ptrparams->eps)) {
		    // std::cout << "Not valid reachable set for x = " << x0 << '\n'
		    // 	      << "Current control input u = ";
		    // for (size_t j = 0; j < _ugrid._dim; ++j)
		    // 	std::cout << _ugrid._data[i][j] << ' ';
		    // std::cout << '\n';
		    // exit(EXIT_FAILURE);
		    return false;
		}
	    }
	    return true;
	}
    
    };


    /**
     * Continuous-time switched system class.
     */
    template<typename F>
    class CTSwSys : public ContinuousTime {
    public:
	std::vector< flowTaylor<F,interval,int>* > _flows;
	
	CTSwSys(const char *name, const double t, const int xd, const int m,
	      const double d, params *p):
	    ContinuousTime(name, t, xd, d, p) {_ugrid._nv = m; _ugrid._dim = 1;}
	
	
	void allocate_flows() {
	    _flows.resize(_ugrid._nv);
	    for (size_t i = 0; i < _flows.size(); ++i) // mode starts from 1.
		_flows[i] = new flowTaylor<F, interval, int> (i+1, _ptrparams, _tau, _delta);
	}

	void release_flows() {
	    if (!_flows.empty()) {
		for (size_t i = 0; i < _flows.size(); ++i) {
		    delete _flows[i];
		}
		_flows.clear();
	    }
	}
    
	bool get_reach_set(std::vector<ivec> &xt, const ivec &x0) {
	    for (size_t i = 0; i < _ugrid._nv; ++i) {
		if(!_flows[i]->reachset_robust(xt[i], x0, _ptrparams->eps)) {
		    // std::cout << "Not valid reachable set for x = " << x0 << '\n'
		    // 	      << "Current control input u = ";
		    // for (size_t j = 0; j < _ugrid._dim; ++j)
		    // 	std::cout << _ugrid._data[i][j] << ' ';
		    // std::cout << '\n';
		    // exit(EXIT_FAILURE);
		    return false;
		}
	    }
	    return true;
	}
    
    };
    

    /**
     * Continuous-time system class without controls.
     */
    template<typename F>
    class CTSys : public ContinuousTime {
    public:
	std::vector< flowTaylor<F,interval,bool>* > _flows;
	
	CTSys(const char *name, const double t, const int xd,
	      const double d, params *p):
	    ContinuousTime(name, t, xd, d, p) {_ugrid._nv = 1; _ugrid._dim = 1;}
	
	
	void allocate_flows() {
		_flows.resize(1);
		_flows[0] = new flowTaylor<F, interval, bool> (_ptrparams, _tau, _delta);
	}

	void release_flows() {
	    if (!_flows.empty()) {
		delete _flows[0];
		_flows.clear();
	    }
	}
    
	bool get_reach_set(std::vector<ivec> &xt, const ivec &x0) {
	    if(!_flows[0]->reachset_robust(xt[0], x0, _ptrparams->eps)) {
		// std::cout << "Not valid reachable set for x = " << x0 << '\n';
		// exit(EXIT_FAILURE);
		return false;
	    }
	    return true;
	}
    
    };




    /* I/O friend functions */
    std::ostream& operator<<(std::ostream &out, const System &sys) {
	out << '<' << sys._name << '>' << ":\n";
	out << "- state dimension: " << sys._xdim << '\n';
	out << "- input dimension: " << sys._udim << '\n';
  
	out << "- state space: " << '\n';
	out << sys._workspace << '\n';
  
	return out;
    }

} // namespace rocs


#endif
