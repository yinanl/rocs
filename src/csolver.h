/**
 *  csolver.h
 *
 *  An interval based control problem solver class.
 *
 *  Created by Yinan Li on Jan. 03, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _csolver_h
#define _csolver_h


#include <climits>
#include <vector>

#include "definitions.h"
#include "interval_paver.h"


namespace rocs {
  
    /**
     * An enum used to specify GOAL and AVOIDANCE.
     */
    enum SPEC
	{
	 GOAL,  // tag = 1
	 AVOID,   // tag = -1
	 FREE  // tag = 0
	};

    /**
     * Bisection type.
     */
    enum BISECT
	{
	 RELMAX, /* i = argmax{width[i]/eps[i]} */
	 ABSMAX /* i = argmax{width[i]} */
	};


    /**
     * A solver of control problems.
     */
    class CSolver
    {
    public:

	SPtree _ctlr;	     /**< A controller. */
	int _xdim;  /**< The state dimension. */
	size_t _nu;  /**< The number of control values. */
	std::vector<ivec> _goal; /**< The goal areas (an array of intervals). */
	std::vector<ivec> _obs;  /**< The avoiding areas (an array of intervals). */
	BISECT _bstype;   /**< The bisection type. */
	size_t _maxiter;  /**< The maximum number of iterations. */
	double _winsize;  /**< The volume of the winning set. */
	
	size_t _fpiter[3];   /**< The number of iterations: max alter depth 3. */
	double _timer;   /**< The time of solving. */


	/**
	 * A constructor.
	 * @param ptrsys the pointer to the system dynamics.
	 * @param bs the bisection type.
	 * @param maxi the maximum number of iterations.
	 */
	template<typename system>
	CSolver(system* ptrsys, BISECT bs, size_t maxi=UINT_MAX):
	    _xdim(ptrsys->_xdim),_nu(ptrsys->_ugrid._nv),
	    _bstype(bs),_maxiter(maxi), 
	    _winsize(0),
	    _fpiter{0,0,0}, _timer(0) {
		
	    SPnode root(ptrsys->_workspace, _nu);
	    _ctlr = SPtree(&root);	/* SPtree copy assignment */
	}
  
	/**
	 * Construct with absolute subdivision.
	 * @param ptrsys the pointer to the system dynamics.
	 * @param maxi the maximum number of iterations.
	 */
	template<typename system>
	CSolver(system* ptrsys, size_t maxi):
	    _xdim(ptrsys->_xdim),_nu(ptrsys->_ugrid._nv),
	    _bstype(ABSMAX), _maxiter(maxi),
	    _winsize(0),
	    _fpiter{0,0,0}, _timer(0) {

		SPnode root(ptrsys->_workspace, _nu);
		_ctlr = SPtree(&root);
	    }

	/**
	 * Construct with absolute subdivision and default maximum number of iterations.
	 * @param ptrsys the pointer to the system dynamics.
	 * @param maxi the maximum number of iterations.
	 */
	template<typename system>
	CSolver(system* ptrsys):
	    _xdim(ptrsys->_xdim), _nu(ptrsys->_ugrid._nv),
	    _bstype(ABSMAX), _maxiter(UINT_MAX),
	    _winsize(0),
	    _fpiter{0,0,0}, _timer(0) {

		SPnode root(ptrsys->_workspace, _nu);
		_ctlr = SPtree(&root);
	    }

	/**
	 * Determine the bisection axis.
	 * @param box the interval vector to be bisected.
	 */
	int bisect_axis(ivec &box, const double eps[]);
  
  
	/* controller initialization */
	/**
	 * Initialize _ctlr by an interval constriant.
	 * Assign a goal or an avoid area (an interval).
	 * @param ap atomic proposition "GOAL" or "AVOID".
	 * @param lb lower bound of the interval.
	 * @param ub upper bound.
	 *
	 * Tagging rules in initialization:
	 * - goal       _tag <- 1;
	 * - free       _tag <- 0;
	 * - free&goal  _tag <- 2;
	 * - avoid      _tag <- -1;
	 * - free&avoid _tag <- -2;
	 * - f, g & a   _tag <- 2;
	 */
	void init(SPEC ap, const double lb[], const double ub[]);
	void init(SPEC ap, ivec &area);
	/**
	 * Initialize by interval constraint.
	 * Mark leaf node by itag if interval matches, otherwise, keep parent's tag.
	 * @param node the node to be splitted w.r.t. box
	 * @param box a given constraint (an interval)
	 * @param itag tag for the constraint 1(goal), -1(avoid)
	 */
	void paver_init(SPnode *node, ivec &box, short itag);
	/**
	 * Refine initialized node.
	 * @see paver_init().
	 */
	void init_refine(SPnode *node, ivec &box, short itag);

  
	/**
	 * Initialize _ctlr by a function constraint f(x)<=0.
	 * @param ap approximation type: inner or outer.
	 * @param f \f$f(x)\leq 0\f$.
	 * @param eps paver precision.
	 */
	void init(SPEC ap, fcst f, const double eps[]);
	void paver_init(SPtree &sp, fcst f, bool inner, short itag,
			const double eps[]);
	/**
	 * CSP w.r.t. a single interval
	 * recursive function call (same as using stacks).
	 * @param sp root of subpaving.
	 * @param ptrnode SPnodes to be refined.
	 * @param cst constraint region (an interval).
	 * @param eps paver precision.
	 * @param inner indicator of inner approximation.
	 * @param itag tag for insiders.
	 */
	void sivia(SPtree &sp, SPnode *ptrnode, ivec &cst, fcst f,
		   bool inner, short itag, const double eps[]);
	void init_refine(SPtree &sp, fcst f, bool inner, short itag,
			 const double eps[]);

	/**
	 * Initialize _goal by collecting leaves with tag 1.
	 * Must be used after init() functions.
	 */
	void init_goal_area();

	/**
	 * Initialize _obs by collecting leaves with tag -1.
	 * Must be used after init() functions.
	 */
	void init_avoid_area();

	/**
	 * Compute the normalized length of current winning set:
	 * this->_winsize = pow(vol, 1/xdim)
	 */
	void compute_winsize();
	

	/**
	 * Test if a box is included in the paving.
	 * @param sp the SPtree with tags (-2, -1, 0), 1, 2.
	 * @param box an interval vector.
	 * @return 0(outside), 1(inside), 2(undetermined).
	 */
	short paver_test(SPtree&, ivec&);

	/**
	 * One-step backward reachable set (tag updates):
	 * The set of states that can reach the target set in one step under some u.
	 * Exist u for all d.
	 *
	 * @param l [inout] pointers of nodes to be tested (empty on return).
	 * @param l0 [inout] outside nodes (updated on return).
	 * @param l1 [inout] inside nodes (as above).
	 * @param l2 [inout] undetermined nodes (as above).
	 * @param fcn [in]function pointer to a vector field.
	 * @param eps [in] minimum paver size.
	 */
	template<typename system>
	void pre_cntl(system* ptrsys,
		      std::stack<SPnode*> &, std::stack<SPnode*> &,
		      std::stack<SPnode*> &, std::stack<SPnode*> &,
		      const double evaleps[]);
	
	/* fixed-point algorithms */
	/**
	 * Initialize queues of leaves for computation
	 * @l0 a stack of leaves with tag=0.
	 * @l1 a stack of leaves with tag=1.
	 * @l2 a stack of leaves with tag=2.
	 */
	void init_leafque(std::stack<SPnode*> &l0,
			  std::stack<SPnode*> &l1,
			  std::stack<SPnode*> &l2);

	/**
	 * Core subroutine for invariance control computation.
	 * @param d the depth of iteration (0 inner most, 2 outer most).
	 * @see init_leafqueue() and invariance_control().
	 */
	template<typename system>
	void inv_compute(system* ptrsys,
			 std::stack<SPnode*> &l0,
			 std::stack<SPnode*> &l1,
			 std::stack<SPnode*> &l2,
			 int d, const double eps[]);

	/**
	 * Core subroutine for reachability control computation.
	 * @see init_leafque() and reachability_control() and inv_compute().
	 */
	template<typename system>
	void reach_compute(system* ptrsys,
			   std::stack<SPnode*> &l0,
			   std::stack<SPnode*> &l1,
			   std::stack<SPnode*> &l2,
			   int d, const double eps[],
			   bool vareps=false, const double emin[]=nullptr);
  
	/**
	 * Maximal controlled invariant sets: \f$\nu X.(\text{Pre}(X)\cap X)\f$.
	 * @param eps absolute paver precision.
	 * @return 0(empty set), 1(non-empty).
	 */
	template<typename system>
	bool invariance_control(system* ptrsys, const double eps[]);
  
	/**
	 * Backward reachable sets: iterating \f$\mu X.(\text{Pre}(X)\cup X)\f$.
	 * @param eps absolute or relative paver precision.
	 * @param epsmin minimum absolute paver size (default=0.001).
	 * @param vareps using variate precision (true-yes, false-no).
	 * @return 0(empty set), 1(non-empty).
	 */
	template<typename system>
	bool reachability_control(system* ptrsys, const double eps[],
				  bool vareps=false, const double emin[]=nullptr);

	/**
	 * Backward reachable set to the maximal controlled invariant set.
	 *
	 * A subset of co-Buchi set. 
	 * @param ei relative precision for invariance.
	 * @param er absolute or relative precision for reachability.
	 * @param ermin minimum absolute paver size (default=0.001).
	 * @param vareps using variate precision (true-yes, false-no).
	 * @return 0(empty set), 1(non-empty).
	 */
	template<typename system>
	bool reachstay_control(system* ptrsys,
			       const double ei[], const double er[],
			       bool vareps=false, const double ermin[]=nullptr);

	/**
	 * Standard coBuchi winning set: \f$\nu Y.\nu X.[\text{Pre}(Y)\cup(B\cap \text{Pre}(X))]\f$.
	 *
	 * Relative and adaptive precisions.
	 * @see reachstay_control().
	 */
	template<typename system>
	bool cobuchi(system* ptrsys,
		     const double ei[], const double er[],
		     bool vareps=false, const double ermin[]=nullptr);


	/**
	 * Standard Buchi winning set: \f$\nu Y.\mu X.[(B\cap\text{Pre}(Y))\cup\text{Pre}(X)]\f$.
	 * @param er absolute or relative precision for reachability.
	 * @param ermin minimum absolute paver size (default=0.001).
	 * @param vareps using variate precision (true-yes, false-no).
	 * @return 0(empty set), 1(non-empty set).
	 */
	template<typename system>
	bool buchi(system* ptrsys, const double eps[],
		   bool vareps=false, const double ermin[]=nullptr);

  

	/* display and save */
	/**
	 * Print controller info to screen.
	 */
	void print_controller_info() const;

	/**
	 * Print a conrol table to screen.
	 */
	void print_controller() const;

	/**
	 * Write (tag=1) leaf nodes of _ctlr to a log file.
	 * @param filename log file name.
	 * @param iter the number of iteration.
	 */
	void log_iterations(const char* filename, int iter);
  
    };


    template<typename system>
    void CSolver::pre_cntl(system* ptrsys,
			   std::stack<SPnode*> &l, std::stack<SPnode*> &l0,
			   std::stack<SPnode*> &l1, std::stack<SPnode*> &l2,
			   const double evaleps[]) {
	SPnode *current;
	std::vector<ivec> y(_nu, ivec(_xdim));
	short t, ut;
	while (!l.empty()) {
	    current = l.top();
	    l.pop();
	    // ptrsys->get_reach_set(y, current->_box);
	    if (!ptrsys->get_reach_set(y, current->_box)) {
		/* If the box is less than the min width, then fail. */
		int axis = bisect_axis(current->_box, evaleps);
		if (current->_box[axis].width() < evaleps[axis]) {
		    std::cout << "CSolver: Fail in computing reachable set for x = "
			      << current->_box << '\n';
		    exit(EXIT_FAILURE);
		} else { /* put the current box into l0 */
		    for (size_t u = 0; u < y.size(); ++ u ) {
			current->_cntl[u] = false;
		    }
		    current->_b0 = true;
		    current->_b1 = false;
		    l0.push(current);
		}
	    } else {
		/* update control info */
		t = 0;	
		for (size_t u = 0; u < y.size(); ++ u ) {
		    ut = paver_test(_ctlr, y[u]); // interval inclusion test
		    if (ut == 1) {
			if (_ctlr._root->_box.isin(y[u])) {
			    if (t != 1)
				t = 1;
			    current->_cntl[u] = true;
			} else {  // same as ut=0
			    current->_cntl[u] = false;
			}
		    } else if (ut == 0) {
			current->_cntl[u] = false;
		    } else {
			if (t != 1)
			    t = 2;
			/* this line is necessary, e.g. for invariant fixed points,
			   an interval can become not controlled invariant even if 
			   it is controlled invariant for the previous iterations. */
			current->_cntl[u] = false;
		    }
		}  // end for (control update)

		/* save new tag in b0 & b1 */
		if (t == 0) {
		    current->_b0 = true;
		    current->_b1 = false;
		    l0.push(current);
		}
		else if (t == 1) {
		    current->_b0 = false;
		    current->_b1 = true;
		    l1.push(current);
		    _winsize += current->_box.volume();
		}
		else {
		    current->_b0 = true;
		    current->_b1 = true;
		
		    /* compute the split axis */
		    int axis = bisect_axis(current->_box, evaleps);
		    if (current->_box[axis].width() < evaleps[axis]) {
			l2.push(current);
		    }
		    else {
			_ctlr.expand(current, axis);
			l.push(current->_left);
			l.push(current->_right);
		    }
		    // double ri, r = 0;
		    // int axis = 0;
		    // for (int i = 0; i < _xdim; ++i) {
		    //     ri = current->_box[i].width()/evaleps[i];
		    //     if (r < ri) {
		    // 	axis = i;
		    // 	r = ri;
		    //     }
		    // }
		    // if (r < 1)
		    //     l2.push(current);
		    // else {
		    //     _ctlr.expand(current, axis);
		    //     l.push(current->_left);
		    //     l.push(current->_right);
		    // }
		
		}  // end if (save new tag in b0 & b1)
	    }  // end if (get_reach_set is successful)
	    
	}  // end while (loop all nodes in l)
    }


    template<typename system>
    void CSolver::inv_compute(system* ptrsys,
			      std::stack<SPnode*> &l0,
			      std::stack<SPnode*> &l1,
			      std::stack<SPnode*> &l2,
			      int d, const double eps[]) {    
	std::stack<SPnode*> l;
	size_t lold;
	bool stop = false;

#ifdef VERBOSE
	std::cout << "inv_compute:: <#iter>:<# of intervals in the winset>,<precision>\n";
#endif
	
	while (!stop) {      	
	    ++_fpiter[d];
	    lold = l0.size() + l2.size();
	    if (!l1.empty()) {
		swap(l, l1);
		pre_cntl(ptrsys, l, l0, l1, l2, eps);
	    }

	    stop = (l0.size() + l2.size()) <= lold;
	    _ctlr.tagging(INNER);  //update the tags
	    
#ifdef VERBOSE
	    std::cout << _fpiter[d] << ": " << l1.size() << ", [";
	    for (int i = 0; i < _xdim; ++i) {
		std::cout << eps[i];
		if (i<_xdim-1)
		    std::cout << ',';
		else
		    std::cout << "]\n";
	    }
#endif
	    
	} // end while
    
    }// end inv_compute


    template<typename system>
    void CSolver::reach_compute(system* ptrsys,
				std::stack<SPnode*> &l0,
				std::stack<SPnode*> &l1,
				std::stack<SPnode*> &l2,
				int d, const double er[],
				bool vareps, const double ermin[]) {
	double *eps = new double[_xdim];
	for (int i = 0; i < _xdim; ++i)
	    eps[i] = er[i];
	
	std::stack<SPnode*> l;
	size_t lold;
	bool stop = false;
	double v = _winsize;
#ifdef VERBOSE
	std::cout << "reach_compute:: <#iter>:<# of intervals in the winset>,<precision>\n";
#endif
	while (!stop && _fpiter[d] < _maxiter) {
	    ++_fpiter[d];
	    lold = l1.size();
	    if (!l2.empty()) {
		swap(l, l2);  // l=l2, l2=empty
		pre_cntl(ptrsys, l, l0, l1, l2, eps);  //l=empty, l012 fill
	    }
	    if (!l0.empty()) {
		swap(l, l0);  // l=l0, l0=empty
		pre_cntl(ptrsys, l, l0, l1, l2, eps);  //l=empty, l012 fill
	    }
	    stop = l1.size() <= lold;
	    _ctlr.tagging(INNER);
	    
#ifdef VERBOSE
	    std::cout << _fpiter[d] << ": " << l1.size() << ", [";
	    for (int i = 0; i < _xdim; ++i) {
		std::cout << eps[i];
		if (i<_xdim-1)
		    std::cout << ',';
		else
		    std::cout << "]\n";
	    }
#endif
	    if (vareps) {  // using adaptive precision
	    	if (stop) {
		    for (int i = 0; i != _xdim; ++i) {
			if (eps[i] > ermin[i]) {
			    eps[i] *= 0.5;
			    stop = false;
			}
		    }
	    	} else {
		    if (_winsize > 2 * v) {
			for (int i = 0; i != _xdim; ++i)
			    eps[i] *= 2;
			v = _winsize;
		    }
	    	} // endif
	    } // endif
	} // endwhile

	delete[] eps;
    }// end reach_compute


    template<typename system>
    bool CSolver::invariance_control(system* ptrsys, const double eps[]) {
	std::stack<SPnode*> l0, l1, l2;
	init_leafque(l0, l1, l2);

	clock_t tb, te;
	tb = clock();
    
	inv_compute(ptrsys, l0, l1, l2, 0, eps);
    
	te = clock();
	_timer = (float)(te - tb)/CLOCKS_PER_SEC;

	if (l1.empty())
	    return false;
	else
	    return true;
    }


    template<typename system>
    bool CSolver::reachability_control(system* ptrsys, const double eps[],
				       bool vareps, const double epsmin[]) {
	std::stack<SPnode*> l0, l1, l2;
	init_leafque(l0, l1, l2);
    
	clock_t tb, te;
	tb = clock();
    
	reach_compute(ptrsys, l0, l1, l2, 0, eps, vareps, epsmin);  // _fpiter[0]
    
	te = clock();
	_timer = (float)(te - tb)/CLOCKS_PER_SEC;
    
	if (l1.empty())
	    return false;
	else
	    return true;
    }


    template<typename system>
    bool CSolver::reachstay_control(system* ptrsys,
				    const double ei[], const double er[],
				    bool vareps, const double ermin[]) {
	double t;
	std::cout << "Start invariance control..." << '\n';
	if ( invariance_control(ptrsys, ei) ) {

	    std::cout << "Start reachability control..." << '\n';
	
	    t = _timer;
	    bool r = reachability_control(ptrsys, er, vareps, ermin);
	    _timer += t;
	    return r;
	
	} else {
	    return false;
	}
    }

    template<typename system>
    bool CSolver::cobuchi(system* ptrsys,
			  const double ei[], const double er[],
			  bool vareps, const double ermin[]) {
	/* do the iterations */
	double *eps = new double[_xdim];
	for (int i = 0; i < _xdim; ++i)
	    eps[i] = er[i];
	
	double v = _winsize;
	bool stop = false;
	size_t Gold1;
    
	clock_t tb, te;
	tb = clock();

	std::stack<SPnode*> G10, G11, G12;
	std::stack<SPnode*> G20, G21, G22;
	std::stack<SPnode*> l, ll0, ll2;
	init_leafque(G10, G21, G12);
	
#ifdef VERBOSE
	std::cout << "cobuchi:: <outer iter>:<# of inner iters>,<current precision parameter>\n";
#endif
	/* outer mu loop */
	while (!stop && _fpiter[1] < _maxiter) {
	    ++ _fpiter[1];
	    // std::cout << _fpiter[1] << ": ";
	    
	    /* inner nu loop */
	    if (!G21.empty()) {
		if (!(G20.empty() && G22.empty()))
		    _ctlr.tagging(INNER);  // mark (Z U G) to be 1
		
		inv_compute(ptrsys, G20, G21, G22, 0, ei);
		
		/* G2<- G20 U G22, G2 starts from empty */
		std::stack<SPnode*> G2;
		while (!G20.empty()) {
		    // G20.top()->_tag = 1;
		    G20.top()->_b0 = false;
		    G20.top()->_b1 = true;
		    G2.push(G20.top());
		    G20.pop();
		}
		while (!G22.empty()) {
		    // G22.top()->_tag = 1;
		    G22.top()->_b0 = false;
		    G22.top()->_b1 = true;
		    G2.push(G22.top());
		    G22.pop();
		}

		// std::cout << "(G2 size: " << G2.size() << ")\n";
		swap(G2, G21);
	    }
	    
	    // std::cout << _fpiter[0] << ", ";
	    
	    /* compute pre(Y) /\ G1 */
	    Gold1 = G11.size();
	    /* G1<- G10 U G12 */
	    if (!G12.empty()) {
		swap(l, G12);
		// pre_cntl(ptrsys, l, G10, G11, G12, eps);
		pre_cntl(ptrsys, l, ll0, G11, ll2, eps);
	    }
	    if (!G10.empty()) {
		swap(l, G10);
		// pre_cntl(ptrsys, l, G10, G11, G12, eps);
		pre_cntl(ptrsys, l, ll0, G11, ll2, eps);
	    }
	    swap(ll0, G10); // G10 is empty as a result of swap(l, G10)
	    swap(ll2, G12);
	    
#ifdef VERBOSE
	    std::cout << _fpiter[1] << ": " << _fpiter[0] << ", [";
	    for (int i = 0; i < _xdim; ++i) {
		std::cout << eps[i];
		if (i<_xdim-1)
		    std::cout << ',';
		else
		    std::cout << "]\n";
	    }
#endif
	    
	    stop = G11.size() <= Gold1;
	    _ctlr.tagging(INNER);	    

	    if (vareps) {  // using adaptive precision
	    	if (stop) {
		    for (int i = 0; i != _xdim; ++i) {
			if (eps[i] > ermin[i]) {
			    eps[i] *= 0.5;
			    stop = false;
			}
		    }
	    	} else {
		    if (_winsize > 2 * v) {
			for (int i = 0; i != _xdim; ++i)
			    eps[i] *= 2;
			v = _winsize;
		    }
	    	} // endif
	    } // endif
	
	} // endwhile

	te = clock();
	_timer = (float)(te - tb)/CLOCKS_PER_SEC;

	return true;
    }
    

    template<typename system>
    bool CSolver::buchi(system* ptrsys, const double er[],
			bool vareps, const double ermin[]) {
	// double eps;
	// if (vareps)
	//     eps = er * _winsize;
	// else
	//     eps = er;
    
	/* do the iterations */
	bool stop = false;
	size_t lold;
	clock_t tb, te;
	tb = clock();
	while (!stop) {
	
	    ++_fpiter[1];
	
	    std::stack<SPnode*> l0, l1, l2;
	    std::stack<SPnode*> x, l;
	    init_leafque(l0, l1, l2);
	    swap(l1, x); // x <- B, l1 <- empty.
	
	    /* inner mu loop */
	    reach_compute(ptrsys,l0, l1, l2, 0, er, vareps, ermin); // l1 does not have B

	    lold = l0.size() + l2.size();
	    swap(l, x);
	    pre_cntl(ptrsys, l, l0, x, l2, er); // x might decrease.
	
	    stop = (l0.size() + l2.size()) <= lold;
	    if (!stop) {
		/* l1 <- 0, l2 <- 0, only x = 1 */
		while (!l1.empty()) {
		    l1.top()->_tag = 0;
		    l1.top()->_b0 = true;
		    l1.top()->_b1 = false;
		    l1.pop();
		}
		while (!l2.empty()) {
		    l2.top()->_tag = 0;
		    l2.top()->_b0 = true;
		    l2.top()->_b1 = false;
		    l2.pop();
		}

		/* reinitialization: retract controller SPtree */
		_ctlr.retract();
	    }

	} // end outer while
	te = clock();
	_timer = (float)(te - tb)/CLOCKS_PER_SEC;
	return true;
    }

} // namespace rocs


#endif
