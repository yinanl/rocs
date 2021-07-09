/**
 *  The source code of the specification-guided engine.
 *
 *  Created by Yinan Li on Jan. 03, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _csolver_h
#define _csolver_h

#include <iostream>
#include <climits>
#include <vector>

#include "definitions.h"
#include "interval_paver.h"


namespace rocs {
  
    /**
     * \brief An enum to specify GOAL and AVOIDANCE.
     */
    enum SPEC
	{
	 GOAL, /**< goal (tag=1) */
	 AVOID, /**< avoid (tag=-1) */
	 FREE  /**< free workspace (tag=0) */
	};

    /**
     * \brief Bisection type.
     */
    enum BISECT
	{
	 RELMAX, /**< i = argmax{width[i]/eps[i]} */
	 ABSMAX /**< i = argmax{width[i]} */
	};


    /**
     * \brief The specification-guided engine solver.
     *
     * The state space of the system is partitioned on-the-fly, which results in a non-uniform partition.
     * Control synthesis w.r.t. an LTL formula relies on its translation to a transition matrix _M.
     * After control synthesis, the controller _ctlr is filled.
     * Computational info is recorded in _fpiter and _timer. 
     */
    class CSolver
    {
    public:

	SPtree _ctlr; /**< A controller */
	int _xdim; /**< The state dimension */
	size_t _nu; /**< The number of control values */
	std::vector<UintSmall> _M; /**< All the transitions from the current %DBA state (q0) to others. If q0-->q1 under a proposition with the label i, then _M[i] gives the index of the %DBA state q1. */
	std::vector<ivec> _goal; /**< The goal areas (an array of intervals) */
	std::vector<ivec> _obs;  /**< The avoiding areas (an array of intervals) */
	BISECT _bstype;   /**< The bisection type */
	size_t _maxiter;  /**< The maximum number of iterations */
	double _winsize;  /**< The volume of the winning set */
	
	size_t _fpiter[3];   /**< The number of iterations: max alter depth 3 */
	double _timer;   /**< The time of solving */


	/**
	 * \brief A constructor.
	 *
	 * Set the maximum number of iterations if you want the iteration to stop earlier.
	 *
	 * @param[in] ptrsys The pointer to the system dynamics 
	 * @param[in] nProps The number of propositions (default=0)
	 * @param[in] bs The bisection type (default=ABSMAX)
	 * @param[in] maxi The maximum number of iterations (default=UINT_MAX)
	 */
	template<typename system>
	CSolver(system* ptrsys, size_t nProps=0, BISECT bs=ABSMAX, size_t maxi=UINT_MAX):
	    _xdim(ptrsys->_xdim), _nu(ptrsys->_ugrid._nv), _M(nProps),
	    _bstype(bs), _maxiter(maxi), _winsize(0),
	    _fpiter{0,0,0}, _timer(0) {
		
	    SPnode root(ptrsys->_workspace, _nu);
	    _ctlr = SPtree(&root);	/* SPtree copy assignment */
	}
  
	// /**
	//  * Construct with relative subdivision.
	//  * @param ptrsys the pointer to the system dynamics.
	//  * @param maxi the maximum number of iterations.
	//  */
	// template<typename system>
	// CSolver(system* ptrsys, int n, size_t maxi):
	//     _xdim(ptrsys->_xdim),_nu(ptrsys->_ugrid._nv),
	//     _bstype(RELMAX), _maxiter(maxi),
	//     _winsize(0),
	//     _fpiter{0,0,0}, _timer(0) {

	// 	SPnode root(ptrsys->_workspace, _nu);
	// 	_ctlr = SPtree(&root);
	//     }

	// /**
	//  * Construct with absolute subdivision and default maximum number of iterations.
	//  * @param ptrsys the pointer to the system dynamics.
	//  * @param maxi the maximum number of iterations.
	//  */
	// template<typename system>
	// CSolver(system* ptrsys):
	//     _xdim(ptrsys->_xdim), _nu(ptrsys->_ugrid._nv),
	//     _bstype(RELMAX), _maxiter(UINT_MAX), _M(0),
	//     _winsize(0),
	//     _fpiter{0,0,0}, _timer(0) {

	// 	SPnode root(ptrsys->_workspace, _nu);
	// 	_ctlr = SPtree(&root);
	//     }

	/**
	 * \brief Set the vector _M.
	 *
	 * @param[in] outedge The lvalue ref to the input transition matrix
	 */
	void set_M(const std::vector<UintSmall> &outedge);

	/**
	 * \brief Determine the bisection axis.
	 *
	 * @param[in,out] box The interval vector to be bisected
	 * @param[in] eps The given precision of the grid, i.e., the minimum size of a grid
	 */
	int bisect_axis(ivec &box, const double eps[]);

	/**
	 * \brief Assign a label to a given hyper-rectangle in the state space.
	 * 
	 * @param[in] lb The lower bound of the interval
	 * @param[in] ub The upper bound of the interval
	 * @param[in] prop The proposition of type UintSmall
	 */
	void labeling(const double lb[], const double ub[], UintSmall prop);
	
	/**
	 * \brief Assign a label to a given hyper-rectangle in the state space.
	 * 
	 * @param[in] area The interval vector
	 * @param[in] prop The proposition of type UintSmall
	 */
	void labeling(ivec &area, UintSmall prop);
	
	void init_label(SPnode *node, ivec &box, UintSmall prop);
	void refine_label(SPnode *node, ivec &box, UintSmall prop);

	/**
	 * \brief Initialize the winning set.
	 * 
	 * Assume that the _ctlr is already partitioned by the labeling function.
	 *
	 */
	void init_winset();

	
	void init_winset(ivec &area);

	/** 
	 * \brief Check if the related S-domains contain targeted areas.
	 *
	 * @param[in] sdoms A vector of the pointers to the S-domains of all %DBA nodes.
	 */
	bool targetset_in_sdoms(std::vector<SPtree*> &sdoms);
  
	/**
	 * \brief Initialize _ctlr by an interval constriant.
	 *
	 * Assign a goal or an avoid area (an interval) to _ctlr by using tags.
	 * Tagging rules are:
	 * - goal       _tag <- 1;
	 * - free       _tag <- 0;
	 * - free&goal  _tag <- 2;
	 * - avoid      _tag <- -1;
	 * - free&avoid _tag <- -2;
	 * - f, g & a   _tag <- 2;
	 *
	 * @param[in] ap A label "GOAL", "AVOID" or "FREE"
	 * @param[in] lb lower bound of the interval
	 * @param[in] ub upper bound of the interval
	 */
	void init(SPEC ap, const double lb[], const double ub[]);

	/**
	 * @see init(SPEC ap, const double lb[], const double ub[])
	 * @param[in] area The given interval
	 */
	void init(SPEC ap, ivec &area);
	
	/**
	 * \brief A subroutin of init().
	 *
	 * Mark leaf node by itag if interval matches, otherwise, keep parent's tag.
	 * @param[in] node The node to be splitted w.r.t. box
	 * @param[in] box A given constraint (an interval)
	 * @param[in] itag Tag for the constraint 1(goal), -1(avoid)
	 */
	void paver_init(SPnode *node, ivec &box, short itag);
	
	/**
	 * \brief A subroutin of init().
	 *
	 * Refine initialized node.
	 * @see paver_init().
	 */
	void init_refine(SPnode *node, ivec &box, short itag);

  
	/**
	 * \brief Initialize _ctlr by a function constraint \f$f(x)\leq 0\f$.
	 *
	 * @param[in] ap A label "GOAL", "AVOID" or "FREE"
	 * @param[in] f \f$f(x)\leq 0\f$
	 * @param[in] eps Paver precision
	 */
	void init(SPEC ap, fcst f, const double eps[]);
	void paver_init(SPtree &sp, fcst f, bool inner, short itag,
			const double eps[]);
	
	/**
	 * \brief The algorithm Set Interval Via Interval Analysis (sivia).
	 *
	 * It is performed by using recursive function calls (same as using stacks).
	 *
	 * @param[in] sp The root of subpaving
	 * @param[in] ptrnode The SPnode to be refined
	 * @param[in] cst The constraint region (an interval)
	 * @param[in] eps The paver precision
	 * @param[in] inner The indicator of inner approximation (True)
	 * @param[in] itag The tag for insiders
	 */
	void sivia(SPtree &sp, SPnode *ptrnode, ivec &cst, fcst f,
		   bool inner, short itag, const double eps[]);
	void init_refine(SPtree &sp, fcst f, bool inner, short itag,
			 const double eps[]);

	/**
	 * \brief Initialize _goal by collecting leaves with tag 1.
	 *
	 * Must be used after init() functions.
	 */
	void init_goal_area();

	/**
	 * \brief Initialize _obs by collecting leaves with tag -1.
	 *
	 * Must be used after init() functions.
	 */
	void init_avoid_area();

	/**
	 * \brief Compute the normalized length of current winning set:
	 * this->_winsize = pow(vol, 1/xdim)
	 */
	void compute_winsize();	

	/**
	 * \brief Test if a box is included in the paving.
	 *
	 * @param sp the SPtree with tags (-2, -1, 0), 1, 2.
	 * @param box an interval vector.
	 * @return 0(outside), 1(inside), 2(undetermined).
	 */
	short paver_test(SPtree&, ivec&);

	/**
	 * \brief Compute the one-step backward reachable set within the scope of the same CSolver.
	 *
	 * The set of states that can reach the target set in one step under some u.
	 * Exist u for all d.
	 *
	 * @param[in,out] l Pointers of nodes to be tested (empty on return)
	 * @param[in,out] l0 Outside nodes (updated on return)
	 * @param[in,out] l1 Inside nodes (as above)
	 * @param[in,out] l2 Undetermined nodes (as above)
	 * @param[in] fcn Function pointer to a vector field
	 * @param[in] eps Minimum paver size
	 */
	template<typename system>
	void pre_cntl(system* ptrsys,
		      std::stack<SPnode*> &, std::stack<SPnode*> &,
		      std::stack<SPnode*> &, std::stack<SPnode*> &,
		      const double evaleps[]);

	/**
	 * \brief Compute the union of one-step backward reachable set across different CSolver objects:
	 *
	 * \f$W_i = \bigcup_j a_{ij}\cap Pre(W_j)\f$
	 *  
	 *  The result (union of predecessors) is represented by the updated tags in _ctlr.
	 * @see pre_cntl()
	 *
	 * @param[in] sdoms A vector of the pointers to S-domains of all %DBA nodes
	 */
	template<typename system>
	void union_of_pres(system* ptrsys, std::vector<SPtree*> &sdoms,
			   std::stack<SPnode*> &l, std::stack<SPnode*> &l0,
			   std::stack<SPnode*> &l1, std::stack<SPnode*> &l2,
			   const double evaleps[]);
	
	/**
	 * \brief Initialize queues of leaves for computation.
	 *
	 * @param[in,out] l0 A stack of leaves with tag=0.
	 * @param[in,out] l1 A stack of leaves with tag=1.
	 * @param[in,out] l2 A stack of leaves with tag=2.
	 */
	void init_leafque(std::stack<SPnode*> &l0,
			  std::stack<SPnode*> &l1,
			  std::stack<SPnode*> &l2);

	/**
	 * \brief Core subroutine for invariance control computation.
	 *
	 * @see init_leafqueue() and invariance_control()
	 * @param d The depth of iteration (0 inner most, 2 outer most)
	 */
	template<typename system>
	void inv_compute(system* ptrsys,
			 std::stack<SPnode*> &l0,
			 std::stack<SPnode*> &l1,
			 std::stack<SPnode*> &l2,
			 int d, const double eps[]);

	/**
	 * \brief Core subroutine for reachability control computation.
	 *
	 * @see init_leafque() and reachability_control() and inv_compute()
	 */
	template<typename system>
	void reach_compute(system* ptrsys,
			   std::stack<SPnode*> &l0,
			   std::stack<SPnode*> &l1,
			   std::stack<SPnode*> &l2,
			   int d, const double eps[],
			   bool vareps=false, const double emin[]=nullptr);
  
	/**
	 * \brief Compute the maximal controlled invariant set: \f$\nu X.(\text{Pre}(X)\cap X)\f$.
	 *
	 * @param[in] eps The paver precision
	 * @return 0(empty set), 1(non-empty).
	 */
	template<typename system>
	bool invariance_control(system* ptrsys, const double eps[]);
  
	/**
	 * \brief Compute the backward reachable set: iterating \f$\mu X.(\text{Pre}(X)\cup X)\f$.
	 *
	 * @param[in] eps The absolute or relative paver precision
	 * @param[in] epsmin The minimum absolute paver size (default=nullptr)
	 * @param[in] vareps An indicator of whether using variate precision (true-yes, false-no(default))
	 * @return 0(empty set), 1(non-empty).
	 */
	template<typename system>
	bool reachability_control(system* ptrsys, const double eps[],
				  bool vareps=false, const double emin[]=nullptr);

	/**
	 * \brief Compute the backward reachable set to the maximal controlled invariant set.
	 *
	 * A subset of co-Buchi set.
	 *
	 * @param[in] ei The relative precision for invariance
	 * @param[in] er The absolute or relative precision for reachability
	 * @param[in] ermin The minimum absolute paver size (default=nullptr)
	 * @param[in] vareps An indicator of whether using variate precision (true-yes, false-no(default)).
	 * @return 0(empty set), 1(non-empty).
	 */
	template<typename system>
	bool reachstay_control(system* ptrsys,
			       const double ei[], const double er[],
			       bool vareps=false, const double ermin[]=nullptr);

	/**
	 * \brief Compute the standard coBuchi winning set: \f$\nu Y.\nu X.[\text{Pre}(Y)\cup(B\cap \text{Pre}(X))]\f$.
	 *
	 * Relative and adaptive precisions.
	 * @see reachstay_control().
	 */
	template<typename system>
	bool cobuchi(system* ptrsys,
		     const double ei[], const double er[],
		     bool vareps=false, const double ermin[]=nullptr);

	/**
	 * \brief Compute the standard Buchi winning set: \f$\nu Y.\mu X.[(B\cap\text{Pre}(Y))\cup\text{Pre}(X)]\f$.
	 *
	 * @see reachstay_control()
	 */
	template<typename system>
	bool buchi(system* ptrsys, const double eps[],
		   bool vareps=false, const double ermin[]=nullptr);

	/**
	 * \brief Print controller info to screen.
	 */
	void print_controller_info() const;

	/**
	 * \brief Print a conrol table to screen.
	 */
	void print_controller() const;

	/**
	 * \brief Write (tag=1) leaf nodes of _ctlr to a log file.
	 *
	 * @param[in] filename The log file name
	 * @param[in] iter The number of iteration
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
	short t;
	while (!l.empty()) {
	    current = l.top();
	    l.pop();
	    // ptrsys->get_reach_set(y, current->_box);
	    if (!ptrsys->get_reach_set(y, current->_box)) {
		/* If the box is less than the min width, then fail. */
		int axis = bisect_axis(current->_box, evaleps);
		if (current->_box[axis].width() < evaleps[axis]) {
		    std::cout << "CSolver::pre_cntl: Fail in computing reachable set for x = "
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
		    switch (paver_test(_ctlr, y[u])) { // interval inclusion test
		    case 0:
			current->_cntl[u] = false; break;
		    case 1:
			// if (_ctlr._root->_box.isin(y[u])) {
			//     if (t != 1)
			// 	t = 1;
			//     current->_cntl[u] = true;
			// } else {  // same as ut=0
			//     current->_cntl[u] = false;
			// }
			// break;
			if (t != 1)
			    t = 1;
			current->_cntl[u] = true; break;
		    case 2:
			if (t != 1)
			    t = 2;
			/* this line is necessary, e.g. for invariant fixed points,
			   an interval can become not controlled invariant even if 
			   it is controlled invariant for the previous iterations. */
			current->_cntl[u] = false; break;
		    default:
			std::cout << "CSolver::pre_cntl: paver_test returns a wrong tag.\n";
			exit(EXIT_FAILURE); break;
		    }
		}  // end for (control update)

		/* save new tag in b0 & b1 */
		switch (t) {
		case 0:
		    current->_b0 = true;
		    current->_b1 = false;
		    l0.push(current);
		    break;
		case 1:
		    current->_b0 = false;
		    current->_b1 = true;
		    l1.push(current);
		    _winsize += current->_box.volume();
		    break;
		case 2:
		    current->_b0 = true;
		    current->_b1 = true;
		    /* compute the split axis */
		    int axis = bisect_axis(current->_box, evaleps);
		    if (current->_box[axis].width() < evaleps[axis]) {
			l2.push(current);
		    } else {
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
		    break;
		}  // end of switch
	    }  // end if (get_reach_set is successful)
	    
	}  // end while (loop all nodes in l)
    }// CSolver::pre_cntl

    template<typename system>
    void CSolver::union_of_pres(system* ptrsys, std::vector<SPtree*> &sdoms,
			   std::stack<SPnode*> &l, std::stack<SPnode*> &l0,
			   std::stack<SPnode*> &l1, std::stack<SPnode*> &l2,
			   const double evaleps[]) {
	// /***** LOGGING  *****/
	// std::ofstream logger("log_dbacontrol.txt", std::ios_base::out | std::ios_base::app);
	// /***** LOGGING  *****/
	
	SPnode *current;
	std::vector<ivec> y(_nu, ivec(_xdim));
	short t;
	while (!l.empty()) {
	    current = l.top();
	    l.pop();
	    // ptrsys->get_reach_set(y, current->_box);
	    if (!ptrsys->get_reach_set(y, current->_box)) {
		/* If the box is less than the min width, then fail. */
		int axis = bisect_axis(current->_box, evaleps);
		if (current->_box[axis].width() < evaleps[axis]) {
		    std::cout << "CSolver::union_of_pres: Fail in computing reachable set for x = "
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
		    switch ( paver_test(*sdoms[_M[current->_label]], y[u]) ) {
		    case 0:
			current->_cntl[u] = false; break;
		    case 1:
			if (t != 1)
			    t = 1;
			current->_cntl[u] = true; break;
		    case 2:
			if (t != 1)
			    t = 2;
			/* this line is necessary, e.g. for invariant fixed points,
			   an interval can become not controlled invariant even if 
			   it is controlled invariant for the previous iterations. */
			current->_cntl[u] = false;
			break;
		    default:
			std::cout << "CSolver::union_of_pres: paver_test returns a wrong tag.\n";
			exit(EXIT_FAILURE); break;
		    
		    }
		}  // end for (control update)

		// /***** LOGGING  *****/
		// if (current->_label == 1)
		//     logger << current->_box << ": test on w" << _M[current->_label] << ", t = " << t << '\n';
		// /***** LOGGING  *****/
		
		/* save new tag in b0 & b1 */
		switch (t) {
		case 0:
		    current->_b0 = true;
		    current->_b1 = false;
		    l0.push(current);
		    break;
		case 1:
		    current->_b0 = false;
		    current->_b1 = true;
		    l1.push(current);
		    _winsize += current->_box.volume();
		    break;
		case 2:
		    current->_b0 = true;
		    current->_b1 = true;
		    /* compute the split axis */
		    int axis = bisect_axis(current->_box, evaleps);
		    if (current->_box[axis].width() < evaleps[axis]) {
			l2.push(current);
		    } else {
			_ctlr.expand(current, axis);
			l.push(current->_left);
			l.push(current->_right);
		    }
		    break;
		} // end of switch
	    }  // end if (get_reach_set is successful)
	    
	}  // end while (loop all nodes in l)
	// /***** LOGGING  *****/
	// logger.close();
	// /***** LOGGING  *****/
    }// CSolver::union_of_pres
    

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


    /* non-member functions */
    /**
     * \brief Control synthesis for DBA objectives.
     * 
     * The resulting winning set and each sub-controllers are recored in the corresponding %CSolver objects.
     * @param[in,out] w A vector of S-domains (CSolvers) (assumption: already initialized)
     * @param[in] ptrsys The pointer to the system
     * @param[in] sdoms The pointers to the S-domains of all %DBA states
     * @param[in] nNodes The number of %DBA states
     * @param[in] isacc A vector of binaries (size nNodes) that marks accepting states
     * @param[in] init_w A function that initializes labeling function
     * @param[in] oid An index mapping between %DBA state indices and the indices  under consideration in the full list of S-domains.
     * @param[in] e A partition precision
     */
    template<typename system, typename F>
    void dba_control(std::vector<CSolver*> &w, system* ptrsys,
		     std::vector<SPtree*> &sdoms,
		     UintSmall nNodes, boost::dynamic_bitset<> &isacc,
		     F init_w, UintSmall oid[],
		     const double e[]) {
	// /***** LOGGING  *****/
	// std::ofstream logger;
	// /***** LOGGING  *****/

	/* Initialize the S-domains */
	std::stack<SPnode*> l;
	std::vector< std::stack<SPnode*> > l0(nNodes);
	std::vector< std::stack<SPnode*> > l1(nNodes);
	std::vector< std::stack<SPnode*> > l2(nNodes);
	for (UintSmall i = 0; i < nNodes; ++i) {
	    init_w(w, i, oid); //initialize the partition 
	    sdoms[oid[i]] = &(w[i]->_ctlr);
	    if (isacc[i]) {
	    	w[i]->init_winset();
	    }
	    w[i]->init_leafque(l0[i], l1[i], l2[i]);
	}

	/* Set maximum iteration for the inner loop (i.e. the reachability loop) */
	size_t maxNoInner = 0;
	for (UintSmall i = 0; i < nNodes; ++i) {
	    if (w[i]->_maxiter > maxNoInner)
		maxNoInner = w[i]->_maxiter;
	}

	boost::dynamic_bitset<> stop(nNodes, false);
	boost::dynamic_bitset<> start(nNodes, false);
	std::vector<size_t> lold(nNodes, 0);
	size_t iter[2] = {0,0};
	std::vector<double> t(nNodes, 0);

	/* Solve */
	// const double e = 0.1;
	clock_t tb, te, cb, ce;
	tb = clock();
	bool outerfp(false), innerfp(false);
	size_t nInner;
	while (!outerfp) {
	    ++iter[1];
	    std::cout << '\n' << "Outer loop " << iter[1] << ":\n";
	    // /***** LOGGING  *****/
	    // logger.open("log_dbacontrol.txt", std::ios_base::out | std::ios_base::app);
	    // logger << "Outer loop " << iter[1] << ":\n";
	    // logger.close();
	    // /***** LOGGING  *****/
	    
	    /* Re-initialize S-domains of non-accepting nodes:
	     * if it is not the first time to call the inner loop,
	     * empty l0, l2, and l0<-l1.
	     */
	    for (rocs::UintSmall i = 0; i < nNodes; ++i) {
		if (!isacc[i] && iter[1] > 1) {
		    delete w[i];
		    init_w(w, i, oid);
		    sdoms[oid[i]] = &(w[i]->_ctlr);
		    l0[i] = std::stack<SPnode*>();
		    l1[i] = std::stack<SPnode*>();
		    l2[i] = std::stack<SPnode*>();
		    w[i]->init_leafque(l0[i], l1[i], l2[i]);
		}
		// if (!l1[i].empty() && !isacc[i] ) {
		//     l0[i] = std::stack<rocs::SPnode*>();
		//     l2[i] = std::stack<rocs::SPnode*>();
		//     swap(l1[i], l0[i]);
		//     w[i]->_ctlr.reset_tags();
		//     // std::cout << "The size of l1[" << i << "]=" << l1[i].size()
		//     // 	      << ". The size of l0[" << i << "]=" << l0[i].size() << '\n';
		// }
	    }

#ifdef VERBOSE
	    std::cout << "Inner loop:\n";
	    std::cout << "<#S-domain>, <#iter>: <#Intervals in the S-domain> <precision>\n";
#endif
	    /* Compute w1, w2, w3 based on current w0 */
	    nInner = 0;
	    while (!innerfp && nInner < maxNoInner) { // loop until non-accepting s-domains terminate
		++iter[0];
		innerfp = true;
		++nInner;
		// /***** LOGGING  *****/
		// logger.open("log_dbacontrol.txt", std::ios_base::out | std::ios_base::app);
		// logger << "Inner loop " << iter[0] << ":\n";
		// logger.close();
		// /***** LOGGING  *****/
		
		/* Check and compute the predecessors */
		for (UintSmall i = 0; i < nNodes; ++i) {
		    if (!isacc[i] ) {
			lold[i] = l1[i].size();
			if (!start[i]) {
			    if (w[i]->targetset_in_sdoms(sdoms)) {
				std::cout << iter[1] << ": start to compute w" << oid[i] << "...\n";
				start[i] = true;
			    }
			}
			if (start[i]) {
			    // /***** LOGGING  *****/
			    // logger.open("log_dbacontrol.txt", std::ios_base::out | std::ios_base::app);
			    // logger << "w" << i << ":\n";
			    // logger.close();
			    // /***** LOGGING  *****/

			    ++w[i]->_fpiter[0];
			    cb = clock();
		
			    if (!l2[i].empty()) {
				swap(l, l2[i]);  // l=l2, l2=empty
				w[i]->union_of_pres(ptrsys, sdoms, l, l0[i], l1[i], l2[i], e);  //l=empty, l012 fill
				// std::cout << "dba_control: computing l2 is complete.\n";
			    }
			    if (!l0[i].empty()) {
				swap(l, l0[i]);  // l=l0, l0=empty
				w[i]->union_of_pres(ptrsys, sdoms, l, l0[i], l1[i], l2[i], e);  //l=empty, l012 fill
				// std::cout << "dba_control: computing l0 is complete.\n";
			    }
			    w[i]->_ctlr.tagging(INNER);

			    ce = clock();
			    t[i] += (double)(ce - cb)/CLOCKS_PER_SEC;
			    // /***** LOGGING  *****/
			    // std::cout << "Winning set of w" << i <<":\n";
			    // w[i]->_ctlr.print_leaves(w[i]->_ctlr._root, 1);
			    // std::cout << '\n';
			    // std::cout << "Partition of w" << i <<":\n";
			    // w[i]->print_controller();
			    // /***** LOGGING  *****/
#ifdef VERBOSE
			    std::cout << 'w' << oid[i] << ", iter " << iter[0] << ": " << l1[i].size() << ", ["; // w[i]->_fpiter[0] <<
			    for (int k = 0; k < ptrsys->_xdim; ++k) {
				std::cout << e[k];
				if (k < ptrsys->_xdim-1)
				    std::cout << ',';
				else
				    std::cout << "]\n";
			    }
#endif
			}
			stop[i] = l1[i].size() <= lold[i];
			innerfp &= stop[i];
		    }// end if !acc[i]
		}// end for loop all non-accepting states
		
	    }// end inner while
	    std::cout << "Outer itration " << iter[1] << ": " << iter[0] << " inner iterations.\n";
	    /* Reset innerloop counter and fixed-point marker */
	    iter[0] = 0;
	    innerfp = false; // forgot to reset in the first version

	    /* Modify w of accepting nodes by w of non-accepting nodes */
	    outerfp = true;
	    for (UintSmall i = 0; i < nNodes; ++i) {
		if (isacc[i]) {
		    lold[i] = l0[i].size() + l2[i].size();
		    if (!start[i]) {
			if (w[i]->targetset_in_sdoms(sdoms)) {
			    std::cout << ": start to compute w" << oid[i] << "...\n";
			    start[i] = true;
			}
		    }
		    if (start[i]) {
			// /***** LOGGING  *****/
			// logger.open("log_dbacontrol.txt", std::ios_base::out | std::ios_base::app);
			// logger << "w" << i << ":\n";
			// logger.close();
			// /***** LOGGING  *****/
			++w[i]->_fpiter[0];
			cb = clock();
	    
			if (!l1[i].empty()) {
			    swap(l, l1[i]);  // l=l2, l2=empty
			    w[i]->union_of_pres(ptrsys, sdoms, l, l0[i], l1[i], l2[i], e);  //l=empty, l012 fill
			}
			w[i]->_ctlr.tagging(INNER);

			ce = clock();
			t[i] += (double)(ce - cb)/CLOCKS_PER_SEC;
		    }
		    stop[i] = (l0[i].size() + l2[i].size()) <= lold[i];
		    outerfp &= stop[i];
#ifdef VERBOSE
		    std::cout << 'w' << oid[i] << ", iter " << w[i]->_fpiter[0] << ": " << l1[i].size() << ", [";
		    for (int k = 0; k < ptrsys->_xdim; ++k) {
			std::cout << e[k];
			if (k < ptrsys->_xdim-1)
			    std::cout << ',';
			else
			    std::cout << "]\n";
		    }
#endif
		}// end if acc[i]
	    }// end for loop all accepting states
	} //end outer while
    
	te = clock();
	std::cout << "Control synthesis stops." << std::endl;

	for (UintSmall i = 0; i < nNodes; ++i) {
	    w[i]->_timer = t[i];
	}
	std::cout << "Total Number of outer iterations: " << iter[1] << std::endl;
	std::cout << "Total time for control synthesis: " << (double)(te - tb)/CLOCKS_PER_SEC << std::endl;
    
    } //dba_control

    

} // namespace rocs

#endif
