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
#include "problem.h"


/**
 * An enum used to specify GOAL and AVOIDANCE.
 */
enum SPEC
  {
    GOAL,  // tag = 1
    AVOID   // tag = -1
  };

/**
 * Bisection type.
 */
enum BISECT
  {
    RELMAXW,
    RELMAXG,
    ABSMAX
  };

typedef ivec (*fcst)(const ivec&);


/**
 * A solver of control problems.
 */
class CSolver
{
 public:

  SPtree _ctlr;	     /**< controller */
  PtrCntlProb _cpb;  /**< the pointer to a problem */
  
  BISECT _bstype;   /**< bisection type */
  size_t _maxiter;  /**< maximum number of iterations */
  double _winsize;  /**< interval covering the current winning set */
  size_t _fpiter[3];   /**< number of iterations: max alter depth 3 */
  double _timer;   /**< time of solving */


  /**
   * Constructor
   * @param prob a control problem.
   * @param bs the bisection type.
   * @param maxi the maximum number of iterations.
   */
  CSolver(PtrCntlProb prob, BISECT bs, size_t maxi) :
  _cpb(prob), _bstype(bs), _maxiter(maxi), _fpiter{0,0,0}, _timer(0) {
    
    SPnode root(_cpb->_workspace, _cpb->_vf->_unum);

    _ctlr = SPtree(&root);	/* SPtree copy assignment */
  }
  
  /**
   * Construct with absolute subdivision.
   * @param prob a control problem.
   * @param maxi the maximum number of iterations.
   */
  CSolver(PtrCntlProb prob, size_t maxi) :
  _cpb(prob), _bstype(ABSMAX), _maxiter(maxi), _fpiter{0,0,0}, _timer(0) {
    
    SPnode root(_cpb->_workspace, _cpb->_vf->_unum);

    _ctlr = SPtree(&root);	/* SPtree copy assignment */
  }

  /**
   * Construct with absolute subdivision and default maximum number of iterations.
   * @param prob a control problem.
   * @param maxi the maximum number of iterations.
   */
  CSolver(PtrCntlProb prob) :
  _cpb(prob), _bstype(ABSMAX), _maxiter(UINT_MAX), _fpiter{0,0,0}, _timer(0) {

    SPnode root(_cpb->_workspace, _cpb->_vf->_unum);

    _ctlr = SPtree(&root);	/* SPtree copy assignment */
  }

  /**
   * Determine the bisection axis.
   * @param box the interval vector to be bisected.
   */
  int bisect_axis(ivec &box);
  
  
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
   * @see paver_init()
   */
  void init_refine(SPnode *node, ivec &box, short itag);

  
  /**
   * Initialize _ctlr by a function constraint f(x)<=0.
   * @param ap approximation type: inner or outer.
   * @param f \f$f(x)\leq 0\f$.
   * @param eps paver precision.
   */
  void init(SPEC ap, fcst f, const double eps = 0.01);
  void paver_init(SPtree &sp, fcst f, bool inner, short itag,
		  const double eps = 0.01);
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
	     bool inner, short itag, const double eps = 0.01);
  void init_refine(SPtree &sp, fcst f, bool inner, short itag,
		   const double eps = 0.01);


  /**
   * Test if a box is included in the paving.
   * @param sp the SPtree with tags (-2, -1, 0), 1, 2.
   * @param box an interval vector.
   * @return 0(outside), 1(inside), 2(undetermined).
  */
  short paver_test(SPtree&, ivec&);

  /**
   * One-step backward reachable set (tag updates).
   * @param l [inout] pointers of nodes to be tested (empty on return).
   * @param l0 [inout] outside nodes (updated on return).
   * @param l1 [inout] inside nodes (as above).
   * @param l2 [inout] undetermined nodes (as above).
   * @param fcn [in]function pointer to a vector field.
   * @param eps [in] minimum paver size.
  */
  void pre_cntl(PtrVF,
		std::stack<SPnode*> &, std::stack<SPnode*> &,
		std::stack<SPnode*> &, std::stack<SPnode*> &,
	        double evaleps = 0.01);

  /**
   * Compute the size of current winning set.
   * a measure of the intervals with tag = 1.
   */
  void compute_winsize();


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
   * @param d the depth of iteration (0 inner most, 2 outer most)
   * @see init_leafqueue() and invariance_control()
   */
  void inv_compute(std::stack<SPnode*> &l0,
		   std::stack<SPnode*> &l1,
		   std::stack<SPnode*> &l2,
		   int d,
		   const double eps, BISECT bs);

  /**
   * Core subroutine for reachability control computation.
   * @see init_leafque() and reachability_control() and inv_compute()
   */
  void reach_compute(std::stack<SPnode*> &l0,
		     std::stack<SPnode*> &l1,
		     std::stack<SPnode*> &l2,
		     int d,
		     const double eps, BISECT bs,
		     const double ermin = 0.001, bool vareps = false);
  
  /**
   * Maximal controlled invariant sets: \f$\nu X.(\text{Pre}(X)\cap X)\f$.
   * @param eps absolute paver precision.
   * @return 0(empty set), 1(non-empty).
   */
  bool invariance_control(double eps, BISECT bs);
  
  /**
   * Backward reachable sets: iterating \f$\mu X.(\text{Pre}(X)\cup X)\f$.
   * @param eps absolute or relative paver precision.
   * @param bs bisection type.
   * @param epsmin minimum absolute paver size (default=0.001).
   * @param vareps using variate precision (true-yes, false-no).
   * @return 0(empty set), 1(non-empty).
   */
  bool reachability_control(const double eps, BISECT bs,
			    const double epsmin = 0.001,
			    bool vareps = false);

  /**
   * Backward reachable set to the maximal controlled invariant set.
   *
   * A subset of co-Buchi set. 
   * @param ei relative precision for invariance.
   * @param bsi bisection type for invariance.
   * @param er absolute or relative precision for reachability.
   * @param ermin minimum absolute paver size (default=0.001).
   * @param vareps using variate precision (true-yes, false-no).
   * @return 0(empty set), 1(non-empty).
   */
  bool reach_stay(const double ei, BISECT bsi,
		  const double er, BISECT bsr,
		  const double ermin = 0.001, bool vareps = false);

  /**
   * Standard Buchi winning set: \f$\nu Y.\mu X.[(B\cap\text{Pre}(Y))\cup\text{Pre}(X)]\f$.
   * @param bs bisection type.
   * @param er absolute or relative precision for reachability.
   * @param ermin minimum absolute paver size (default=0.001).
   * @param vareps using variate precision (true-yes, false-no).
   * @return 0(empty set), 1(non-empty set)
   */
  bool buchi(BISECT bs, const double er,
	     const double ermin = 0.01, bool vareps = false);

  /**
   * Standard coBuchi winning set: \f$\nu Y.\nu X.[\text{Pre}(Y)\cup(B\cap \text{Pre}(X))]\f$.
   *
   * Relative and adaptive precisions.
   * @see reach_stay()
   */
  bool cobuchi(const double ei, BISECT bsi,
	       const double er, BISECT bsr,
	       const double ermin = 0.001, bool vareps = false);
  

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
   * Write a control table (admissible leaf boxes) to a .mat file.
   * - pavings: non-uniform intervals;
   * - ctlr: a logical matrix indicating feasible controls;
   * - tag: indicating whether or not in the winning set.
   */
  void write2mat_controller(const char* filename);
  
  /**
   * Write a control table to a .txt file.
   * @see write2mat_controller()
   */
  void write2txt_controller(const char* filename);

  /**
   * Write a control tree to a .mat file (as arrays)
   * - ctree: index-searching array (2^H-1 x [split axis(1), intervals(2dim)]).
   * - cindex: comparison array (# of nodes x intervals(2dim)).
   * - cvalue: control array (# of leaves x [index(1), # of control values]).
   * @param filename
   */
  void serialize_controller(const char* filename); /* store _ctlr in arrays (keep tree structure) */

  /**
   * Write (tag=1) leaf nodes of _ctlr to a log file.
   * @param filename log file name.
   * @param iter the number of iteration.
   */
  void log_iterations(const char* filename, int iter);
  
};




#endif
