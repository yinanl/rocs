/**
 * interval_paver.h 
 *
 * A class of interval paver.
 *
 * Created by Yinan Li on Aug. 16, 2016.
 *
 * Hybrid Systems Group, University of Waterloo.
 */


#ifndef _interval_paver_h_
#define _interval_paver_h_


#include "intervals/interval.h"
#include "intervals/interval_vector.h"
#include <queue>
#include <stack>
#include <boost/dynamic_bitset.hpp>
#include <mat.h>
#include <matrix.h>
#include <fstream>


namespace rocs {
  
/**
 * An enum defines set approximation type.
 */
enum APPROXTYPE {

  INNER, OUTER, EXACT
};

/**
 * A struct for subpaver node.
 */
struct SPnode
{
  ivec _box;			/**< interval vector */

  boost::dynamic_bitset<> _cntl; /**< control indicator */

  short _tag;			/**< marking in or out */
  
  bool _b0, _b1;		/**< outside(b0=1), inside(b1=1) */
  
  SPnode *_left, *_right;	/**< pointers to left & right childern */

  int _split;			/**< splitting dimension */
  

  /**
   * constructors 
   */
  SPnode() {}
  SPnode(const ivec &b, const int cn) {

    _box = b;
    _cntl.resize(cn, false);
    
    _tag = 0;			/* free workspace */
    _b0 = true;
    _b1 = false;
    
    _left =  NULL;
    _right = NULL;
    _split = -1;		/* no splitting */
  }
  
  SPnode(const ivec &b, const boost::dynamic_bitset<> &c,
       const short t, const bool b0, const bool b1) : 
  _box(b), _cntl(c), _tag(t), _b0(b0), _b1(b1), _left(NULL), _right(NULL) {

    _split = -1;
  } 

  /**
   * copy constructor
   */
  SPnode(const SPnode &other) :
  _box(other._box), _cntl(other._cntl),
  _tag(other._tag), _b0(other._b0), _b1(other._b1),
  _left(other._left), _right(other._right), _split(other._split) {}


  bool isempty() const { return _box.isempty(); }
  
};

/**
 * Print SPnode 
 */
std::ostream& operator<<(std::ostream&, const SPnode&);


/**
 * A subpaving tree class.
 */
class SPtree
{
 public:

  SPnode *_root;		/**< the root of the SPtree */
  
  /**
   * Constructors
   * @param root the pointer to the node to be copied.
   */
  SPtree(): _root(NULL) {}
  SPtree(const SPnode *root) { _root = new SPnode(*root); }

  /**
   * Copy constructor (recursive)
   * @param other the SPtree to be copied.
   */
  SPtree(const SPtree &r);
  SPnode* copyHelper(const SPnode*);

  /**
   * Destructor
   */
  ~SPtree() { release(_root);}

  /**
   * Copy assignment
   */
  SPtree& operator=(const SPtree &other);

  /**
   * Check if a node is a leaf.
   * @param ptrn the pointer to the node.
   * @return 1(is leaf), 0(otherwise).
   */
  bool isleaf(SPnode *ptrn) const;
  /**
   * Check if a node is a leaf.
   * @return 1(is empty), 0(otherwise).
   */
  bool isempty() const { return _root == NULL ? true : false;}


  /**
   * Release node memory (recursively).
   * @param ptrn pointer to the subtree to be released.
   */
  void release(SPnode *ptrn);
  
  /**
   * Expand a leaf by bisection.
   * @param leaf the leaf to be expanded by bisection.
   * @param axis the splitting axis.
   */
  void expand(SPnode *leaf, int axis);
  /**
   * Expand a leaf by connection.
   * @param leaf the leaf to be expanded by bisection.
   * @param lcld left child.
   * @param rcld right child.
   */
  void expand(SPnode *leaf, SPnode *lcld, SPnode *rcld);

  /**
   * Refine a node by another SPtree target. (not used)
   * Assumption: \f$x\in\f$ sp, \f$x\f$ is a leaf of this tree
   * @param node pointer to the node.
   * @param ptrn pointer to the node to be refined.
   * @param ptrsp pointer to the SPtree.
   * @param tag 
   */
  void refine_leaf(SPnode *ptrn, SPnode *ptrsp, short tag);

  /**
   * Retract the tree.
   * Traverse from bottom to top and merge the leaves tagged same.
   */
  void retract();
  
  /**
   * Tagging the tree (update _tag by _b0 and _b1)
   * - INNER & OUTER : leaves are either 0, 1, or -1
   * - EXACT: leaves can be 0, 1, 2, and -1 (unchanged)
   *
   * Tagging rules (all nodes):
   * - left._tag=0 & right._tag=0 -> node._tag=0
   * - left._tag=1 & right._tag=1 -> node._tag=1
   * - others                     -> node._tag=2
   *
   * Inner approximation (leaves only):
   * - b0 = 0 & b1 = 1 -> node._tag = 1
   * - others          -> node._tag = 0
   *
   * Outer approximation (leaves only):
   * - b1 = 1          -> node._tag = 1
   * - others          -> node._tag = 0
   *
   * @param approx approximation method (INNER, OUTER or EXACT).
   */
  void tagging(APPROXTYPE approx);
  void tagging();		/* to be discarded */

  /**
   * Tag negation (1-->0, 0-->1).
   */
  void negate();

  /**
   * Count leaves in preorder (depth-first).
   * @param ptrn pointer to the node from which the counted leaves branch.
   * @return the number of leaves under the node.
   */
  size_t leafcount(SPnode *ptrn) const;
  
  /**
   * Count all nodes in preorder (depth-first).
   * @param ptrn pointer to the node from which the counted nodes branch.
   * @return the number of nodes under the given node.
   */
  size_t nodecount(SPnode *ptrn) const;

  /**
   * Calculate the height of a node in a tree
   * - leaf height is 1.
   * - number of nodes of a balanced BTree: 2^H - 1.
   * @param ptrn pointer to the node.
   * @return the node height.
   */
  size_t height(SPnode *ptrn) const;

  /**
   * All leaves under a node.
   * @param ptrn pointer to the node from which the leaves branch.
   * @return a list of leaves.
   */
  std::vector<SPnode*> leaves(SPnode *ptrn) const;

  /**
   * Print SPtree level-by-level (Breadth-first)
   * @param ptrn pointer to the node.
   */
  void print_subtree(SPnode *ptrn) const;
  /**
   * Print leaves in preorder (Depth-first).
   * @param ptrn pointer to the node.
   */
  void print_leaves(SPnode *ptrn) const;
  /**
   * Print leaves with a specified tag in preorder (Depth-first).
   * @param ptrn pointer to the node.
   * @param tag the given tag.
   */
  void print_leaves(SPnode *ptrn, short tag) const;

  /**
   * Write leaves under a node to a .mat file (Depth-first).
   * @param ptrn pointer to the given node.
   * @param filename the file name.
   */
  void write2mat_leaves(SPnode *ptrn, const char *filename);
  /**
   * Write leaves under a node to a .txt file (bigger).
   * @see write2mat_leaves().
   */
  void write2txt_leaves(SPnode *ptrn, const char *filename);
};


} // namespace rocs

#endif
