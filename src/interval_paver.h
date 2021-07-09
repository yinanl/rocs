/**
 * A class of interval paver.
 *
 * Created by Yinan Li on Aug. 16, 2016.
 *
 * Hybrid Systems Group, University of Waterloo.
 */

#ifndef _interval_paver_h_
#define _interval_paver_h_

#include <queue>
#include <stack>
#include <boost/dynamic_bitset.hpp>
#include <fstream>

#include "definitions.h"
#include "interval_vector.h"


namespace rocs {
  
    /**
     * \brief An enum defines set approximation type.
     */
    enum APPROXTYPE {

		     INNER, /**< Treat the undetermined intervals outsiders. */
		     OUTER, /**< Treat the undetermined intervals insiders. */
		     EXACT /**< Keey the undetermined intervals undetermined. */
    };

    /**
     * \brief A struct for subpaver node.
     */
    struct SPnode
    {
	ivec _box;			/**< The interval vector */
	boost::dynamic_bitset<> _cntl; /**< The control indicator */
	short _tag;			/**< Marking in or out */
	bool _b0, _b1;		/**< outside(b0=1), inside(b1=1) */
	UintSmall _label;  /**< The label(default=0) of the node (assigned by a labeling function) */
	int _split;			/**< The splitting dimension */
	SPnode *_left, *_right;	/**< Pointers to left & right childern */
  

	/**
	 * \brief The default constructor.
	 */
	SPnode() {}

	/**
	 * \brief A constructor.
	 *
	 * @param[in] b An interval vector
	 * @param[in] cn The number of control inputs
	 */
	SPnode(const ivec &b, const int cn) :
	    _box(b), _tag(0), _b0(true), _b1(false), _label(0),
	    _split(-1), _left(NULL), _right(NULL) {
	    
	    _cntl.resize(cn, false);
	}

	/**
	 * \brief A constructor.
	 *
	 * @param[in] b An interval vector
	 * @param[in] cn The control table
	 * @param[in] t The tagg
	 * @param[in] b0 1 if intersect with free space, and 0 otherwise
	 * @param[in] b1 1 if intersect with goal, and 0 otherwise
	 * @param[in] pp The label for the %SPnode
	 */
	SPnode(const ivec &b, const boost::dynamic_bitset<> &c,
	       const short t, const bool b0, const bool b1, const UintSmall pp) : 
	    _box(b), _cntl(c), _tag(t), _b0(b0), _b1(b1), _label(pp),
	    _split(-1), _left(NULL), _right(NULL) {} 

	/**
	 * \brief A copy constructor.
	 */
	SPnode(const SPnode &other) :
	    _box(other._box), _cntl(other._cntl),
	    _tag(other._tag), _b0(other._b0), _b1(other._b1), _label(other._label),
	    _split(other._split), _left(other._left), _right(other._right) {}

	/**
	 * \brief Test if the node is empty.
	 */
	bool isempty() const { return _box.isempty(); }
  
    };

    /**
     * \brief Print the information of a %SPnode.
     */
    std::ostream& operator<<(std::ostream&, const SPnode&);


    /**
     * \brief A subpaving tree class.
     * 
     * This is the data structure for organizing non-uniform grids. Nodes are splitted only when it is necessary.
     */
    class SPtree
    {
    public:

	SPnode *_root;		/**< The root of the SPtree */
  
	/**
	 * \brief The default constructor.
	 */
	SPtree(): _root(NULL) {}

	/**
	 * \brief A constructor.
	 *
	 * @param[in] root The pointer to the node to be copied
	 */
	SPtree(const SPnode *root) { _root = new SPnode(*root); }

	/**
	 * \brief A copy constructor (recursive).
	 *
	 * @param r The %SPtree to be copied
	 */
	SPtree(const SPtree &r);
	SPnode* copyHelper(const SPnode*);

	/**
	 * \brief A destructor.
	 */
	~SPtree() { release(_root);}

	/**
	 * \brief A copy assignment.
	 *
	 * @param other The SPtree to be copied
	 */
	SPtree& operator=(const SPtree &other);

	/**
	 * \brief Check if a node is a leaf.
	 *
	 * @param ptrn The pointer to the node
	 * @return 1(is leaf), 0(otherwise).
	 */
	bool isleaf(SPnode *ptrn) const;
	
	/**
	 * \brief Check if a node is empty.
	 *
	 * @return 1(is empty), 0(otherwise).
	 */
	bool isempty() const { return _root == NULL ? true : false;}


	/**
	 * \brief Release node memory (recursively).
	 *
	 * @param ptrn The pointer to the subtree to be released
	 */
	void release(SPnode *ptrn);
  
	/**
	 * \brief Expand a leaf by bisection, keeping all the values of the leaf.
	 *
	 * @param leaf The leaf to be expanded by bisection
	 * @param axis The splitting axis
	 */
	void expand(SPnode *leaf, int axis);
	
	/**
	 * \brief Expand a leaf by connection.
	 *
	 * @param leaf The leaf to be expanded by bisection
	 * @param lcld The left child
	 * @param rcld The right child
	 */
	void expand(SPnode *leaf, SPnode *lcld, SPnode *rcld);

	/**
	 * \brief Refine a node by another SPtree target. (deprecated)
	 *
	 * Assumption: \f$x\in\f$ sp, \f$x\f$ is a leaf of this tree.
	 * @param node The pointer to the node
	 * @param ptrn The pointer to the node to be refined
	 * @param ptrsp The pointer to the SPtree
	 * @param tag 
	 */
	void refine_leaf(SPnode *ptrn, SPnode *ptrsp, short tag);

	/**
	 * \brief Retract the tree.
	 *
	 * Traverse from bottom to top and merge the leaves tagged same.
	 */
	void retract();
  
	/**
	 * \brief Tagging the tree (update _tag by _b0 and _b1)
	 *
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
	 * @param approx The approximation method (INNER, OUTER or EXACT)
	 */
	void tagging(APPROXTYPE approx);

	/**
	 * \brief Tag negation (1-->0, 0-->1, 2-->2, -1-->-1).
	 */
	void negate();

	/**
	 * \brief Reset the tags of every nodes of the tree back to 0, 
	 * except for the ones tagged -1.
	 */
	void reset_tags();

	/**
	 * \brief Count leaves in preorder (depth-first).
	 *
	 * @param ptrn The pointer to the node from which the counted leaves branch
	 * @return the number of leaves under the node.
	 */
	size_t leafcount(SPnode *ptrn) const;
  
	/**
	 * \brief Count all nodes in preorder (depth-first).
	 *
	 * @param ptrn The pointer to the node from which the counted nodes branch
	 * @return the number of nodes under the given node.
	 */
	size_t nodecount(SPnode *ptrn) const;

	/**
	 * \brief Calculate the height of a node in a tree.
	 *
	 * - leaf height is 1.
	 * - number of nodes of a balanced BTree: 2^H - 1.
	 * @param ptrn The pointer to the node
	 * @return the node height.
	 */
	size_t height(SPnode *ptrn) const;

	/**
	 * \brief All leaves under a node.
	 *
	 * @param ptrn The pointer to the node from which the leaves branch
	 * @return a list of leaves.
	 */
	std::vector<SPnode*> leaves(SPnode *ptrn) const;

	/**
	 * \brief Print SPtree level-by-level (Breadth-first).
	 *
	 * @param ptrn The pointer to the node
	 */
	void print_subtree(SPnode *ptrn) const;
	
	/**
	 * \brief Print leaves in preorder (Depth-first).
	 *
	 * @param ptrn The pointer to the node
	 */
	void print_leaves(SPnode *ptrn) const;
	
	/**
	 * \brief Print leaves with a specified tag in preorder (Depth-first).
	 *
	 * @param ptrn The pointer to the node
	 * @param tag The given tag
	 */
	void print_leaves(SPnode *ptrn, short tag) const;
  
    };


} // namespace rocs

#endif
