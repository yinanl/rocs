/**
 * interval_paver.cpp
 *
 * A class of interval paver.
 *
 * Created by Yinan Li on Aug. 16, 2016.
 *
 * Hybrid Systems Group, University of Waterloo.
 */




#include "interval_paver.h"


namespace rocs {

    /* SPnode */
    std::ostream& operator<<(std::ostream &out, const SPnode &node) {

	if (node.isempty())
	    return out << "()";
    
	out << "(" << node._box << ", "  << node._label << ", " << node._tag << ", " << node._split;

	out << ", {";

	// bool first = true;
	// for (size_t i = 0; i < node._cntl.size(); ++i) {

	//     if (node._cntl[i]) {

	// 	if (first) {

	// 	    out << i;
	// 	    first = false;
	// 	} else {
		    
	// 	    out << ", " << i;
	// 	}
	//     }
	// }
	out << "})";

	return out;
    }


    /* SPtree */
    SPtree::SPtree(const SPtree &other) {

	_root = copyHelper(other._root);
    }

    SPnode* SPtree::copyHelper(const SPnode *other) {

	if ((*other).isempty()) {
	
	    return NULL;
	} else {

	    SPnode *node = new SPnode(*other);

	    if (other->_left) {

		node->_left = copyHelper(other->_left);
	    } else {

		node->_left = NULL;
	    }

	    if (other->_right) {

		node->_right = copyHelper(other->_right);
	    } else {

		node->_right = NULL;
	    }

	    return node;
	}
    }

    SPtree& SPtree::operator=(const SPtree &other) {

	if (this != &other)
	    _root = copyHelper(other._root);
	
	return (*this);
    }


    bool SPtree::isleaf(SPnode *node) const {

	if (node != NULL) {

	    return node->_left == NULL && node->_right == NULL;
	}

	return false;
    }


    void SPtree::release(SPnode *node) {

	if (node != NULL) {

	    release(node->_left);
	    release(node->_right);

	    delete node;
	}
    }

    void SPtree::expand(SPnode *leaf, int axis) {

	if (!isleaf(leaf))
	    return;

	/* left right childern inherit _tag and _cntl */
	leaf->_left = new SPnode(lowerhalf(leaf->_box, axis), leaf->_cntl,
				 leaf->_tag, leaf->_b0, leaf->_b1, leaf->_label);
	leaf->_right = new SPnode(upperhalf(leaf->_box, axis), leaf->_cntl,
				  leaf->_tag, leaf->_b0, leaf->_b1, leaf->_label);
	leaf->_split = axis;
    }

    void SPtree::expand(SPnode *leaf, SPnode *lcld, SPnode *rcld) {

	if (!isleaf(leaf))
	    return;

	leaf->_left = lcld;
	leaf->_right = rcld;
    }


    void SPtree::refine_leaf(SPnode *x, SPnode *sp, short tag) {

	assert(sp->_box.isin(x->_box));

	if (isleaf(sp)) {

	    if (sp->_tag == 1) {

		x->_tag = tag;
		x->_b0 = false;
		x->_b1 = true;
	    }
	    else {

		x->_tag = 0;
		x->_b0 = true;
		x->_b1 = false;
	    }
	}
	else {
	
	    if (sp->_left->_box.isin(x->_box)) {
	
		refine_leaf(x, sp->_left, tag);
	    }
	    else if (sp->_right->_box.isin(x->_box)) {

		refine_leaf(x, sp->_right, tag);
	    }
	    else {
		/* split x box into two w.r.t. sp tree structure */
		x->_split = sp->_split;
	    
		SPnode *ptrlx = new SPnode(*x);
		SPnode *ptrrx = new SPnode(*x);
		ptrlx->_box[x->_split].setsup(sp->_left->_box[x->_split].getsup());
		ptrrx->_box[x->_split].setinf(sp->_right->_box[x->_split].getinf());
	    

		expand(x, ptrlx, ptrrx);  // connect two childern to x
	    

		refine_leaf(x->_left, sp->_left, tag);
		refine_leaf(x->_right, sp->_right, tag);
	    
	    }
	}
    
    }


    /* Update tree: tagging, retracting wrt tags */
    void SPtree::retract() {

	/* store internal nodes in reverse order */
	std::stack<SPnode*> stk;
	std::queue<SPnode*> que;
	que.push(_root);

	SPnode *current;
	while (!que.empty()) {

	    current = que.front();
	    que.pop();

	    if (current->_right || current->_left)
		stk.push(current);  // only store internal

	    if (current->_right)
		que.push(current->_right);

	    if (current->_left)
		que.push(current->_left);
	}

	while (!stk.empty()) {

	    current = stk.top();
	    stk.pop();
	
	    /* retract nodes only if the childern are leaves */
	    if (current->_left->_tag == current->_right->_tag &&
		isleaf(current->_left) && isleaf(current->_right) &&
		current->_left->_label == current->_right->_label) {

		current->_tag = current->_left->_tag;
		release(current->_left);
		release(current->_right);
		current->_left = NULL;
		current->_right = NULL;
	    }
	}

	return;
    }


    void SPtree::tagging(APPROXTYPE approx) {

	std::stack<SPnode*> stk;
    
	std::queue<SPnode*> que;
	que.push(_root);

	SPnode *current;
	while (!que.empty()) {

	    current = que.front();
	    que.pop();

	    stk.push(current);  // store nodes in reverse order

	    if (current->_right)
		que.push(current->_right);

	    if (current->_left)
		que.push(current->_left);
	}

    
	while (!stk.empty()) {  //tagging by its childern

	    current = stk.top();
	    stk.pop();

	    if (isleaf(current)) {

		if (current->_tag != -1) {  // only update those not -1
		    /* update _tag by _b0 and _b1 */
		    switch (approx) {

		    case INNER:
			current->_tag = (!current->_b0 && current->_b1) ? 1 : 0;
			break;

		    case OUTER:
			current->_tag = current->_b1 ? 1 : 0;
			break;
		
		    case EXACT:
			if (current->_b0 && !current->_b1)
			    current->_tag = 0;
			else if (!current->_b0 && current->_b1)
			    current->_tag = 1;
			else
			    current->_tag = 2;
			break;
		
		    }  // end switch
		}  // end if
	    
	    }
	    else {

		if (current->_left->_tag == current->_right->_tag)
		    current->_tag = current->_left->_tag;
		else if (current->_left->_tag < 1 && current->_right->_tag < 1)
		    current->_tag = -2;
		else
		    current->_tag = 2;
	    }
	}
    }


    void SPtree::negate() {

	if (!isempty()) {

	    std::vector<SPnode*> lev = leaves(_root);
	    for (size_t i = 0; i < lev.size(); ++i) {
		if (lev[i]->_tag > 0 && lev[i]->_tag < 2)
		    lev[i]->_tag ^= 1;  // flip between 0 and 1
	    }

	    tagging(EXACT);
	}
    
    }


    void SPtree::reset_tags() {
	std::stack<SPnode*> stk;
	stk.push(_root);
	SPnode *current;
	while(!stk.empty()) {
	    current = stk.top();
	    stk.pop();
	    if(current->_tag != -1) {
		current->_b0 = true;
		current->_b1 = false;
		current->_tag = 0;
	    }
	    if(current->_right)
		stk.push(current->_right);
	    if(current->_left)
		stk.push(current->_left);
	}

	tagging(EXACT);
    }

    /* Count numbers */
    size_t SPtree::leafcount(SPnode *node) const {
    
	if (node == NULL) 
	    return 0;
    
	if (isleaf(node)) 
	    return 1;
	else
	    return leafcount(node->_left) + leafcount(node->_right);
    }

    size_t SPtree::nodecount(SPnode *node) const {
    
	if (node == NULL) 
	    return 0;
    
	if (isleaf(node)) 
	    return 1;
	else
	    return 1 + nodecount(node->_left) + nodecount(node->_right);
    }

    size_t SPtree::height(SPnode *node) const {
    
	if (node == NULL) 
	    return 0;

	size_t lh = height(node->_left);
	size_t rh = height(node->_right);
    
	if (lh >= rh)
	    return lh + 1;
	else
	    return rh + 1;
    }

    std::vector<SPnode*> SPtree::leaves(SPnode *node) const {

	std::vector<SPnode*> lef;
	
	if (node->isempty())
	    return lef;
    
	std::stack<SPnode*> stk;
	stk.push(node);

	SPnode *current;
	while (!stk.empty()) {

	    current = stk.top();
	    stk.pop();

	    if (isleaf(current)) {

		lef.push_back(current);
	    }
	    else {

		if (current->_right)
		    stk.push(current->_right);

		if (current->_left)
		    stk.push(current->_left);
	    }
	}

	return lef;
    }

    /* display functions */
    void SPtree::print_subtree(SPnode *root) const {

	if (root == NULL)
	    return;

	std::queue<SPnode*> q;
	std::queue<int> lq;
    
	SPnode *current = root;
	int level, levelold;

	q.push(current);
	lq.push(0);

	levelold = -1;
	while (!q.empty()) {

	    current = q.front();
	    level = lq.front();
	
	    /* do the printing */
	    if (levelold != -1 && level > levelold)
		std::cout << std::endl;
	
	    std::cout << *current << ", ";
	

	    if (current->_left != NULL) {
	    
		q.push(current->_left);
		lq.push(level + 1);
	    }

	    if (current->_right != NULL) {
	    
		q.push(current->_right);
		lq.push(level + 1);
	    }

	    q.pop();
	    lq.pop();
	
	    levelold = level;
	}
    
    }

    void SPtree::print_leaves(SPnode *node) const {
    
	if (node == NULL) 
	    return;

	if (isleaf(node)) {

	    std::cout << *node << ", ";
	} else {
	
	    if (node->_left != NULL) //recursions
		print_leaves(node->_left);
    
	    if (node->_right != NULL)
		print_leaves(node->_right);

	}
    
    }

    void SPtree::print_leaves(SPnode *node, short tag) const {
    
	if (node == NULL) 
	    return;

	if (isleaf(node)) {

	    if (node->_tag == tag)
		std::cout << *node << ", ";
	
	} else {
	
	    if (node->_left != NULL) //recursions
		print_leaves(node->_left, tag);
    
	    if (node->_right != NULL)
		print_leaves(node->_right, tag);

	}
    
    }
    
} // namespace rocs
