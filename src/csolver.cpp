/** 
 *  csolver.cpp
 *  
 *  An interval based control problem solver class.
 *
 *  Created by Yinan Li on Jan. 03, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <time.h>

#include "csolver.h"


namespace rocs {

    void CSolver::set_M(const std::vector<UintSmall> &outedge) {
	UintSmall NofProps = outedge.size();
	if (_M.size() != NofProps)
	    _M.resize(NofProps);
	for (UintSmall i = 0; i < NofProps; ++i)
	    _M[i] = outedge[i];
    }

    int CSolver::bisect_axis(ivec &box, const double eps[]) {

	int axis = 0;
	double ri, r = 0;
    
	switch (_bstype) {
	case ABSMAX:
	    axis = box.maxdim();
	    break;
	case RELMAX:
	    for (int i = 0; i < _xdim; ++i) {
		ri = box[i].width()/eps[i];
		if (ri > r) {
		    axis = i;
		    r = ri;
		}
	    }
	    break;
	default:
	    break;
	}

	return axis;
    }

    void CSolver::labeling(const double lb[], const double ub[], UintSmall prop) {
	ivec area(_xdim);
	for (int i = 0; i < _xdim; ++i)
	    area[i] = interval(lb[i], ub[i]);
	
	labeling(area, prop);
    }
    void CSolver::labeling(ivec &area, UintSmall prop) {
	assert(!_ctlr.isempty());
	if (_ctlr.isleaf(_ctlr._root)) {
	    init_label(_ctlr._root, area, prop);
	} else {
	    refine_label(_ctlr._root, area, prop);
	}
    }
    void CSolver::init_label(SPnode *current, ivec &cbox, UintSmall prop) {
	if (cbox.isempty() || current->_tag == -1)
	    return;

	ivec rbox, lbox;
	for (int i = 0; i < cbox.getdim(); ++ i) {
	    rbox = current->_box;
	    lbox = current->_box;

	    if (cbox[i].getsup() < rbox[i].getsup()) {
		rbox[i].setinf(cbox[i].getsup());
		lbox[i].setsup(cbox[i].getsup());
		/* keep parent's tag, and no childern */
		current->_right = new SPnode(rbox, current->_cntl, current->_tag,
					     current->_b0, current->_b1, current->_label); 

		current->_left = new SPnode(lbox, current->_cntl, current->_tag,
					    current->_b0, current->_b1, current->_label);
		current->_split = i;
		current = current->_left; // expand the left child
		rbox = current->_box;
		lbox = current->_box;
	    } /* otherwise, no right node */

	    if (cbox[i].getinf() > lbox[i].getinf()) {
		rbox[i].setinf(cbox[i].getinf());
		lbox[i].setsup(cbox[i].getinf());
		if (i < cbox.getdim() - 1) {
		    current->_right = new SPnode(rbox, current->_cntl, current->_tag,
						 current->_b0, current->_b1, current->_label);
		} else {
		    /* assign the label if i is the last dimension */
		    current->_right = new SPnode(rbox, current->_cntl, current->_tag,
						 current->_b0, current->_b1, prop);
		    // std::cout << "logging init_label:\n" << current->_right->_label << std::endl;
		}
		current->_left = new SPnode(lbox, current->_cntl, current->_tag,
					    current->_b0, current->_b1, current->_label);
		current->_split = i;
		current = current->_right; /* expand the right child */
	    } else {
		/* no need to split */
		if (i >= cbox.getdim() - 1) {
		    /* if i is the last dimension, the whole node is a leaf and assign the label. */
		    current->_label= prop;
		    // std::cout << "logging init_label:\n" << current->_label << std::endl;
		} //end if
	    }//end if
	
	}  //end for
    }
    void CSolver::refine_label(SPnode *node, ivec &box, UintSmall prop) {
	assert(node->_box.isin(box));  // assume that box inside node

	if (_ctlr.isleaf(node)) {
	    // std::cout << "refine_label 0:" << node->_box << std::endl;
	    init_label(node, box, prop);
	} else {
	    if (node->_left->_box.isin(box)) {
		// std::cout << "refine_label 1:" << node->_left->_box << std::endl;
		refine_label(node->_left, box, prop);
	    } else if (node->_right->_box.isin(box)) {
		// std::cout << "refine_label 2:" << node->_right->_box << std::endl;
		refine_label(node->_right, box, prop);
	    } else {  
		/* split box along the _split */
		ivec lbox(box), rbox(box);
		lbox[node->_split].setsup(node->_left->_box[node->_split].getsup());
		rbox[node->_split].setinf(node->_right->_box[node->_split].getinf());
		// std::cout << "refine_label 3:" << lbox << std::endl;
		// std::cout << "refine_label 3:" << rbox << std::endl;
		refine_label(node->_left, lbox, prop);
		refine_label(node->_right, rbox, prop);
	    }// end if
	}// end if
    }

    void CSolver::init_winset() {
	/* The winning set is initialized to the whole state space:
	 * set _tag of every leaf to be 1 if it is not an avoid node 
	 */
	std::stack<SPnode*> stk;
	stk.push(_ctlr._root);
	SPnode *current;
	while (!stk.empty()) {
	    current = stk.top();
	    stk.pop();
	    if (_ctlr.isleaf(current)) {
		if (current->_tag >= 0) {
		    current->_b0 = false;
		    current->_b1 = true;
		    current->_tag = 1;
		}
	    } else {
		if (current->_right)
		    stk.push(current->_right);
		if (current->_left)
		    stk.push(current->_left);
	    }
	}
	_ctlr.tagging(EXACT);
    }
    void CSolver::init_winset(ivec &area) {
	assert(!_ctlr.isempty());
	if (area == _ctlr._root->_box) {
	    init_winset();
	} else {
	    init(GOAL, area);
	}
    }

    bool CSolver::targetset_in_sdoms(std::vector<SPtree*> &sdoms) {
	assert(!_M.empty());
	for (size_t i = 0; i < _M.size(); ++i) {
	    if (_M[i] < sdoms.size()) {
		/* Assume that root tag of the controller _ctlr has been updated.
		 * _tag = 1 or 2 if there are some leaves with tag = 1.
		 */
		if (sdoms[_M[i]]->_root->_tag > 0)
		    return true;
	    }
	}
	return false;
    }

    void CSolver::init(SPEC ap, const double lb[], const double ub[]) {
	assert(!_ctlr.isempty());
	ivec area(_xdim);
	for (int i = 0; i < _xdim; ++i) {      
	    area[i] = interval(lb[i], ub[i]);
	}
	init(ap, area);
    }
    void CSolver::init(SPEC ap, ivec &area) {
	assert(!_ctlr.isempty());
	short itag;
	switch(ap) {
	case GOAL:
	    itag = 1; break;
	case AVOID:
	    itag = -1; break;
	case FREE:
	    itag = 0; break;
	default:
	    std::cout << "CSolver::init: Specification should be one of (GOAL, AVOID, FREE).\n";
	    exit(EXIT_FAILURE);
	}
	// short itag = (ap == GOAL) ? 1 : -1;
	if (_ctlr.isleaf(_ctlr._root)) { // first time
	    paver_init(_ctlr._root, area, itag);
	} else {  // has been initialized, so refine _ctlr
	    init_refine(_ctlr._root, area, itag);
	}
	_ctlr.tagging(EXACT);

	compute_winsize();  // get the initial volume of the winning set
    }


    void CSolver::paver_init(SPnode *current, ivec &cbox, short itag) {

	if (cbox.isempty())
	    return;

	bool b1 = (itag == 1) ? true : false;
	bool b0 = (itag == 0) ? true : false;

	ivec rbox, lbox;
    
	for (int i = 0; i < cbox.getdim(); ++ i) {

	    rbox = current->_box;
	    lbox = current->_box;

	    if (cbox[i].getsup() < rbox[i].getsup()) {
	    
		rbox[i].setinf(cbox[i].getsup());
		lbox[i].setsup(cbox[i].getsup());

		/* keep parent's tag, and no childern */
		current->_right = new SPnode(rbox, current->_cntl, current->_tag,
					     current->_b0, current->_b1, current->_label); 

		current->_left = new SPnode(lbox, current->_cntl, current->_tag,
					    current->_b0, current->_b1, current->_label);
		current->_split = i;

		current = current->_left; // expand the left child
		rbox = current->_box;
		lbox = current->_box;
	    }

	    if (cbox[i].getinf() > lbox[i].getinf()) {
	    
		rbox[i].setinf(cbox[i].getinf());
		lbox[i].setsup(cbox[i].getinf());

		if (i < cbox.getdim() - 1) {

		    current->_right = new SPnode(rbox, current->_cntl, current->_tag,
						 current->_b0, current->_b1, current->_label);
		}
		else {

		    current->_right = new SPnode(rbox, current->_cntl, itag, b0, b1, current->_label);
		}

		current->_left = new SPnode(lbox, current->_cntl, current->_tag,
					    current->_b0, current->_b1, current->_label);
		current->_split = i;
	    
		current = current->_right; //expand the right child
	    }
	    else {

		if (i >= cbox.getdim() - 1) {

		    current->_tag = itag;
		    current->_b0 = b0;
		    current->_b1 = b1;
		}
	    
	    }//end if
	
	}  //end for
    
    }


    void CSolver::init_refine(SPnode *node, ivec &box, short itag) {

	assert(node->_box.isin(box));  // assume that box inside node

	if (_ctlr.isleaf(node)) {
	
	    paver_init(node, box, itag);
	}
	else {

	    if (node->_left->_box.isin(box)) {
	
		init_refine(node->_left, box, itag);
	    }
	    else if (node->_right->_box.isin(box)) {

		init_refine(node->_right, box, itag);
	    }
	    else {  // split box along the right axis

		/* find the split dimension */
		int axis = node->_split;
		ivec lbox = box;
		ivec rbox = box;

		lbox[axis].setsup(node->_left->_box[axis].getsup());
		rbox[axis].setinf(node->_right->_box[axis].getinf());

		init_refine(node->_left, lbox, itag);
		init_refine(node->_right, rbox, itag);
	    
	    }// end if
	
	}// end if
    }


    void CSolver::init(SPEC ap, fcst f, const double eps[]) {
	assert(!_ctlr.isempty());

	// bool inner = (ap == GOAL) ? true : false;
	// short itag = (ap == GOAL) ? 1 : -1;
	bool inner;
	short itag;
	switch(ap) {
	case GOAL:
	    itag = 1; inner = true; break;
	case AVOID:
	    itag = -1; inner = false; break;
	case FREE:
	    itag = 0; inner = false; break;
	default:
	    std::cout << "CSolver::init: Specification should be one of (GOAL, AVOID, FREE).\n";
	    exit(EXIT_FAILURE);
	}

	if (_ctlr.isleaf(_ctlr._root)) {
	
	    paver_init(_ctlr, f, inner, itag, eps);
	}
	else {
	
	    init_refine(_ctlr, f, inner, itag, eps);
	}

	if (inner) // inner approximation if goal
	    _ctlr.tagging(INNER);
	else // outer if obstacle
	    _ctlr.tagging(OUTER);

	compute_winsize();
    }

    /* Initialize function constraint:
     *
     * f -- the constraint (defined by a function) 
     * f(x) <= 0
     */
    void CSolver::paver_init(SPtree &sp, fcst f, bool inner, short itag,
			     const double eps[]) {

	if (sp.isempty())
	    return;

	/* Y = [-oo, 0] */
	ivec y = f(sp._root->_box); // calculate f in order to get a correct y dimension
	for (int i = 0; i < y.getdim(); ++i) {
	    y[i] = interval(NINF, 0);
	}

	/* expand SPtree by f(x)\in y */
	sivia(sp, sp._root, y, f, inner, itag, eps);
    }

    /* refine node in _ctlr by new constraint function */
    void CSolver::init_refine(SPtree &sp, fcst f, bool inner, short itag,
			      const double eps[]) {

	if (sp.isempty())
	    return;

	std::vector<SPnode*> lev = sp.leaves(sp._root);
    
	ivec y = f(sp._root->_box);
	for (int i = 0; i < y.getdim(); ++i) {
	    y[i] = interval(NINF, 0);
	}

	/* refine each leaf */
	size_t nl = sp.leafcount(sp._root);
	for (size_t i = 0; i < nl; ++i) {
	    sivia(sp, lev[i], y, f, inner, itag, eps);
	}
    
    }

    void CSolver::sivia(SPtree &sp, SPnode *ptrnode, ivec &cst, fcst fcn,
			bool inner, short itag, const double eps[]) {
	/* the node (*ptrnode) is not processed if it is empty/null or an avoid region */
	if (sp.isempty() || ptrnode == NULL || ptrnode->_tag == -1)
	    return;
    
	bool b1 = (itag == 1) ? true : false;
	bool b0 = (itag == 0) ? true : false;

	/* evaluate mapping */
	ivec box = ptrnode->_box;
	ivec fbox = fcn(box);
    
	if (cst.isin(fbox)) { // inside change _tag
	    ptrnode->_tag = itag;
	    ptrnode->_b0 = b0;
	    ptrnode->_b1 = b1;
	    return;
	}

	if (cst.isout(fbox)) // outside _tag unchanged
	    return;

	/* compute the split axis */
	int axis = bisect_axis(box, eps);
	if (box[axis].width() < eps[axis]) {
	    if (!inner) { // outer approximation
		ptrnode->_tag = itag;
		ptrnode->_b0 = b0;
		ptrnode->_b1 = b1;
	    } // _tag unchanged if inner
	    return;
	}
	sp.expand(ptrnode, axis);
	sivia(sp, ptrnode->_left, cst, fcn, inner, itag, eps);
	sivia(sp, ptrnode->_right, cst, fcn, inner, itag, eps);
    
    }

    
    void CSolver::init_goal_area() {
	std::stack<SPnode*> stk;
	stk.push(_ctlr._root);
    
	SPnode *current;
	while (!stk.empty()) {

	    current = stk.top();
	    stk.pop();

	    if (_ctlr.isleaf(current)) {

		if (current->_tag == 1) {
		    _goal.push_back(current->_box);
		}
	    }
	    else {

		if (current->_right)
		    stk.push(current->_right);

		if (current->_left)
		    stk.push(current->_left);
	    }
	} // end while
    }

    void CSolver::init_avoid_area() {
	std::stack<SPnode*> stk;
	stk.push(_ctlr._root);
    
	SPnode *current;
	while (!stk.empty()) {

	    current = stk.top();
	    stk.pop();

	    if (_ctlr.isleaf(current)) {

		if (current->_tag == -1) {
		    _obs.push_back(current->_box);
		}
	    }
	    else {

		if (current->_right)
		    stk.push(current->_right);

		if (current->_left)
		    stk.push(current->_left);
	    }
	} // end while
    }

    void CSolver::compute_winsize() {
    	double v, vol = 0;
    	std::vector<double> w(_xdim);
    	std::stack<SPnode*> stk;
    	stk.push(_ctlr._root);
    	SPnode *current;
    	while (!stk.empty()) {
    	    current = stk.top();
    	    stk.pop();
    	    v = 1;
    	    if (current->_tag == 1) {
    		w = current->_box.width();
    		for (int i = 0; i < _xdim; ++i) {
    		    v *= w[i];
    		}
    		vol += v;
    	    } else {
    		// if (!_ctlr.isleaf(current)) {
    		//     if (current->_right)
    		// 	stk.push(current->_right);
    		//     if (current->_left)
    		// 	stk.push(current->_left);
    		// }
		if (current->_right)
		    stk.push(current->_right);
		if (current->_left)
		    stk.push(current->_left);
    	    }  // end if
    	}  // end while
    	// _winsize = std::pow(vol, 1.0/_xdim);
	_winsize = v;
    }

    short CSolver::paver_test(SPtree &sp, ivec &box) {
	
	std::queue<SPnode*> q;
	q.push(sp._root);

	bool t0(false), t1(false);

	/**
	 * Handle out-of-domain scenarios
	 */
	/* If box is fully out of domain */
	if (_ctlr._root->_box.isout(box))
	    return 0;
	/* If box is partially out of domain */
	if(!_ctlr._root->_box.isin(box))
	    t0 = true;
	
	SPnode *current;
	while (!q.empty()) {

	    current = q.front();
	    q.pop();

	    switch (current->_tag) {
	    case 1 :
		if (current->_box.isin(box))  // box must be inside domain
		    return 1;
		if (!current->_box.isout(box))
		    t1 = true;
		break;

	    case 2 :
		// only follow the branch that has intersection
		if (!current->_box.isout(box)) {
		    if (current->_left == NULL && current->_right == NULL)  // if leaves can only be 1 or 0,
			// then this line is not necessary
			return 2;
		    else {
			if (current->_right)
			    q.push(current->_right);
			if (current->_left)
			    q.push(current->_left);
		    }
		}
		break;

	    default :  // for 0, -1, -2
		if (current->_box.isin(box))  // box inside current => inside domain
		    return 0;
		if (!current->_box.isout(box))  // box intersect current
		    t0 = true;
		break;
	    }

	    if (t0 && t1)
		return 2;  // no need to continue the loop
	}

	if (!t1)
	    return 0;  // t0 & t1 all false: box is out of domain => out
	else if (!t0)
	    return 1;
	else  // t1=1, t0=1
	    return 2;
	
    }


    /* fixed-point algorithms */
    void CSolver::init_leafque(std::stack<SPnode*> &l0,
			       std::stack<SPnode*> &l1,
			       std::stack<SPnode*> &l2) {

	std::stack<SPnode*> stk;
	stk.push(_ctlr._root);
    
	SPnode *current;
	while (!stk.empty()) {

	    current = stk.top();
	    stk.pop();

	    if (_ctlr.isleaf(current)) {

		if (current->_tag == 0) {
		    l0.push(current);
		}

		if (current->_tag == 1) {
		    l1.push(current);
		}

		if (current->_tag == 2) {
		    l2.push(current);
		}
	    }
	    else {

		if (current->_right)
		    stk.push(current->_right);

		if (current->_left)
		    stk.push(current->_left);
	    }
	} // end while
    
    }


    void CSolver::print_controller() const {

	assert( !_ctlr.isempty() );
    
	_ctlr.print_leaves(_ctlr._root);
	std::cout << '\n';
    
    }


    void CSolver::print_controller_info() const {

	if (_ctlr.isempty()) {
	    std::cout << "No controller generated.\n";
	} else {
	    std::cout << "Number of partitions: "
		      << _ctlr.leafcount(_ctlr._root) << '\n';
	    std::cout << "Controller tree height: "
		      << _ctlr.height(_ctlr._root) << '\n';
	    std::cout << "Total number of nodes in Interval Paver: "
		      << _ctlr.nodecount(_ctlr._root) << '\n';
	    std::cout << "Number of iterations:"
		      << _fpiter[0] << ','
		      << _fpiter[1] << ','
		      << _fpiter[2] << ',' << '\n';
	
	    std::cout << "Time of solving: "
		      << _timer << '\n';
	}
    }


    void CSolver::log_iterations(const char* filename, int iter) {

	std::fstream logfile;
	logfile.open(filename, std::ios::ate | std::ios::app); // append data

	logfile << iter << ":\n";
    
	std::stack<SPnode*> stk;
	stk.push(_ctlr._root);
	SPnode *current;
	std::vector<double> lower, upper;
	while (!stk.empty()) {

	    current = stk.top();
	    stk.pop();

	    if (_ctlr.isleaf(current) && current->_tag == 1) { //write leaves

		lower = current->_box.getinf();
		upper = current->_box.getsup();

		for (int i = 0; i < _xdim; ++ i)
		    logfile << lower[i] << ' ' << upper[i] << ' ';

		logfile << '\n';
	    }

	    if (current->_left)
		stk.push(current->_left);

	    if (current->_right)
		stk.push(current->_right);
	}
	logfile << '\n';
    
	logfile.close();
    }

    
} // namespace rocs
