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

    int CSolver::bisect_axis(ivec &box) {

	int axis = 0;
	double ri, r = 0;
    
	switch (_bstype) {
	case ABSMAX:
	    axis = box.maxdim();
	    break;

	case RELMAXW:
	
	    for (int i = 0; i < _xdim; ++i) {
		ri = box[i].width()/_ctlr._root->_box.width()[i];
		if (ri > r) {
		    axis = i;
		    r = ri;
		}  // end if
	    }  // end for
	    break;
	case RELMAXG:
	    for (int i = 0; i < _xdim; ++i) {
		ri = box[i].width()/_goalsize[i];
		if (ri > r) {
		    axis = i;
		    r = ri;
		}  // end if
	    }  // end for
	    break;
	default:
	    break;
	}

	return axis;
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
	short itag = (ap == GOAL) ? 1 : -1;
	if (_ctlr.isleaf(_ctlr._root)) { // first time
	
	    paver_init(_ctlr._root, area, itag);
	}
	else {  // has been initialized, so refine _ctlr
	
	    init_refine(_ctlr._root, area, itag);
	}

	_ctlr.tagging(EXACT);
    
    }


    void CSolver::paver_init(SPnode *current, ivec &cbox, short itag) {

	if (cbox.isempty())
	    return;

	bool b1 = (itag == 1) ? true : false;

	ivec rbox, lbox;
    
	for (int i = 0; i < cbox.getdim(); ++ i) {

	    rbox = current->_box;
	    lbox = current->_box;

	    if (cbox[i].getsup() < rbox[i].getsup()) {
	    
		rbox[i].setinf(cbox[i].getsup());
		lbox[i].setsup(cbox[i].getsup());

		/* keep parent's tag, and no childern */
		current->_right = new SPnode(rbox, current->_cntl, current->_tag,
					     current->_b0, current->_b1); 

		current->_left = new SPnode(lbox, current->_cntl, current->_tag,
					    current->_b0, current->_b1);
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
						 current->_b0, current->_b1);
		}
		else {

		    current->_right = new SPnode(rbox, current->_cntl, itag, false, b1);
		}

		current->_left = new SPnode(lbox, current->_cntl, current->_tag,
					    current->_b0, current->_b1);
		current->_split = i;
	    
		current = current->_right; //expand the right child
	    }
	    else {

		if (i >= cbox.getdim() - 1) {

		    current->_tag = itag;
		    current->_b0 = false;
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


    void CSolver::init(SPEC ap, fcst f, const double eps/* 0.01 */) {

	assert(!_ctlr.isempty());

	bool inner = (ap == GOAL) ? true : false;
	short itag = (ap == GOAL) ? 1 : -1;

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

	// compute_winsize();
    }

    /* Initialize function constraint:
     *
     * f -- the constraint (defined by a function) 
     * f(x) <= 0
     */
    void CSolver::paver_init(SPtree &sp, fcst f, bool inner, short itag,
			     const double eps/* 0.01 */) {

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


    void CSolver::sivia(SPtree &sp, SPnode *ptrnode, ivec &cst, fcst fcn,
			bool inner, short itag, const double eps/* 0.01 */) {

	if (sp.isempty() || ptrnode == NULL)
	    return;
    
	bool b1 = (itag == 1) ? true : false;

	/* evaluate mapping */
	ivec box = ptrnode->_box;
	ivec fbox = fcn(box);
    
	if (cst.isin(fbox)) { // inside change _tag
	    if (ptrnode->_tag != -1) {
		ptrnode->_tag = itag;
		ptrnode->_b0 = false;
		ptrnode->_b1 = b1;
	    }
	    return;
	}

	if (cst.isout(fbox)) // outside _tag unchanged
	    return;
	
	if (box.maxwidth() < eps) { // undetermined
	    if (!inner && ptrnode->_tag != -1) { // outer approximation
		ptrnode->_tag = itag;
		ptrnode->_b0 = false;
		ptrnode->_b1 = b1;
	    } // _tag unchanged if inner
	    return;
	}

	/* precision not satisfied: expansion */
	sp.expand(ptrnode, bisect_axis(box));  // bisect relatively maximum axis

	sivia(sp, ptrnode->_left, cst, fcn, inner, itag, eps);
	sivia(sp, ptrnode->_right, cst, fcn, inner, itag, eps);
    
    }

    /* refine node in _ctlr by new constraint function */
    void CSolver::init_refine(SPtree &sp, fcst f, bool inner, short itag,
			      const double eps/* 0.01 */) {

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
	compute_goalsize();
	compute_winsize();
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

    void CSolver::compute_goalsize() {
        std::vector<double> lb(_xdim, PINF);
	std::vector<double> ub(_xdim, NINF);

        std::stack<SPnode*> stk;
        stk.push(_ctlr._root);

        SPnode *current;
        std::vector<double> ubs, lbs;
        while (!stk.empty()) {

	    current = stk.top();
	    stk.pop();

	    if (current->_tag == 1) {
	    	lbs = current->_box.getinf();
	    	ubs = current->_box.getsup();
	    	for (int i = 0; i < _xdim; ++i) {

	    	    lb[i] = lb[i] > lbs[i] ? lbs[i] : lb[i];
	    	    ub[i] = ub[i] < ubs[i] ? ubs[i] : ub[i];
	    	}
	    }
	    else {
	    	if (!_ctlr.isleaf(current)) {
		
	    	    if (current->_right)
	    		stk.push(current->_right);

	    	    if (current->_left)
	    		stk.push(current->_left);
	    	}
	    
	    }
        }

	// std::cout << _goalsize.size() << "\n";
	// std::cout << _xdim << '\n';
	// std::cout << "The goal size is (";
        for (int i = 0; i < _goalsize.size(); ++i) {
	    _goalsize[i] = ub[i] - lb[i];
	    // std::cout << _goalsize[i] << ", ";
        }

        // std::cout << ")\n";
    }


    short CSolver::paver_test(SPtree &sp, ivec &box) {

	short t = 1;

	std::queue<SPnode*> q;
	q.push(sp._root);

	bool t0 = false;
	bool t1 = false;
	SPnode *current;
	while (!q.empty()) {

	    current = q.front();
	    q.pop();

	    switch (current->_tag) {

	    case 1 :
		if (current->_box.isin(box))
		    return 1;
		if (!current->_box.isout(box))
		    t1 = true;
		break;

	    case 2 :
		// only follow the branch that has intersection
		if (!current->_box.isout(box)) {

		    if (sp.isleaf(current))  // if leaves can only be 1 or 0,
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
		if (current->_box.isin(box))  // box inside current
		    return 0;
		if (!current->_box.isout(box))  // box intersect current
		    t0 = true;
		break;
	    }

	    if (t0 && t1)
		break;  // no need to continue the loop
	}

	if (!t1)
	    t = 0;  // t0 & t1 all false: box is out of domain => out
	else if (!t0)
	    t = 1;
	else  // t1=1, t0=1
	    t = 2;

    
	return t;
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
	    }
	    else {

		if (!_ctlr.isleaf(current)) {
		
		    if (current->_right)
			stk.push(current->_right);

		    if (current->_left)
			stk.push(current->_left);
		}
	    
	    }  // end if
	
	}  // end while

	_winsize = std::pow(vol, 1.0/_xdim);
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
    
    }


    void CSolver::print_controller_info() const {

	if (_ctlr.isempty()) {

	    std::cout << "No controller generated.\n";
	}
	else {

	    std::cout << "Number of partitions: "
		      << _ctlr.leafcount(_ctlr._root) << '\n';
	
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
	logfile.open("log.txt", std::ios::ate | std::ios::app); // append data

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
