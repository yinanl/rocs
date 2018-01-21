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
	
	for (int i = 0; i < _cpb->_xdim; ++i) {

	    ri = box[i].width()/_cpb->_workspace[i].width();
	    if (ri > r) {
		axis = i;
		r = ri;
	    }  // end if
	}  // end for
	break;
    case RELMAXG:
	for (int i = 0; i < _cpb->_xdim; ++i) {

	    ri = box[i].width()/_cpb->_goal[i].width();
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

    ivec area(_cpb->_xdim);
    for (int i = 0; i < _cpb->_xdim; ++i) {
      
      area[i] = interval(lb[i], ub[i]);
    }

    // short itag = (ap == GOAL) ? 1 : -1;
    short itag;
    switch (ap) {

    case GOAL:
    	itag = 1;
	_cpb->_goal = area;
	_winsize = area.maxwidth(); // winset=goal
    	break;

    case AVOID:
    	itag = -1;
	_cpb->_obs.push_back(area);
    	break;
	
    default:
	
    	break;
    }
    
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
	    // int axis = 0;
	    // for (int i = 0; i < _cpb->_xdim; ++i) {

	    // 	if (node->_left->_box[i] != node->_right->_box[i]) {
		    
	    // 	    axis = i;
	    // 	    break;
	    // 	}
	    // }// end for

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
    
    _ctlr.tagging(EXACT);

    compute_winsize();
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
    ivec y = f(sp._root->_box);
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

    // short utag = (itag == -1) ? -2 : 2;
    bool b1 = (itag == 1) ? true : false;

    /* evaluate mapping */
    ivec box = ptrnode->_box;
    ivec fbox = fcn(box);
    
    if (cst.isin(fbox)) { // inside
	ptrnode->_tag = itag;
	ptrnode->_b0 = false;
	ptrnode->_b1 = b1;
	
	return;
    }

    if (cst.isout(fbox)) { // outside
	// ptrnode->_tag = 0;
	// ptrnode->_b0 = true;
	// ptrnode->_b1 = false;
	return;
    }
	
    if (box.maxwidth() < eps) { // undetermined

	// std::cout << box.maxwidth() << ",  " << box << '\n';
	if (!inner) { // outer approximation
	    ptrnode->_tag = 1;  
	    ptrnode->_b0 = false;
	    ptrnode->_b1 = b1;
	}
	// else { // inner approximation

	//     ptrnode->_tag = 0;  
	//     ptrnode->_b0 = true;
	//     ptrnode->_b1 = false;
	// }

	return;
    }

    /* precision not satisfied: expansion */
    // ptrnode->_tag = utag;
    // ptrnode->_b0 = true;
    // ptrnode->_b1 = true;
    
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
    for (size_t i = 0; i < sp.leafcount(sp._root); ++i) {

	sivia(sp, lev[i], y, f, inner, itag, eps);
    }
    
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

	// case 0 :
	//     if (current->_box.isin(box))  // box inside current
	// 	return 0;
	//     if (!current->_box.isout(box))  // box intersect current
	// 	t0 = true;
	//     break;

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


void CSolver::pre_cntl(PtrVF pf,
		       std::stack<SPnode*> &l, std::stack<SPnode*> &l0,
		       std::stack<SPnode*> &l1, std::stack<SPnode*> &l2,
		       double evaleps /* 0.01 */) {

    // /*********** logging ***********/
    // std::fstream logfile;
    // logfile.open("cpre.log", std::ios::out | std::ios::ate | std::ios::app);
    // logfile << _fpiter << ":\n";
    // /*********** logging ***********/


    SPnode *current;
    ivec box;
    std::vector<ivec> fbox;

    short t, ut;

    while (!l.empty()) {

	current = l.top();
	l.pop();

	box = current->_box;  // leaves
	fbox = (*pf)(box);  // a vector of boxes corresponding to controls

	// /*********** logging ***********/
	// logfile << box << ", " << current->_tag << ", ";
	// for (int k = 0; k < fbox.size(); k++) {
	//     logfile << current->_cntl[k] << ", ";
	// }
	// /*********** logging ***********/

	/* update control info */
	t = 0;
	
	for (int u = 0; u < fbox.size(); ++ u ) {
	    
	    ut = paver_test(_ctlr, fbox[u]); // interval inclusion test

	    // /*********** logging ***********/
	    // logfile << u << ": " << fbox[u] << ut << '\n';
	    // /*********** logging ***********/

	    if (ut == 1) {

		if (_cpb->_workspace.isin(fbox[u])) {
		    
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
	

	// /*********** logging ***********/
	// logfile << "t=" << t << '\n';
	// /*********** logging ***********/

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
	}
	else {

	    current->_b0 = true;
	    current->_b1 = true;

	    if (box.maxwidth() < evaleps) {

		l2.push(current);
	    }
	    else {

		_ctlr.expand(current, bisect_axis(box));
		l.push(current->_left);
		l.push(current->_right);
	    }
	}  // end if (save new tag in b0 & b1)

	// /*********** logging ***********/
	// logfile << box << ", " << current->_split << ", ";
	// for (int k = 0; k < fbox.size(); k++) {
	//     logfile << current->_cntl[k] << ", ";
	// }
	// logfile << '\n';
	// /*********** logging ***********/

    }  // end while (loop all nodes in l)
    

    // logfile.close();
}


// void CSolver::compute_winsize() {

//     std::vector<double> lb = _cpb->_goal.getinf();
//     std::vector<double> ub = _cpb->_goal.getsup();

//     std::stack<SPnode*> stk;
//     stk.push(_ctlr._root);

//     SPnode *current;
//     std::vector<double> ubs, lbs;
//     while (!stk.empty()) {

//     	current = stk.top();
//     	stk.pop();

//     	if (current->_tag == 1) {

//     	    lbs = current->_box.getinf();
// 	    ubs = current->_box.getsup();
// 	    for (int i = 0; i < _cpb->_xdim; ++i) {

// 		lb[i] = lb[i] > lbs[i] ? lbs[i] : lb[i];
// 		ub[i] = ub[i] < ubs[i] ? ubs[i] : ub[i];
// 	    }
//     	}
//     	else {

// 	    if (!_ctlr.isleaf(current)) {
		
// 		if (current->_right)
// 		    stk.push(current->_right);

// 		if (current->_left)
// 		    stk.push(current->_left);
// 	    }
	    
//     	}
//     }

//     for (int i = 0; i < _cpb->_xdim; ++i) {

// 	if (_winsize < ub[i]-lb[i])
// 	    _winsize = ub[i]-lb[i];
//     }

//     std::cout << _winsize << '(winset width)\n';
// }


void CSolver::compute_winsize() {

    double v, vol = 0;
    std::vector<double> w(_cpb->_xdim);
    
    std::stack<SPnode*> stk;
    stk.push(_ctlr._root);

    SPnode *current;
    while (!stk.empty()) {

    	current = stk.top();
    	stk.pop();

	v = 1;
    	if (current->_tag == 1) {
	    w = current->_box.width();
	    for (int i = 0; i < _cpb->_xdim; ++i) {

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

    _winsize = pow(vol, 1.0/_cpb->_xdim);

    // std::cout << _winsize << "(winset radius)\n";
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

    	    if (current->_tag == 0)
    		l0.push(current);

    	    if (current->_tag == 1)
    		l1.push(current);

	    if (current->_tag == 2)
    		l2.push(current);
    	}
    	else {

    	    if (current->_right)
    		stk.push(current->_right);

    	    if (current->_left)
    		stk.push(current->_left);
    	}
    } // end while
    
}


void CSolver::inv_compute(std::stack<SPnode*> &l0,
			  std::stack<SPnode*> &l1,
			  std::stack<SPnode*> &l2,
			  int d,
			  const double eps, BISECT bs) {

    _bstype = bs;
    
    std::stack<SPnode*> l;
    int lold;
    bool stop = false;

    while (!stop) {

	// logging_iters("log.txt", _fpiter);
	
	++_fpiter[d];

	lold = l0.size() + l2.size();
	
	if (!l1.empty()) {

	    swap(l, l1);	
	    pre_cntl(_cpb->_vf, l, l0, l1, l2, eps);
	}

	stop = (l0.size() + l2.size()) <= lold;

	_ctlr.tagging(INNER);  //update the tags

    } // end while
    
}


void CSolver::reach_compute(std::stack<SPnode*> &l0,
			    std::stack<SPnode*> &l1,
			    std::stack<SPnode*> &l2,
			    int d,
			    const double er, BISECT bs,
			    const double ermin, bool vareps) {
    _bstype = bs;

    double eps;
    if (vareps)
	// std::cout << _winsize << '\n';
	eps = er * _winsize;
    else
	eps = er;
   
    std::stack<SPnode*> l;
    int lold;
    bool stop = false;
    bool sf = false;
    while (!stop && _fpiter[d] <= _maxiter) {

	++_fpiter[d];
	
	lold = l1.size();

	if (!l2.empty()) {

	    swap(l, l2);  // l=l2, l2=empty
	    pre_cntl(_cpb->_vf, l, l0, l1, l2, eps);  //l=empty, l012 fill
	}

	if (!l0.empty()) {

	    swap(l, l0);  // l=l0, l0=empty
	    pre_cntl(_cpb->_vf, l, l0, l1, l2, eps);  //l=empty, l012 fill
	}

	stop = l1.size() <= lold;
	_ctlr.tagging(INNER);

	if (vareps) {  // using adaptive precision

	    if (stop) {

		if (!sf)
		    sf = true;
		
		if (eps > ermin) {
		    eps /= 2;
		    stop = false;
		}
	    
	    } else {
	    
		if (!sf) {
		    compute_winsize();
		    eps = er * _winsize;
		}
	    
	    } // endif
	} // endif

	// if (vareps) {  // using adaptive precision

	//     if (stop) {
	// 	if (eps > ermin) {
	// 	    eps /= 2;
	// 	    stop = false;
	// 	}
	    
	//     } else {
	    
	// 	if (_winsize < _cpb->_workspace.maxwidth()) {
	// 	    compute_winsize();
	// 	    eps = er * _winsize;
	// 	}
	    
	//     } // endif
	// } // endif
	
    } // endwhile
}


bool CSolver::invariance_control(const double eps, BISECT bs) {

    std::stack<SPnode*> l0, l1, l2;
    init_leafque(l0, l1, l2);

    clock_t tb, te;
    tb = clock();
    
    inv_compute(l0, l1, l2, 0, eps, bs);
    
    te = clock();
    _timer = (float)(te - tb)/CLOCKS_PER_SEC;

    if (l1.empty())
	return false;
    else
	return true;
}


bool CSolver::reachability_control(const double eps, BISECT bs,
				   const double epsmin, bool vareps) {

    
    std::stack<SPnode*> l0, l1, l2;
    init_leafque(l0, l1, l2);
    
    clock_t tb, te;
    tb = clock();
    
    reach_compute(l0, l1, l2, 0, eps, bs, epsmin, vareps);  // _fpiter[0]
    
    te = clock();
    _timer = (float)(te - tb)/CLOCKS_PER_SEC;
    
    if (l1.empty())
	return false;
    else
	return true;
}


bool CSolver::reach_stay(const double ei, BISECT bsi,
			 const double er, BISECT bsr,
			 const double ermin, bool vareps) {

    double t;
    std::cout << "Start invariance control..." << '\n';
    if ( invariance_control(ei, bsi) ) {

	std::cout << "Start reachability control..." << '\n';
	
	t = _timer;
	bool r = reachability_control(er, bsr, ermin, vareps);
	_timer += t;
	return r;
	
    } else {

	return false;
    }
    
}


bool CSolver::buchi(BISECT bs, const double er,
		    const double ermin, bool vareps) {

    double eps;
    if (vareps)
	eps = er * _winsize;
    else
	eps = er;
    
    /* do the iterations */
    bool stop = false;
    int lold;
    clock_t tb, te;
    tb = clock();
    while (!stop) {
	
	++_fpiter[1];
	
	std::stack<SPnode*> l0, l1, l2;
	std::stack<SPnode*> x, l;
	init_leafque(l0, l1, l2);
	swap(l1, x); // x <- B, l1 <- empty.
	
	/* inner mu loop */
        reach_compute(l0, l1, l2, 0, er, bs, ermin, vareps); // l1 does not have B

	lold = l0.size() + l2.size();
	swap(l, x);
	pre_cntl(_cpb->_vf, l, l0, x, l2, eps); // x might decrease.
	
	stop = (l0.size() + l2.size()) <= lold;
	if (!stop) {
	    /* l1 <- 0, l2 <- 0, only x = 1 */
	    while (!l1.empty()) {
		l1.top()->_tag = 0;
		l1.pop();
	    }
	    while (!l2.empty()) {
		l2.top()->_tag = 0;
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


// bool CSolver::cobuchi(BISECT bs, const double eps) {
    
//     /* do the iterations */
//     bool stop = false;
//     int Gold1;
//     clock_t tb, te;
//     tb = clock();

//     std::stack<SPnode*> G10, G11, G12;
//     std::stack<SPnode*> G20, G21, G22;
//     std::stack<SPnode*> l;
//     init_leafque(G10, G21, G12);
    
//     while (!stop) {

// 	++ _fpiter[1];

// 	std::cout << _fpiter[1] << ": ";
// 	/* inner nu loop */
// 	_fpiter[0] = 0;
// 	if (!G21.empty()) {

// 	    if (!(G20.empty() && G22.empty()))
// 		_ctlr.tagging(INNER);  // mark (Z U G) to be 1

// 	    inv_compute(G20, G21, G22, 0, eps, bs);
// 	    std::cout << _fpiter[0] << ", ";

// 	    /* G2<- G20 U G22, G2 starts from empty */
// 	    std::stack<SPnode*> G2;
// 	    while (!G20.empty()) {
		
// 		G20.top()->_tag = 1;
// 		G2.push(G20.top());
// 		G20.pop();
// 	    }
// 	    while (!G22.empty()) {
		
// 		G22.top()->_tag = 1;
// 		G2.push(G22.top());
// 		G22.pop();
// 	    }

// 	    std::cout << G2.size() << '\n';
// 	    swap(G2, G21);
// 	}

// 	/* compute pre(Y) /\ G1 */
// 	Gold1 = G11.size();
// 	/* G1<- G10 U G12 */
// 	if (!G12.empty()) {
	    
// 	    swap(l, G12);
// 	    pre_cntl(_cpb->_vf, l, G10, G11, G12, eps);
// 	}
// 	if (!G10.empty()) {
	    
// 	    swap(l, G10);
// 	    pre_cntl(_cpb->_vf, l, G10, G11, G12, eps);
// 	}
// 	_ctlr.tagging(INNER);

// 	stop = G11.size() <= Gold1;
	
//     }

//     te = clock();
//     _timer = (float)(te - tb)/CLOCKS_PER_SEC;

//     return true;
// }


bool CSolver::cobuchi(const double ei, BISECT bsi,
		      const double er, BISECT bsr,
		      const double ermin, bool vareps) {
    
    /* do the iterations */
    double eps;  // absolute epsilon for repeated reach computationdouble eps;
    if (vareps)
	// std::cout << _winsize << '\n';
	eps = er * _winsize;
    else
	eps = er;
    
    bool stop = false;
    bool sf = false;
    int Gold1;
    
    clock_t tb, te;
    tb = clock();

    std::stack<SPnode*> G10, G11, G12;
    std::stack<SPnode*> G20, G21, G22;
    std::stack<SPnode*> l;
    init_leafque(G10, G21, G12);

    std::cout << "<outer iter>:<# of inner iters>,<current precision parameter>\n";
    
    while (!stop && _fpiter[1] <= _maxiter) {

	++ _fpiter[1];

	std::cout << _fpiter[1] << ": ";
	/* inner nu loop */
	_fpiter[0] = 0;
	if (!G21.empty()) {

	    if (!(G20.empty() && G22.empty()))
		_ctlr.tagging(INNER);  // mark (Z U G) to be 1

	    inv_compute(G20, G21, G22, 0, ei, bsi);
	    std::cout << _fpiter[0] << ", ";

	    /* G2<- G20 U G22, G2 starts from empty */
	    std::stack<SPnode*> G2;
	    while (!G20.empty()) {
		
		G20.top()->_tag = 1;
		G2.push(G20.top());
		G20.pop();
	    }
	    while (!G22.empty()) {
		
		G22.top()->_tag = 1;
		G2.push(G22.top());
		G22.pop();
	    }

	    // std::cout << "(G2 size: " << G2.size() << ")\n";
	    swap(G2, G21);
	} else {
	    std::cout << _fpiter[0] << ", ";
  	}

	/* compute pre(Y) /\ G1 */
	_bstype = bsr;
	Gold1 = G11.size();
	/* G1<- G10 U G12 */
	if (!G12.empty()) {
	    
	    swap(l, G12);
	    pre_cntl(_cpb->_vf, l, G10, G11, G12, eps);
	}
	if (!G10.empty()) {
	    
	    swap(l, G10);
	    pre_cntl(_cpb->_vf, l, G10, G11, G12, eps);
	}

	std::cout << eps << '\n';

	stop = G11.size() <= Gold1;
	_ctlr.tagging(INNER);

	if (vareps) {  // using adaptive precision

	    if (stop) {

		if (!sf)
		    sf = true;
		
		if (eps > ermin) {
		    eps /= 2;
		    stop = false;
		}
	    
	    } else {
	    
		if (!sf) {
		    compute_winsize();
		    eps = er * _winsize;
		}
	    
	    } // endif
	} // endif
	
    } // endwhile

    te = clock();
    _timer = (float)(te - tb)/CLOCKS_PER_SEC;

    return true;
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
	// std::cout << "Maximal number of iterations: " << _maxiter << '\n';
	
	std::cout << "Time of solving: "
		  << _timer << '\n';
    }
}


void CSolver::write2mat_controller(const char *filename) {

    assert( !_ctlr.isempty() );

    _ctlr.write2mat_leaves(_ctlr._root, filename);
}


void CSolver::write2txt_controller(const char *filename) {

    assert( !_ctlr.isempty() );

    _ctlr.write2txt_leaves(_ctlr._root, filename);
}


void CSolver::serialize_controller(const char *filename) {

    if (_ctlr.isempty()) {

	std::cout << "No Controller Generated.\n";
	return;
    }

    MATFile *pmat;
    pmat = matOpen(filename, "w");
    if (pmat == NULL) {
	std::cout << "Error creating file" << filename;
	return;
    }

    size_t H = _ctlr.height(_ctlr._root);
    // mwSize N = _ctlr.nodecount(_ctlr._root);
    mwSize L = _ctlr.leafcount(_ctlr._root);
    mwSize M = pow(2,H)-1;
    
    int nx = _cpb->_xdim;
    int nu = _cpb->_vf->_ugrid._nv;
    
    /* define outputs */
    mxArray *mats = mxCreateDoubleMatrix(M, 2*nx+1, mxREAL);  // Btree
    double *ptrs = mxGetPr(mats);
    mxArray *matu = mxCreateLogicalMatrix(L, nu);  // control table
    mxLogical *ptru = mxGetLogicals(matu);
    mxArray *mati = mxCreateDoubleMatrix(L, 2, mxREAL);  // indices & tags of control entries
    double *ptri = mxGetPr(mati);

    /* write data to buffers using stacks */
    std::queue<SPnode*> qnode;
    std::queue<size_t> qid;
    qnode.push(_ctlr._root);
    qid.push(0);

    std::vector<double> lower, upper;
    SPnode *current;
    size_t i;
    size_t lc = 0;
    while(!qnode.empty()) {

	current = qnode.front();
	qnode.pop();
	i = qid.front();
	qid.pop();

	ptrs[i] = current->_split;  
	
	lower = current->_box.getinf();  
	upper = current->_box.getsup();
	for (int j = 0; j < nx; ++j) {
		
	    ptrs[i + (2*j+1)*M] = lower[j];
	    ptrs[i + (2*j+2)*M] = upper[j];
	}

	if (_ctlr.isleaf(current)) {
	    
	    for (int k = 0; k < nu; ++k) { // control array

		ptru[lc + k*L] = current->_cntl[k];
	    }
	    ptri[lc] = i + 1;
	    ptri[lc + L] = current->_tag;

	    ++lc;
	}
	else {

	    if (current->_left != NULL) {

		qnode.push(current->_left);
		qid.push(2*i+1);
	    }

	    if (current->_right != NULL) {

		qnode.push(current->_right);
		qid.push(2*i+2);
	    }
	} // end if
    } // end while

    assert(i < M);
    assert(lc-1 < L);

    /* write to variables in .mat */
    if (matPutVariable(pmat, "ctree", mats) ||
	matPutVariable(pmat, "cindex", mati) ||
	matPutVariable(pmat, "cvalue", matu) != 0) {
	
	std::cout << "Error outputting variables.\n";
        return;
    }
    
    mxDestroyArray(mats);
    mxDestroyArray(matu);
    mxDestroyArray(mati);
    
    if (EOF == matClose(pmat)) {
        return;
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

	    for (int i = 0; i < _cpb->_xdim; ++ i)
		logfile << lower[i] << ' ' << upper[i] << ' ';

	    // for (int u = 0; u < node->_cntl.size(); ++ u)
	    // 	logfile << ' ' << current->_cntl[u];

	    logfile << '\n';
	}

	if (current->_left)
	    stk.push(current->_left);

	if (current->_right)
	    stk.push(current->_right);

	// std::cout << std::endl; //logging
    }
    logfile << '\n';
    
    logfile.close();
}

    
} // namespace rocs
