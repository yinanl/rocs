/**
 *  patcher.cpp
 *
 *  Implementation of the patcher class.
 *
 *  Created by Yinan Li on Oct. 17, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include "patcher.h"
#include <fstream>


namespace rocs {

    void Patcher::initialize_winning_graph(HEAD &winprod) {
	SCOPE *winset = winprod.nts_product.scope;
	NODE_PRODUCT *graph = winprod.nts_product.graph;
	_na = winprod.action;
	_nwin = winprod.ctrlr.wp;

	/* Fill up _idmap */
	_idmap.resize(winprod.nts_product.n, -1); //np
	size_t idwin = 0;
	for(SCOPE *iter = winset; iter->next; ++idwin) {
	    iter += iter->next;
	    _idmap[iter->x] = idwin;
	}
	assert(idwin == _nwin);

	/* Fill up _winfts */
	_winfts.init(_nwin, _na);
	NODE_PRODUCT *pnode;
	size_t ntrans = 0;
	size_t ntrans1 = 0; //debug only
	for(SCOPE *iter_winset = winset; iter_winset->next;) {
	    iter_winset += iter_winset->next;
	    /* Fill up _winfts._npost (or the outdeg table) */
	    pnode = graph + iter_winset->x;  //point to a post node
	    // if(pnode->outdeg_reset) {
	    for(OUTDEG_RESET *iter = pnode->outdeg_reset; iter->next;) {
		iter += iter->next;
		_winfts._npost[_idmap[iter_winset->x]*_na+iter->a] = iter->outdeg;
		ntrans += iter->outdeg;  //accumulate the # of transitions
	    }
	    // }
	    /* Fill up _winfts._npre */
	    if(pnode->in_edge) { //A winning node can have no in edges.
		for(EDGE *edgein = pnode->in_edge; edgein->next;) {
		    edgein += edgein->next;
		    ++_winfts._npre[_idmap[iter_winset->x]*_na+edgein->a];
		    ++ntrans1; //debug only
		}
	    }
	}
	assert(ntrans == ntrans1);

	/* Fill up _winfts._ptrpre and _ptrpost */
	size_t sum_post = 0;
	size_t sum_pre = 0;
	for(size_t i = 0; i < _nwin; ++i) {
	    for(size_t j = 0; j < _na; ++j) {
		if(_winfts._npost[i*_na+j]) {
		    _winfts._ptrpost[i*_na+j] = sum_post;
		    sum_post += _winfts._npost[i*_na+j];
		}
		if(_winfts._npre[i*_na+j]) {
		    _winfts._ptrpre[i*_na+j] = sum_pre;
		    sum_pre += _winfts._npre[i*_na+j];
		}
	    }
	}
	/* Fill up _winfts._idpre and _idpost */
	_winfts._idpost.resize(ntrans);
	_winfts._idpre.resize(ntrans);
	std::vector<int> incpre(_nwin*_na, 0);
	std::vector<int> incpost(_nwin*_na, 0);
	size_t idxpre, idxpost;
	for(SCOPE *iter=winset; iter->next;) {
	    iter += iter->next;
	    pnode = graph + iter->x;
	    if(pnode->in_edge) {
		for(EDGE *edgein=pnode->in_edge; edgein->next;) {
		    edgein += edgein->next;
		    /* Fill up _idpre */
		    idxpre = _idmap[iter->x]*_na + edgein->a;
		    _winfts._idpre[_winfts._ptrpre[idxpre] + incpre[idxpre]++] = edgein->x;
		    /* Fill up _idpost */
		    idxpost = _idmap[edgein->x]*_na + edgein->a;
		    _winfts._idpost[_winfts._ptrpost[idxpost] + incpost[idxpost]++] = iter->x;
		}
	    }
	}
	assert(incpre == _winfts._npre);
	assert(incpost == _winfts._npost);

	/* Fill up _reachstep by backward reachability */
	_reachstep.resize(_nwin, -1);
	std::vector<size_t> queue(_nwin);
	size_t p = 0;  //point to the tail of queue
	/* Put accepting nodes into a queue */
	for(SCOPE *iter=winset; iter->next;) {
	    iter += iter->next;
	    if(winprod.nts_product.acc[iter->x]) {
		queue[p] = iter->x;
		_reachstep[_idmap[iter->x]] = 0;
		++p;
	    }
	}
	size_t outidx, e = p;
	size_t step = 1;
	for(size_t f = 0; f < e; ++f) {
	    if(f == p) {
		++step;
		p = e;
	    }
	    /* Loop all the inedges of node queue[f] */
	    pnode = graph+queue[f];
	    if(pnode->in_edge) {
		for(EDGE *it = pnode->in_edge; it->next;) {
		    it += it->next;
		    outidx = _idmap[it->x]*_na + it->a;
		    if(incpost[outidx]) { //out degree>0(can be 0 or positive)
			--incpost[outidx];
			if(!incpost[outidx]) { //out degree=0
			    if(_reachstep[_idmap[it->x]] < 1) {
				//The reachstep of it->x hasn't been marked yet.
				//Possible that reachstep of a node can be marked multiple times.
				//The reachstep should record the smallest one.
				if(_reachstep[_idmap[it->x]] < 0) { //it->x isn't in queue
				    queue[e] = it->x; //push it->x to queue
				    ++e; //move e forward
				}
				/* step>=1 for all winning states including accepting ones */
				_reachstep[_idmap[it->x]] = step;
			    }
			}
		    }
		} //end for loop of all the inedges of queue[f]
	    }
	} //end BFS for backward reachability
	std::vector<int> test(_nwin*_na, 0);
	assert(incpost == test); //incpost should be all 0, since every transition in _winfts is winning

	/* Fill up _encode[n0xn1]->np, _decode[np]->n0xn1 */
	size_t n0 = winprod.nts_pre.n; //the # of states in the original FTS
	size_t n1 = winprod.dba.n; //the # of DBA nodes
	size_t np = winprod.nts_product.n;
	size_t x, q;
	_encode.resize(n0*n1, -1);
	_decode.resize(np);
	for(size_t i = 0; i < np; ++i) {
	    x = *(winprod.decode2+i) / n1;
	    q = *(winprod.decode2+i) % n1;
	    _decode[i] = *(winprod.decode1+x)*n1 + q;
	    _encode[*(winprod.decode1+x)*n1 + q] = i;
	}

    } //end initialize_winning_graph


    int Patcher::solve_local_reachability(size_t v, int h,
    					  boost::dynamic_bitset<> &safecnt) {
    	_ctlr = -1; //re-initialization
    	boost::dynamic_bitset<> vers(_nwin, false); //default false
    	std::vector<size_t> queue(_nwin);
    	queue[0] = _encode[v];
    	vers[_idmap[queue[0]]] = true;

    	boost::dynamic_bitset<> visited(_nwin, false);
    	visited[_idmap[queue[0]]] = true;
    	/**
    	 * Collect target nodes by forward propagation from v:
    	 * The nodes propagated at the last step are treated as sinks,
    	 * e.g., their outgoing edges are not included.
    	 **/
    	size_t f, e, p; //f points the top, e the end, p the end of one horizon
    	size_t row, idx, pidstart;
    	int i, ntrans = 0;
    	// f = 0;
    	// e = p = 1;
    	for(i = 0, f = 0, e = p = 1; f < e; ++f) { //h is the propagation steps
    	    if(f == p) {
    		++i;
    		p = e;
    	    }
    	    if(i >= h)
    		break;
    	    row = _idmap[queue[f]];
    	    for(size_t j = 0; j < _na; ++j) {
    		if(!f) { //only for the initial node _encode[v]
    		    if(!safecnt[j]) { //disable unsafe actions
    			continue;
    		    }
    		}
    		pidstart = _winfts._ptrpost[row*_na+j];
    		for(size_t pid = pidstart; pid < pidstart+_winfts._npost[row*_na+j]; ++pid) {
    		    idx = _winfts._idpost[pid];
    		    if(!vers[_idmap[idx]]) { //if not in queue, enqueue
    			queue[e] = idx;
    			vers[_idmap[idx]] = true; //mark inqueue
    			++e;
    		    }
    		    ++ntrans;
    		}
    	    }
    	    visited[row] = true;
    	}
    	// for(;f < e; ++f) {//calculate the # of transitions in subgraph
    	//     row = _idmap[queue[f]];
    	//     for(size_t j = 0; j < _na; ++j) {
    	// 	pidstart = _winfts._ptrpost[row*_na+j];
    	// 	for(size_t pid = pidstart; pid<pidstart+_winfts._npost[row*_na+j]; ++pid) {
    	// 	    idx = _winfts._idpost[pid];
    	// 	    if(vers[_idmap[idx]])
    	// 		++ntrans;
    	// 	}
    	//     }
    	// }

    	// /********** Logging ***********/
    	// std::cout << "# of post transitions in the subgraph: " << ntrans << '\n';
    	// std::cout << "# of nodes in the subgraph: " << e << '\n';
    	// /********** Logging ***********/

    	/* Construct a subgraph with inedge info */
    	std::vector<long long> keyidx(_nwin, -1); //index mapping from np-based to e-based
    	std::vector<int> outdeg(e*_na, 0);
    	std::vector<size_t> npre_sub(e*_na, 0);
    	std::vector<size_t> ptrpre_sub(e*_na);
    	std::vector<size_t> idpre_sub(ntrans);
    	int idtrans = 0;
    	int minstep = _reachstep[_idmap[queue[0]]];
    	int maxstep = minstep;
    	for(f = 0; f < e; ++f) { //loop all the expanded nodes in queue
    	    row = _idmap[queue[f]]; //the row in _winfts._ptrpre table
    	    /* Fill up the index mapping keyidx */
    	    keyidx[row] = f;
    	    for(size_t j = 0; j < _na; ++j) { //loop all the actions
    		/* Fill up npre_sub, ptrpre_sub, idpre_sub */
    		pidstart = _winfts._ptrpre[row*_na+j];
    		for(size_t pid = pidstart; pid < pidstart+_winfts._npre[row*_na+j]; ++pid) {
    		    idx = _winfts._idpre[pid];
    		    if(idx==queue[0] && !safecnt[j]) {
    			continue;
    		    } else {
    			if(visited[_idmap[idx]]) { //idx is in queue and expanded
    			    if(!npre_sub[f*_na+j])
    				ptrpre_sub[f*_na+j] = idtrans; //only store the ptr of the first pre
    			    ++npre_sub[f*_na+j];
    			    idpre_sub[idtrans] = idx;
    			    ++idtrans;
    			}
    		    }
    		}
    		/* Fill up outdeg for queue[1 ~ p-1] */
    		if(f > 0 && f < p) { //for queue[p->e], outdeg = 0, no need to fill
    		    outdeg[f*_na+j] = _winfts._npost[row*_na+j];
    		}
    	    }
    	    /* Search for the least and greatest reachstep */
    	    if(_reachstep[_idmap[queue[f]]] < minstep)
    		minstep = _reachstep[_idmap[queue[f]]];
    	    if(_reachstep[_idmap[queue[f]]] > maxstep)
    	    	maxstep = _reachstep[_idmap[queue[f]]];
    	}
    	/* Fill up outdeg for f=0 */
    	for(size_t j = 0; j < _na; ++j) {
    	    if(safecnt[j])
    		outdeg[j] = _winfts._npost[_idmap[queue[0]]*_na+j];
    	}
    	if(!(minstep < _reachstep[_idmap[queue[0]]])) {
    	    std::cout << "Patcher::solve_local_reachability: Cannot make progress in planner.\n";
    	    return 0;
    	}
    	assert(ntrans == idtrans);

    	// /********** Logging ***********/
    	// std::cout << "# of pre transitions: " << idtrans << '\n';
    	// std::cout << "minstep: " << minstep << ", maxstep: " << maxstep << '\n';
    	// /********** Logging ***********/

    	/**
    	 * Solve for local controller by backward reachability:
    	 * pick the nodes with smallest value in _reachstep as targets,
    	 * at least smaller than _reachstep[_idmap[queue[0]]].
    	 **/
    	// /********** Logging ***********/
    	// std::cout << "Start local backward reachability.\n";
    	// /********** Logging ***********/
    	size_t ib, ie;
    	std::vector<size_t> rque(e);
    	boost::dynamic_bitset<> reach(e, false);
    	for(f = 0, ie = 0; f < e; ++f) { //enqueue targets
    	    if(_reachstep[_idmap[queue[f]]]>=minstep &&
    	       _reachstep[_idmap[queue[f]]]<=minstep+h/3) {
    	    // if(_reachstep[_idmap[queue[f]]]>=minstep &&
    	    //    _reachstep[_idmap[queue[f]]]<maxstep) { //maxstep=_reachstep[queue[0]]
    		reach[f] = true;
    		rque[ie] = queue[f];
    		++ie;
    	    }
    	} //end targets enqueue

    	/********** Logging ***********/
    	// std::cout << "# of target nodes: " << ie << '\n';
    	std::cout << "minstep= " << minstep << ", "
    		  << "maxstep= " << maxstep << ", "
    		  << "# of targets= " << ie << '\n';
    	/********** Logging ***********/

    	/* Solve by BFS */
    	size_t rpre, rpost;
    	for(ib = 0; ib < ie; ++ib) { //ie=# of elements in rque from the last loop
    	    rpost = keyidx[_idmap[rque[ib]]]; //the row in npre_sub, ptrpre_sub
    	    for(size_t j = 0; j < _na; ++j) {
    		/* Loop all the inedges of node rque[ib] with control j */
    		for(size_t i = 0; i < npre_sub[rpost*_na+j]; ++i) {
    		    idx = idpre_sub[ptrpre_sub[rpost*_na+j]+i];
    		    rpre = keyidx[_idmap[idx]]; //the row of post node
    		    if(outdeg[rpre*_na+j]) {
    			--outdeg[rpre*_na+j];
    			if(!outdeg[rpre*_na+j]) { //out degree = 0
    			    if(idx == queue[0] && _ctlr < 0) {
    				_ctlr = j;
    				break;
    			    }
    			    if(!reach[rpre]) { //pre node not in reach set
    				reach[rpre] = true;
    				rque[ie] = idx; //enque pre node
    				++ie;
    			    }
    			}
    		    }
    		} //end loop all the inedges
    		if(_ctlr >= 0)
    		    break;
    	    } //end loop all actions
    	    if(_ctlr >= 0)
    		break;
    	} //end loop rque

    	if(_ctlr < 0) //didn't reach back to queue[0]
    	    return 0;

    	return 1;
    } //end solve_local_reachability


    std::vector< std::pair<size_t, int> >
    Patcher::solve_local_reachavoid(size_t v, size_t N,
				    std::vector<long long> &avoid,
				    std::vector<long long> &target,
				    boost::dynamic_bitset<> &safecnt) {
    	assert(N <= target.size());

    	/* Mark avoid and target nodes
	 * Assumption: target and avoid may contain nodes not on the _winfts
	 */
    	boost::dynamic_bitset<> isavoid(_nwin, false);
    	boost::dynamic_bitset<> isgoal(_nwin, false);
    	for(auto &item : avoid) {
	    if(_idmap[item]>=0)
		isavoid[_idmap[item]] = true;
	}
    	for(auto &item : target) {
	    if(_idmap[item]>=0)
		isgoal[_idmap[item]] = true;
	}

    	/* Create a copy of _winfts._npost */
    	std::vector<int> win_npost(_winfts._npost);
    	/* Assign corresponding outdegs of predecessors of avoid to 0 */
	long long idx, row;
    	size_t pidstart;
    	for(auto &item : avoid) {
    	    row = _idmap[item];
	    if(row>=0) {
		for(size_t j = 0; j < _na; ++j) {
		    pidstart = _winfts._ptrpre[row*_na+j];
		    for(size_t pid=pidstart;pid<pidstart+_winfts._npre[row*_na+j];++pid) {
			idx = _winfts._idpre[pid];
			win_npost[_idmap[idx]*_na+j] = 0;
		    }
		}
	    }
    	}
	
	
    	/**
    	 * Forward propagation from v until N target nodes are collected.
    	 **/
    	/* Create a queue for the subgraph nodes */
    	std::vector<size_t> queue(_nwin);
    	queue[0] = _encode[v];
	boost::dynamic_bitset<> inque(_nwin, false); //mark if a node is in queue
    	inque[_idmap[queue[0]]] = true;
	boost::dynamic_bitset<> visited(_nwin, false); //mark the node that has been expanded
	/********** Logging ***********/
    	std::cout << "# of winning nodes: " << _nwin << '\n';
    	std::cout << "Current node in the winning graph: " << queue[0] << '\n';
    	/********** Logging ***********/
    	
    	/* Enqueue all safe nodes on the _winfts. */
    	size_t b, e, p; //b(begin), e(end)
    	// size_t row, idx, pidstart;
    	size_t i, ntrans = 0;
    	for(i = 0, b = 0, e = 1; b < e; ++b) {
    	    if(i > N) {
    		p = b;
    		break;
    	    }
    	    row = _idmap[queue[b]];
    	    if(isgoal[row]) {
		// std::cout << row << ", is in queue? " << inque[row] << '\n';
    		++i; //count target nodes
    		// continue; // no need to expand target node (a target is a sink)
    	    }
	    visited[row] = true; //mark visited (or expanded)
    	    for(size_t j = 0; j < _na; ++j) {
    		if(!b) { //only for the initial node _encode[v]
    		    if(!safecnt[j]) { //disable unsafe actions
    			continue;
    		    }
    		}
    		pidstart = _winfts._ptrpost[row*_na+j];
    		// for(size_t pid = pidstart; pid < pidstart+_winfts._npost[row*_na+j]; ++pid) {
    		//     idx = _winfts._idpost[pid];
    		//     if(isavoid[_idmap[idx]]) {
    		// 	win_npost[row*_na+j] = 0; //disable control j
    		// 	break;
    		//     }
    		// }
    		for(size_t pid=pidstart; pid<pidstart+win_npost[row*_na+j];++pid) {
    		    idx = _winfts._idpost[pid];
    		    if(!inque[_idmap[idx]]) { //if not in queue, enqueue
    			queue[e] = idx;
    			inque[_idmap[idx]] = true; //mark inqueue
    			++e;
    		    }
    		    ++ntrans;
    		}
    	    }
    	}
    	if(b==e)
    	    p = e;

	/********** Logging ***********/
    	std::cout << "# of post transitions in the subgraph: " << ntrans << '\n';
    	std::cout << "# of nodes in the subgraph: " << e << '\n';
	std::cout << "# of target nodes: " << i << '\n';
	std::cout << "i= " << i << ", b=" << b << ", p=" << p << ", e=" << e << '\n';
    	/********** Logging ***********/

	std::vector< std::pair<size_t, int> > ctlrtable(e);
	// _replan.resize(e);
	if(i < 10) {
	    std::cout << "Patcher::solve_local_reachavoid: not enough target points propagated.\n";
	    return std::move(ctlrtable);
	}

    	/* Construct a subgraph with inedge info */
	// /********** Logging ***********/
	// logger.open("logs_subpre.txt", std::ios::out);
    	// /********** Logging ***********/
    	std::vector<long long> keyidx(_nwin, -1); //index mapping from np-based to e-based
    	std::vector<int> outdeg(e*_na, 0);
    	std::vector<size_t> npre_sub(e*_na, 0);
    	std::vector<size_t> ptrpre_sub(e*_na);
    	std::vector<size_t> idpre_sub(ntrans);
    	size_t idtrans = 0;
    	for(b = 0; b < e; ++b) {
    	    row = _idmap[queue[b]]; //the row in _winfts._ptrpre table
    	    keyidx[row] = b;
    	    for(size_t j = 0; j < _na; ++j) { //loop all the actions
		/********** Logging ***********/
		// std::cout << "b= " << b << ", node id= " << queue[b] << ", ctlr id= " << j
		// 	  << ", idtrans= " << idtrans << '\n';
		// if(b == 64718 && j == 16) {
		//     std::cout << "b=64781, j=16: \nNo of pres= " << _winfts._npre[row*_na+j]
		// 	      << ", size of idpre_sub(before j)= " << idpre_sub.size() << '\n';
		// }
		/********** Logging ***********/
    		/* Fill up npre_sub, ptrpre_sub, idpre_sub */
    		pidstart = _winfts._ptrpre[row*_na+j];
    		for(size_t pid=pidstart;pid<pidstart+_winfts._npre[row*_na+j];++pid) {
    		    idx = _winfts._idpre[pid];
    		    if(idx==queue[0] && !safecnt[j]) {
    			continue;
    		    } else {
    			if(visited[_idmap[idx]] && win_npost[_idmap[idx]*_na+j]) { //idx is in queue
    			    if(!npre_sub[b*_na+j])
    				ptrpre_sub[b*_na+j] = idtrans; //only store the ptr of the first pre
    			    ++npre_sub[b*_na+j];
    			    idpre_sub[idtrans] = idx;
    			    ++idtrans;
			    // /********** Logging ***********/
			    // logger << idx << ',' << j << ',' << queue[b] << '\n';
			    // /********** Logging ***********/
    			}
    		    }
    		}
		/********** Logging ***********/
		// if(b == 64718 && j == 16) {
		//     std::cout << "size of idpre_sub(after j)= " << idpre_sub.size() << '\n';
		// }
		/********** Logging ***********/
    		/* Fill up outdeg for queue[1 ~ p-1] */
    		if(!isgoal[row] && b > 0 && b < p) { //for queue[p~e], outdeg = 0, no need to fill
    		    outdeg[b*_na+j] = win_npost[row*_na+j];
    		}
    	    }
    	}

    	/* Fill up outdeg for b=0 (or v) */
    	for(size_t j = 0; j < _na; ++j) {
    	    if(safecnt[j])
    		outdeg[j] = win_npost[_idmap[queue[0]]*_na+j];
    	}

	/********** Logging ***********/
	// logger.close();
    	std::cout << "# of pre transitions: " << idtrans << '\n';
    	/********** Logging ***********/
	assert(ntrans == idtrans);

    	/**
    	 * Solve for local controller by backward reachability.
    	 **/
    	size_t ib, ie;
    	std::vector<size_t> rque(e);
    	boost::dynamic_bitset<> reach(e, false);
    	for(b = 0, ie = 0; b < e; ++b) {
    	    if(isgoal[_idmap[queue[b]]]) { //enqueue target nodes
    		reach[b] = true;
    		rque[ie] = queue[b];
    		++ie;
    	    }
    	}
    	/* Solve by BFS */
    	size_t rpre, rpost;
    	for(ib = 0; ib < ie; ++ib) {
    	    rpost = keyidx[_idmap[rque[ib]]];
    	    for(size_t j = 0; j < _na; ++j) {
    		for(size_t i = 0; i < npre_sub[rpost*_na+j]; ++i) {
    		    idx = idpre_sub[ptrpre_sub[rpost*_na+j]+i];
    		    rpre = keyidx[_idmap[idx]];
    		    if(outdeg[rpre*_na+j]) {
    			--outdeg[rpre*_na+j];
    			if(!outdeg[rpre*_na+j]) {
    			    if(!reach[rpre]) { //pre node not in reach set
    				reach[rpre] = true;
    				rque[ie] = idx; //enque pre node
    				ctlrtable[rpre] = std::make_pair(_decode[idx], j);
				// _replan[rpre] = std::make_pair(_decode[idx], j); //save (n0xn1)-based node id
    				++ie;
    			    }
    			}
    		    }
    		}
    	    }
    	}

    	return std::move(ctlrtable);
	// return _replan[0].first;
    }//solve_local_reachavoid


    std::vector<long long> Patcher::replan_region(const rocs::Rn &xo, const double ro,
						  rocs::grid &x_grid,
						  const rocs::UintSmall &q,
						  const rocs::UintSmall &nNodes) {
	/* collect the grid points intersect with the box */
	rocs::ivec box{rocs::interval(xo[0]-ro, xo[0]+ro),
		       rocs::interval(xo[1]-ro, xo[1]+ro),
		       x_grid._bds[2]};
	// std::cout << "The box enclose the region: " << box << '\n';
	std::vector<size_t> ids = x_grid.subset(box, 0, 0);
	std::vector<long long> ix(ids.begin(), ids.end());
			    
	/* determine the grid cell intersect with the r-ball around xo */
	double xl, xr, yr, yl, xsqr, ysqr;
	double rx = x_grid._gw[0]/2.0;
	double ry = x_grid._gw[1]/2.0;
	std::vector<double> x(x_grid._dim);
	for(auto &id : ix) {
	    x_grid.id_to_val(x, id);
	    // std::cout << x[0] << ' ' << x[1] << ' ' << x[2] <<'\n';
	    xl = x[0]-xo[0] - rx;
	    xr = x[0]-xo[0] + rx;
	    yl = x[1]-xo[1] - ry;
	    yr = x[1]-xo[1] + ry;
	    xsqr = (xr*xr) > (xl*xl) ? (xl*xl) : (xr*xr);
	    ysqr = (yr*yr) > (yl*yl) ? (yl*yl) : (yr*yr);
	    if(xsqr+ysqr <= ro*ro) {
		// std::cout << id <<',';
		id = _encode[id*nNodes+q]; //n0xn1->np
	    } else
		id = -1;
	}
			    
	return std::move(ix);
    }//replan_region

    
}//namespace rocs
