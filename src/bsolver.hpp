/**
 *  bsolver.hpp
 *
 *  A class for solving a Buchi game on product system of an NTS and a DBA.
 *
 *  Created by Yinan Li on Oct. 26, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _bsolver_h_
#define _bsolver_h_


#include <climits>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <memory>

#include "config.h"
#include "abstraction.hpp"
#include "hdf5io.h"

extern "C" {
#include "buchi.h"
}


namespace rocs {
  
    /**
     * A control synthesis solver class over a (finite) transition system.
     */
    class BSolver
    {
    public:
	HEAD _sol;

	BSolver() {initialization(&_sol);}
	BSolver(const BSolver&) = delete;
	BSolver(BSolver&&) = delete;
	BSolver& operator=(const BSolver&) = delete;
	BSolver& operator=(BSolver&&) = delete;
  
	~BSolver() {free_memory(&_sol, 1);}

	/**
	 * Construct the DBA struct in _sol.
	 */
	void construct_dba(int nAP, int nNodes, int q0,
			   std::vector<rocs::UintSmall> &acc,
			   std::vector<std::vector<rocs::UintSmall>> arrayM);

	/**
	 * Construct a graph with in-going edges in _sol.
	 */
	template<typename S>
	void load_abstraction(abstraction<S> &abst);

	/**
	 * Take product of the forward graph and the DBA.
	 */
	template<typename S>
	void generate_product(abstraction<S> &abst);

	/**
	 * Solve the buchi game on the product.
	 */
	void solve_buchigame_on_product() {
	    buchi_and_controller(&_sol);
	}

	/**
	 * Write controller to a file
	 */
	void write_controller_to_txt(char *filename) {
	    write_controller(&_sol, filename);
	}
	// void write_controller_to_h5(const std::string filename) {
	//     h5FileHandler wtr(filename, H5F_ACC_TRUNC);
	//     /* Write W_X0: the list of winning nodes */
	// }
    };


    /**
     * Construct the DBA struct in _sol: assign values to
     * - n: # of DBA states.
     * - k: # of atomic propositions.
     * - ini: the initial DBA state.
     * - acc: an array (k elements) of BOOL denoting accepting or not.
     * - q_prime: the lookup table of size (2^k x n), recording transition relation.
     * - q_prime_size: the total # of elements in q_prime (2^k x n).
     */
    inline void BSolver::construct_dba(int nAP, int nNodes, int q0,
				       std::vector<rocs::UintSmall> &acc,
				       std::vector<std::vector<rocs::UintSmall>> arrayM) {
	DBA *dba = &(_sol.dba);
	
	dba->k = nAP;
	dba->n = nNodes;
	dba->ini = q0;
	dba->acc = (BOOL*)calloc(nNodes, sizeof(BOOL)); // initialize to zero
	for(rocs::UintSmall i = 0; i < acc.size(); ++i)
	    *(dba->acc + acc[i]) = 1;

	/* Assign q_prime */
	size_t R = arrayM.size();  // R = nNodes
	size_t C = arrayM[0].size();  // C = nProps = 2^nAP
	dba->q_prime_size = (int)std::pow(2,nAP) * nNodes;
	dba->q_prime = (int*)malloc(dba->q_prime_size * sizeof(int));
	for(size_t i = 0; i < R; ++i)
	    for(size_t j = 0; j < C; ++j)
		*(dba->q_prime + j*R + i) = arrayM[i][j];
    } //end construct_dba
    

    /**
     * Construct a graph with in-going edges in _sol:
     * - nts_pre.n: # of graph states.
     * - nts_pre.m: # of edges.
     * - nts_pre.graph: an array of reverse graph nodes
     * - nts_pre.outdeg: an out degree table of the graph nts_pre
     * - nts_pre.in_edge: an array of in-going edges
     * - action: # of actions.
     */
    template<typename S>
    void BSolver::load_abstraction(abstraction<S> &abst) {
	NODE_PRE *nts_pre;
	int a,*outdeg;
	long long x1,x2,x,n0,m0;
	EDGE *p5;
	
	fts nts = abst._ts;
	_sol.nts_pre.n = n0 = nts._nx;
	_sol.nts_pre.m = m0 = nts._ntrans;
	_sol.action = nts._nu;
	_sol.nts_pre.graph=nts_pre=(NODE_PRE*)calloc(n0, sizeof(NODE_PRE));
	_sol.nts_pre.outdeg=outdeg=(int*)calloc(n0 * _sol.action, sizeof(int));
	p5=_sol.nts_pre.in_edge=(EDGE*)malloc((m0 + n0)*sizeof(EDGE));

	x=-1;
	for (long long i = 0; i < n0; ++i) {
	    for (long long j = 0; j < _sol.action; ++j) {
		for (const auto & item : nts.get_pre(i, j)) {
		    x1 = item; a = j; x2 = i;
		    (*(outdeg + x1*_sol.action + a))++;
		    if(x != x2) {
			if(x!=-1)
			    p5++->next=0;
			(nts_pre+x2)->in_edge=p5;
			x = x2;
		    }
		    p5++->next=1; p5->x=x1; p5->a=a;
		}
	    }
	}
	p5->next=0;
    }


    /**
     * 1) pre_saConstruct a backward graph for NTS, and do safety processing.
     * 2) Convert the backward graph to a forward graph.
     * 4) Take product of the forward graph and the DBA.
     * 5) Solve the buchi game on the product.
     */
    template<typename S>
    void BSolver::generate_product(abstraction<S> &abst) {
	int *labels = new int[abst._labels.size()];
	for(size_t i = 0; i < abst._labels.size(); ++i)
	    labels[i] = abst._labels[i];

	/**
	 * Construct and assign memory to
	 * 1)head->nts_pre, 2)encode1[n0], 3)decode1[n0_1].
	 */
	safety_pre(&_sol);

	/**
	 * Construct and assign memory to
	 * head->nts_post by converting from head->nts_pre. 
	 * 
	 * Release memory of head->nts_pre.
	 */
	pre2post(&_sol, labels);

	/**
	 * Construct and assign memory to
	 * 1)head->nts_product, 2)encode2[n0_1*n1], 3)decode2[np]
	 * 
	 * Release memory of 
	 * 1)encode2, 2)dba->acc
	 */
	post2product(&_sol);
	
	delete[] labels;
    }


} // namespace rocs

#endif
