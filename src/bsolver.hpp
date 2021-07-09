/**
 *  A class for solving a Buchi game on product system of an %NTS and a %DBA.
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
#include "buchi.h"


namespace rocs {
  
    /**
     * \brief The abstraction-based engine solver 
     * 
     * It is a wrapper class of the core data structures and algorithms defined in C files buchi.h and buchi.c
     *
     */
    class BSolver
    {
    public:
	HEAD _sol; /**< A HEAD of solver info */

	/**
	 * \brief The default constructor.
	 *
	 * No other constructors are allowed.
	 */
	BSolver() {initialization(&_sol);}
	BSolver(const BSolver&) = delete;
	BSolver(BSolver&&) = delete;
	BSolver& operator=(const BSolver&) = delete;
	BSolver& operator=(BSolver&&) = delete;

	/**
	 * \brief A destructor to release memory.
	 */
	~BSolver() {free_memory(&_sol, 1);}

	/**
	 * \brief Construct the DBA struct in _sol.
	 *
	 * Assign values to
	 * - n: # of %DBA states.
	 * - k: # of atomic propositions.
	 * - ini: the initial %DBA state.
	 * - acc: an array (k elements) of BOOL denoting accepting or not.
	 * - q_prime: the lookup table of size (2^k x n), recording transition relation.
	 * - q_prime_size: the total # of elements in q_prime (2^k x n).
	 *
	 * @param[in] nAP # of atomic propositions
	 * @param[in] nNodes # of %DBA states
	 * @param[in] q0 The initial state (UintSmall type)
	 * @param[in] acc The accepting states (a vector of UintSmall variables)
	 * @param[in] arrayM arrayM The transition matrix (a 2d array of UintSmall variables)
	 */
	void construct_dba(int nAP, int nNodes, int q0,
			   std::vector<rocs::UintSmall> &acc,
			   std::vector<std::vector<rocs::UintSmall>> arrayM);

	/**
	 * \brief Construct a graph with in-going edges in _sol.
	 *
	 * - nts_pre.n: # of graph states.
	 * - nts_pre.m: # of edges.
	 * - nts_pre.graph: an array of reverse graph nodes
	 * - nts_pre.outdeg: an out degree table of the graph nts_pre
	 * - nts_pre.in_edge: an array of in-going edges
	 * - action: # of actions.
	 *
	 * @param[in] abst An abstraction object
	 */
	template<typename S>
	void load_abstraction(abstraction<S> &abst);

	/**
	 * \brief Take product of the forward graph and the %DBA.
	 *
	 * - Construct a backward graph for %NTS, and do safety processing.
	 * - Convert the backward graph to a forward graph.
	 * - Take product of the forward graph and the %DBA.
	 *
	 * @param[in] abst An abstraction object
	 */
	template<typename S>
	void generate_product(abstraction<S> &abst);

	/**
	 * \brief Solve the buchi game on the product.
	 */
	void solve_buchigame_on_product() {
	    buchi_and_controller(&_sol);
	}

	/**
	 * \brief Write controller to a file
	 * 
	 * @param[in] filename The name of the output file
	 */
	void write_controller_to_txt(char *filename) {
	    write_controller(&_sol, filename);
	}
	
    };

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
    }//end load_abstraction
    
    template<typename S>
    void BSolver::generate_product(abstraction<S> &abst) {
	int *labels = new int[abst._labels.size()];
	for(size_t i = 0; i < abst._labels.size(); ++i)
	    labels[i] = abst._labels[i];

	/*
	 * Construct and assign memory to
	 * 1)head->nts_pre, 2)encode1[n0], 3)decode1[n0_1].
	 */
	safety_pre(&_sol);

	/*
	 * Construct and assign memory to
	 * head->nts_post by converting from head->nts_pre. 
	 * 
	 * Release memory of head->nts_pre.
	 */
	pre2post(&_sol, labels);

	/*
	 * Construct and assign memory to
	 * 1)head->nts_product, 2)encode2[n0_1*n1], 3)decode2[np]
	 * 
	 * Release memory of 
	 * 1)encode2, 2)dba->acc
	 */
	post2product(&_sol);
	
	delete[] labels;
    }//end generate_product


} // namespace rocs

#endif
