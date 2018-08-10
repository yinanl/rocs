/**
 *  testDSolver.cpp
 *
 *  Test DSolver for reachability and invariance games.
 *
 *  Created by Yinan Li on Sept. 05, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TESTDSOLVER
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>

#include <cmath>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include "src/dsolver.h"
#include "src/matlabio.h"



/* ----------------- test predicate ----------------- */
template <class T>
boost::test_tools::predicate_result
my_check_permutation(std::vector<T> act, std::vector<T> exp) {
    
    boost::test_tools::predicate_result res(false);
    res.message() << "Elements [exp:";
    for (typename std::vector<T>::iterator iter = exp.begin();
         iter != exp.end(); ++iter) {
        res.message() << *iter << " ";
    }
    res.message() << "act:";
    for (typename std::vector<T>::iterator iter = act.begin();
         iter != act.end(); ++iter) {
        res.message() << *iter << " ";
    }
    res.message() << "].";
    
    
    if (act.size() != exp.size()) {
        res.message() << "Different size.";
        return res;
        
    } else {
        
        bool permute = std::is_permutation(act.begin(), act.end(), exp.begin());
        if (!permute) {
            res.message() << "Different elements.";
            return res;
            
        } else {
            return true;
        }
    }
    
}


/* ----------------- test fixture ----------------- */
struct abstFixture{
    
    /* member variables */
    std::string spec;
    int nU;
    int nV;
    rocs::fts transitions;
    std::vector<size_t> goalSet, avoidSet, safeSet;
    std::vector<size_t> expReachSet, expInvSet;

    /* constructor */
    abstFixture(std::string gc):spec(gc), nU(0), nV(0) {
        BOOST_TEST_MESSAGE("Test setup...");
        
        /* read a specification */
        std::ifstream specfile;
        specfile.open(spec);
        if (specfile.is_open()) {
            
            /* start parsing file: read per line */
            std::string fline, name;
            std::string delims("=,");
	    
	    /* get number of states, inputs and target/avoid sets */
            while (std::getline(specfile, fline)) {
                std::vector<std::string> tokens;
		boost::split(tokens, fline, boost::is_any_of(delims));
		for( std::vector<std::string>::iterator iter = tokens.begin()+1;
		     iter != tokens.end(); ++iter) {
                                
		    if (tokens[0] == "nU") {
			this->nU = std::stoi(*iter);
		    } else if (tokens[0] == "nV") {
			this->nV = std::stoi(*iter);
		    } else if (tokens[0] == "goals") {
			this->goalSet.push_back( std::stoi(*iter) );
		    } else if (tokens[0] == "avoids") {
			this->avoidSet.push_back( std::stoi(*iter) );
		    } else if (tokens[0] == "exp_reach") {
			this->expReachSet.push_back(std::stoi(*iter));
		    } else if (tokens[0] == "exp_inv") {
			this->expInvSet.push_back(std::stoi(*iter));
		    } else {
			continue;
		    }
		} // end token loop	
               
	    } // end while

	    /* get transitions */
	    transitions.init(nV, nU);
	    /* assign npost, ptrpost and idpost */
	    for (int i = 0; i < nV; ++i) {
		for (int j = 0; j < nU; ++j) {
		    bool first = 1;
		    specfile.clear();
		    specfile.seekg(0, std::ios::beg);
		    while (std::getline(specfile, fline)) {
			std::vector<std::string> tokens;
			boost::split(tokens, fline, boost::is_any_of(delims));
			if (tokens[0]=="trans") {
			    if (first) {
				transitions._ptrpost[i*nU + j] = transitions._ntrans;
				first = 0;
			    }
			    if(std::stoi(tokens[1]) == i && std::stoi(tokens[3]) == j) {
				// std::cout << fline << '\n';
				transitions._npost[i*nU + j] ++;
				transitions._ntrans ++;
				transitions._idpost.push_back(std::stoi(tokens[2]));
				transitions._cost.push_back(std::stof(tokens[5]));

				transitions._npre[std::stoi(tokens[2])*nU + j] ++;
			    }
			} // end assignment of one transition from i under j
		    }  // end assignment of all transitions from i under j
		} // end input loop
	    } // end state loop

	    // /********** logging **********/
	    // std::cout << "Number of transitions: " << transitions._ntrans << '\n';
	    // std::cout << "Number of states: " << transitions._nx << '\n';
	    // std::cout << "Number of inputs: " << transitions._nu << '\n';
	    // for (size_t i = 0; i < nV; ++i) {
	    // 	for (size_t j = 0; j < nU; ++j) {
	    // 	    for (int k = 0; k < transitions._npost[i*nU+j]; ++k) {
	    // 		std::cout << i << ','
	    // 			  << transitions._idpost[transitions._ptrpost[i*nU+j]+k] << ','
	    // 			  << j << ','
	    // 			  << transitions._cost[transitions._ptrpost[i*nU+j]+k]
	    // 			  << '\n';
	    // 	    }
	    // 	}
	    // }
	    // /********** logging **********/

	    /* assign npre, ptrpre, and idpre */
	    transitions._idpre.resize(transitions._ntrans);  // initialize the size of pre's

	    /* assign _ptrpre */
	    size_t sum = 0;
	    for (size_t i = 0; i < nV; ++i) {
		for (size_t j = 0; j < nU; ++j) {
		    transitions._ptrpre[i*nU + j] = sum;
		    sum += transitions._npre[i*nU + j];
		}
	    }
	    /* assign _idpre */
	    std::vector<size_t> precount(nU*nV, 0);
	    size_t idtspre;
	    for (size_t i = 0; i < nV; ++i) {
		for (size_t j = 0; j < nU; ++j) {
		    for (int ip = 0; ip < transitions._npost[i*nU + j]; ++ip) {
		
			idtspre = transitions._idpost[transitions._ptrpost[i*nU+j] + ip]*nU + j;
			transitions._idpre[transitions._ptrpre[idtspre] + precount[idtspre]] = i;
			precount[idtspre] ++;
		    }
		}
	    }

	    /********** logging **********/
	    for (size_t i = 0; i < nV; ++i) {
	    	for (size_t j = 0; j < nU; ++j) {
	    	    for (int k = 0; k < transitions._npre[i*nU+j]; ++k) {
	    		std::cout << i << ','
	    			  << transitions._idpre[transitions._ptrpre[i*nU+j]+k] << ','
	    			  << j << '\n';
	    	    }
	    	}
	    }
	    /********** logging **********/
            
            specfile.close();
        } // end if
        
    }  // abstFixture
    
    
    // destructor
    ~abstFixture(){
        BOOST_TEST_MESSAGE("Test teardown...");
    }
    
};

struct GF_ReachInv:abstFixture {
    GF_ReachInv():abstFixture("./test_graphs/graph_1"){}
};

struct GF_mReach:abstFixture {
    GF_mReach():abstFixture("./test_graphs/graph_2"){}
};

struct GF_cReach:abstFixture {
    GF_cReach():abstFixture("./test_graphs/graph_3"){}
};

struct GF_cInv:abstFixture {
    GF_cInv():abstFixture("./test_graphs/graph_4"){}
};




/*
  test case 1: graph_3, controlled reach set
*/
BOOST_FIXTURE_TEST_CASE(test_controlled_reach, GF_cReach)
{
    rocs::DSolver solver(&transitions);
    solver.reachability(goalSet, avoidSet);
    
    std::vector<size_t> winset;
    std::cout << "Winning set: ";
    for (size_t i = 0; i < nV; ++i) {
	if (solver._winset[i]) {
	    winset.push_back(i);
	    std::cout << i << ", ";
	}
    }
    std::cout << '\n';
    
    BOOST_CHECK(my_check_permutation(winset, expReachSet));  // check winset
    rocs::matWriter wtr("data_testgraph3.mat");
    wtr.open();
    wtr.write_discrete_controller(solver, "leastctlr", "optctlr");
    wtr.close();
}


/*
  test case 2: graph_4, controlled invariant set
*/
BOOST_FIXTURE_TEST_CASE(test_controlled_invariant, GF_cInv)
{
    rocs::DSolver solver(&transitions);
    solver.invariance(goalSet);
    std::vector<size_t> winset;
    for (size_t i = 0; i < nV; ++i) {
	if (solver._winset[i])
	    winset.push_back(i);
    }
    
    BOOST_CHECK(my_check_permutation(winset, expInvSet));
}


/*
  test case 3: graph_1, controlled reach avoid set
*/
BOOST_FIXTURE_TEST_CASE(test_controlled_reach_avoid_stay, GF_ReachInv)
{
    BOOST_CHECK(true);  // check winset
    
    
}
