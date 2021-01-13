/**
 *  testWinGraph.cpp
 *
 *  Test the winning graph generation in Patcher.
 *
 *  Created by Yinan Li on Nov. 9, 2020.
 *  Hybrid Systems Group, University of Waterloo.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <algorithm>
#include <cmath>
#include <sys/stat.h>

#include "src/grid.h"
#include "src/definitions.h"

#include "src/abstraction.hpp"

#include "src/bsolver.hpp"
#include "src/patcher.h"

#include "src/hdf5io.h"


int main()
{
    clock_t tb, te;
    
    /* Launch patcher */
    rocs::Patcher local;
    std::string graphfile = "graph_winning.h5";
    struct stat buffer;
    if(stat(graphfile.c_str(), &buffer) == 0) {
    	/* Read from a file */
    	std::cout << "\nReading winning graph...\n";
    	rocs::h5FileHandler graphRdr(graphfile, H5F_ACC_RDONLY);
    	tb = clock();
    	graphRdr.read_winning_graph(local);
    	te = clock();
    	std::cout << "Time of reading graph: " << (float)(te - tb)/CLOCKS_PER_SEC << '\n';
    } else {
    	std::cout << "Graph file doesn't exist.\n";
    	return 1;
    }


    /* Test post-pre consistency */
    size_t na = local._na;
    size_t si, sk, post, pre, k;
    bool suc = 0;
    for(size_t i = 0; i < local._nwin; ++i) {
	for(size_t j = 0; j < na; ++j) {
	    si = local._winfts._ptrpost[i*na+j];
	    for(size_t p=si; p<si+local._winfts._npost[i*na+j]; ++p) {
		post = local._winfts._idpost[p];
		k = local._idmap[post];
		/* Test if the pre of post by j contains i */
		suc = 0;
		sk = local._winfts._ptrpre[k*na+j];
		for(size_t pp=local._winfts._ptrpre[k*na+j];
		    pp<sk+local._winfts._npre[k*na+j]; ++pp) {
		    pre = local._winfts._idpre[pp];
		    if(local._idmap[pre] == i) {
			suc = 1;
			break;
		    }
		}
		if(!suc) {//two cases: npre(k,j)=0 or no i in npre(k, j)
		    std::cout << "Post and pre transitions are inconsistent "
			      << i << "->(" << j << ")->" << k << '\n';
		    return 1;
		}
		    
	    }
	}
    }
    std::cout << "Every post transition has its corresponding pre transition.\n";

    for(size_t i = 0; i < local._nwin; ++i) {
	for(size_t j = 0; j < na; ++j) {
	    si = local._winfts._ptrpre[i*na+j];
	    for(size_t p=si; p<si+local._winfts._npre[i*na+j]; ++p) {
		pre = local._winfts._idpre[p];
		k = local._idmap[pre];
		/* Test if the post of pre by j contains i */
		suc = 0;
		sk = local._winfts._ptrpost[k*na+j];
		for(size_t pp=local._winfts._ptrpost[k*na+j];
		    pp<sk+local._winfts._npost[k*na+j]; ++pp) {
		    post = local._winfts._idpost[pp];
		    if(local._idmap[post] == i) {
			suc = 1;
			break;
		    }
		}
		if(!suc) {//two cases: npost(k,j)=0 or no i in npost(k, j)
		    std::cout << "Post and pre transitions are inconsistent "
			      << k << "->(" << j << ")->" << i << '\n';
		    return 1;
		}
		    
	    }
	}
    }
    std::cout << "Every pre transition has its corresponding post transition.\n";

    return 0;
}
