/**
 *  txtfileio.cpp
 *
 *  Definition of a text file input/output class.
 *
 *  Created by Yinan Li on August 09, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include "txtfileio.h"

namespace rocs {
    
    void txtWriter::write_uniform_grids(const grid &g) {
	if (!_txtfile.is_open()) {open();}
	if (g._data.empty()) {
	    std::cout << "txtWriter: No valid data in the grid.\n";
	    return;
	}
	for (size_t row = 0; row < g._nv; ++row) {
	    _txtfile << row << " ";
	    for (int col = 0; col < g._dim; ++col) {
		_txtfile << g._data[row][col] << ' ';
	    }
	    _txtfile << '\n';
	}
    }//txtWriter::write_uniform_grids

    void txtWriter::write_sptree_leaves(const SPtree &c, SPnode *node) {

	if (node == NULL) {
	    std::cout << "txtWriter: Nothing Generated." << std::endl;
	    return;
	}

	/* write data to the file line by line */
	std::stack<SPnode*> stk;
	SPnode *current = node;
	stk.push(current);

	std::vector<double> lower, upper;
	while (!stk.empty()) {
	    current = stk.top();
	    stk.pop();
	    // std::cout<< j << ", " << current->_box << ", "; //logging

	    if (c.isleaf(current)) { //write leaves
		lower = current->_box.getinf();
		upper = current->_box.getsup();
		// std::cout<< j << ", " << current->_box << std::endl; //for logging
		for (int i = 0; i < node->_box.getdim(); ++ i)
		    _txtfile << lower[i] << ' ' << upper[i] << ' ';

		_txtfile << " " << current->_tag << ' ';

		for (size_t u = 0; u < node->_cntl.size(); ++ u)
		    _txtfile << ' ' << current->_cntl[u];

		_txtfile << '\n';
	    }

	    if (current->_left)
		stk.push(current->_left);

	    if (current->_right)
		stk.push(current->_right);

	    // std::cout << std::endl; //logging
	}
    }//txtWriter::write_sptree_leaves

    void txtWriter::write_post_transitions(const fts &nts, const grid &g,
					   const double xlb[], const double xub[]) {
	if (!_txtfile.is_open()) {open();}

	size_t n = nts._nx;
	size_t m = nts._nu;
	size_t i = 0;
	bool inx = true;
	for (size_t row = 0; row < n; ++row) {
	    inx = true;
	    for (int k = 0; k < g._dim; ++k) {
		if (g._data[row][k] > xub[k] || g._data[row][k] < xlb[k]) {
		    inx = false;
		    break;
		}
	    }
	    if (inx) {
		// std::cout << row << ':';
		// for (int k = 0; k < g._dim; ++k)
		//     std::cout << g._data[row][k] << ',';
		// std::cout << '\n';
		for (size_t col = 0; col < m; ++col) {
		    // std::cout << nts._npost[row*m+col] << '\n';
		    for (int l = 0; l < nts._npost[row*m+col]; ++l) {
			_txtfile << ++i << ' ';
			_txtfile << row << ' ' << col << ' ' << nts._idpost[nts._ptrpost[row*m+col]+l];
			_txtfile << '\n';
		    }
		}
	    }
	}
    }

    void txtWriter::write_pre_transitions(const fts &nts, const grid &g,
					  const double xlb[], const double xub[]) {
	if (!_txtfile.is_open()) {open();}

	size_t n = nts._nx;
	size_t m = nts._nu;
	size_t i = 0;
	bool inx = true;
	for (size_t row = 0; row < n; ++row) {
	    inx = true;
	    for (int k = 0; k < g._dim; ++k) {
		if (g._data[row][k] > xub[k] || g._data[row][k] < xlb[k]) {
		    inx = false;
		    break;
		}
	    }
	    if (inx) {
		std::cout << row << ':';
		for (int k = 0; k < g._dim; ++k)
		    std::cout << g._data[row][k] << ',';
		std::cout << '\n';
		for (size_t col = 0; col < m; ++col) {
		    // std::cout << nts._npost[row*m+col] << '\n';
		    for (int l = 0; l < nts._npre[row*m+col]; ++l) {
			_txtfile << ++i << ' ';
			_txtfile << nts._idpre[nts._ptrpre[row*m+col]+l] << ' ' << col << ' ' << row;
			_txtfile << '\n';
		    }
		}
	    }
	}
    }


}//namespace rocs
