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
    
    void txtFileHandler::write_uniform_grids(const grid &g) {
	if (!_txtfile.is_open()) {open(std::ios::out);}
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
    }//txtFileHandler::write_uniform_grids

    void txtFileHandler::write_sptree_leaves(const SPtree &c, SPnode *node) {

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
    }//txtFileHandler::write_sptree_leaves

    void txtFileHandler::write_post_transitions(const fts &nts, const grid &g,
					   const double xlb[], const double xub[]) {
	if (!_txtfile.is_open()) {open(std::ios::out);}

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

    void txtFileHandler::write_pre_transitions(const fts &nts, const grid &g,
					  const double xlb[], const double xub[]) {
	if (!_txtfile.is_open()) {open(std::ios::out);}

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

    void txtFileHandler:: read_discrete_controller(std::vector<long long> &w_x0,
						   std::vector<long long> &encode3,
						   std::vector<NODE_POST> &nts_ctrlr,
						   std::vector<CTRL> &ctrl,
						   std::vector<int> &q_prime) {
	if (!_txtfile.is_open())
	    open(std::ifstream::in | std::ifstream::binary);

	/* Read the winning states to w_x0 */
	long long w0;
	_txtfile >> w0; //read the # of winning set in X
	std::cout << "W_X0 size: " << w0 << '\n';
	w_x0.resize(w0);
	for(long long i = 0; i < w0; ++i) {
	    _txtfile >> w_x0[i];
	}
	std::cout << w_x0[0] << ',' << w_x0[1] << '\n';
	/* Read encode3:n0->n0_2 for getting the index in nts_ctlr */
	long long n0;
	_txtfile >> n0;
	std::cout << "encode3 size: " << n0 << '\n';
	encode3.resize(n0);
	for(size_t i = 0; i < encode3.size(); ++i) {
	    _txtfile >> encode3[i];
	}
	std::cout << encode3[0] << ',' << encode3[1] << '\n';
	/* Read nts_ctrlr */
	long long n0_2;
	_txtfile >> n0_2;
	std::cout << "nts_ctrlr size: " << n0_2 << '\n';
	nts_ctrlr.resize(n0_2);
	for(size_t i = 0; i < nts_ctrlr.size(); ++i) {
	    _txtfile >> nts_ctrlr[i].num_a;
	    _txtfile >> nts_ctrlr[i].label;
	    _txtfile >> nts_ctrlr[i].pos;
	}
	for(long long i = 0; i < 3; ++i)
	    std::cout << '(' << nts_ctrlr[i].num_a << ','
		      << nts_ctrlr[i].label << ','
		      << nts_ctrlr[i].pos << ") ";
	std::cout << '\n';
	/* Read ctrl: stores the control input u by searching for q. */
	long long wp;
	_txtfile >> wp;
	ctrl.resize(wp);
	std::cout << "Product controller size: " << wp << '\n';
	for(long long i = 0; i < wp; ++i) {
	    _txtfile >> ctrl[i].q;
	    _txtfile >> ctrl[i].u;
	}
	for(long long i = 0; i < 3; ++i)
	    std::cout << '(' << ctrl[i].q << ',' << ctrl[i].u << ") ";
	std::cout << '\n';
	/* Read q_prime: the DBA transition matrix */
	long long q_prime_size;
	_txtfile >> q_prime_size;
	std::cout << "q_prime size: " << q_prime_size << '\n';
	q_prime.resize(q_prime_size);
	for(size_t i = 0; i < q_prime.size(); ++i) {
	    _txtfile >> q_prime[i];
	}
    }


}//namespace rocs
