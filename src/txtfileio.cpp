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

	    if (c.isleaf(current)) { //write leaves
		lower = current->_box.getinf();
		upper = current->_box.getsup();
		for (int i = 0; i < node->_box.getdim(); ++ i)
		    _txtfile << lower[i] << ' ' << upper[i] << ' ';

		_txtfile << " " << current->_tag << ' ';

		for (int u = 0; u < node->_cntl.size(); ++ u)
		    _txtfile << ' ' << current->_cntl[u];

		_txtfile << '\n';
	    }

	    if (current->_left)
		stk.push(current->_left);

	    if (current->_right)
		stk.push(current->_right);
	    
	}
    }


}//namespace rocs
