/**
 *  hdf5io.h
 *
 *  Input/output classes to .h5 files.
 *
 *  Created by Yinan Li on Aug 16, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _hdf5io_h
#define _hdf5io_h

#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <vector>

#include <H5Cpp.h>

#include "grid.h"
#include "csolver.h"
#include "transition.hpp"


namespace rocs {

    class h5FileHandler {
    private:
	const H5std_string _filename;
	const H5::H5File _h5file;

    public:
	h5FileHandler() = delete;
	h5FileHandler(std::string f, unsigned int flag) :
	    _filename(f), _h5file(_filename, flag) {}

	/* Write data */
	int write_real_number(const double x, const std::string varname);
	int write_real_array(const std::vector<double> &x,
			     const std::string varname);
	int write_2d_real_array(const std::vector< std::vector<double> > &x,
				const std::string varname);
	int write_uint_number(const size_t d, const std::string varname);
	int write_uint_array(const std::vector<size_t> &arr,
			     const std::string varname);
	int write_state_space(const ivec &ws, const std::string varname);
	int write_input_values(const grid &ugrid, const std::string varname);
	int write_ivec_array(const std::vector<ivec> &arr, const std::string varname);

	int write_uniform_grids(const grid &g, const std::string varname);
	int write_transitions(const fts &trans,
			      const std::string vtranspost,
			      const std::string vpost,
			      const std::string vptrpost,
			      const std::string vtranspre,
			      const std::string vpre,
			      const std::string vptrpre);


	int write_sptree_leaves(const SPtree &ctlr, SPnode *ptrn);
	int write_sptree_controller(const CSolver &sol) {
	    write_sptree_leaves(sol._ctlr, sol._ctlr._root);
	    return 0;
	}

	int write_transitions(const fts &trans);

	template<typename S>
	void write_problem_setting(const S &sys, const CSolver &sol) {
	    write_real_number(sys._tau, "ts");
	    write_state_space(sys._workspace, "X");
	    write_input_values(sys._ugrid, "U");
	    write_ivec_array(sol._goal, "G");
	    write_ivec_array(sol._obs, "xobs");
	}


	/* Read data */
	int read_uint_number(size_t *ptrd, const std::string varname);
	int read_uint_array(std::vector<size_t> &arr,
			    const std::string varname);
	int read_real_array(std::vector<double> &arr,
			    const std::string varname);
	int read_transitions(fts &trans);

    };//h5FileHandler class


    /**
     * Write control problem setup and synthesized controller into a .h5 file.
     * @param specfile the specification file name.
     */
    template<typename S>
    void write_csolvers_to_h5(const S &sys, std::string specfile,
			      std::vector<rocs::CSolver*> &w) {
	std::string datafile;
	std::vector<std::string> tokens;
	boost::split(tokens, specfile, boost::is_any_of("."));
	for (size_t i = 0; i < w.size(); ++i) {
	    // w[i]->_timer = t[i];
	    // std::cout << std::endl << "S-domain for q" << std::to_string(i) << '\n';
	    // w[i]->print_controller_info();

	    datafile = "data_" + tokens[0] + "_w" + std::to_string(i) + ".h5";
	    rocs::h5FileHandler h5f(datafile, H5F_ACC_TRUNC);
	    h5f.write_problem_setting(sys, *(w[i]));
	    h5f.write_sptree_controller(*(w[i]));
	}
    }

}//namespace rocs


#endif
