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
#include <boost/dynamic_bitset.hpp>
#include <vector>

#include <H5Cpp.h>

#include "grid.h"
#include "csolver.h"
#include "patcher.h"
#include "dsolver.h"
// #include "bsolver.hpp" //this will give a linking error
#include "buchi.h"

#include "transition.hpp"


namespace rocs {
    
    template<typename T>
    H5::PredType get_datatype();
    

    class h5FileHandler {
    private:
	const H5std_string _filename;
	const H5::H5File _h5file;

    public:
	h5FileHandler() = delete;
	/**
	 * Constructor. flag=
	 * for read: H5F_ACC_RDONLY, for write: H5F_ACC_TRUNC
	 */
	h5FileHandler(std::string f, unsigned int flag) :
	    _filename(f), _h5file(_filename, flag) {}

	/* Write data */
	void create_group(const std::string groupname) {
	    H5::Group group(_h5file.createGroup(groupname));
	}

	/**
	 * Write/Read a number/array into/from a .h5 file:
	 * can be double, int, long long, unsigned int, etc.
	 * @param d/arr the integer number/array
	 * @param varname the variable name to be written/read.
	 */
	template<typename T>
	int write_number(const T x, const std::string varname);
	
	template<typename T>
	int write_array(const std::vector<T> &arr, const std::string varname) {
	    if(arr.empty()) {
		std::cout << "hdf5FileHandler::write_array: Input array is empty. Writing is abandoned.\n ";
		return 1;
	    } else {
		return write_array<T>(&(arr[0]), arr.size(), varname);
	    }
	}
	template<typename T>
	int write_array(const T *arr, const size_t len, const std::string varname);
	
	template<typename T>
	int write_2d_array(const std::vector< std::vector<T> > &arr,
			   const std::string varname);
	template<typename T>
	int write_2d_array(const std::vector<T> &arr, const size_t *len,
			  const std::string varname);
	template<typename T>
	int write_2d_array(const T *arr, const size_t *len,
			   const std::string varname);

	int write_ivec_array(const std::vector<ivec> &arr, const std::string varname);

	/**
	 * Write system settings
	 */
	int write_state_space(const ivec &ws, const std::string varname);
	int write_input_values(const grid &ugrid, const std::string varname);
	template<typename S>
	void write_problem_setting(const S &sys) {//, const CSolver &sol) {
	    write_number<double>(sys._tau, "ts");
	    write_state_space(sys._workspace, "X");
	    write_input_values(sys._ugrid, "U");
	    // write_ivec_array(sol._goal, "G");
	    // write_ivec_array(sol._obs, "xobs");
	}
	
	/**
	 * Write DSolver into a .h5 file.
	 * @param filename the name of the file to be saved.
	 */
	int write_discrete_controller(HEAD *sol);
	int write_discrete_controller(const DSolver &dsol);
	int write_transitions(const fts &trans);
	int write_winning_graph(const Patcher &patcher);

	/**
	 * Write CSolver into a .h5 file.
	 * @param filename the name of the file to be saved.
	 */
	int write_sptree_leaves(const SPtree &ctlr, SPnode *ptrn);
	int write_sptree_controller(const CSolver &sol) {
	    write_sptree_leaves(sol._ctlr, sol._ctlr._root);
	    return 0;
	}


	/* Read data */
	template<typename T>
	int read_number(T *ptrd, const std::string varname);
	template<typename T>
	int read_array(std::vector<T> &arr, const std::string varname);
	template<typename T>
	int read_2d_array(std::vector<T> &arr, size_t *len,
			  const std::string varname);
	
	int read_transitions(fts &trans);
	int read_winning_graph(Patcher &patcher);
	int read_discrete_controller(boost::dynamic_bitset<> &win,
				     boost::dynamic_bitset<> &lsctlr,
				     size_t *cdims,
				     std::vector<size_t> &optctlr,
				     std::vector<double> &value);
	int read_discrete_controller(std::vector<long long> &w_x0,
				     std::vector<long long> &encode3,
				     std::vector<NODE_POST> &nts_ctrlr,
				     std::vector<CTRL> &ctrl,
				     std::vector<int> &q_prime);
	int read_sptree_controller(std::vector<double> &pavings, size_t *pdims,
				   std::vector<int> &tag,
				   boost::dynamic_bitset<> &cntl, size_t *cdims);

    };//h5FileHandler class


    template<typename T>
    int h5FileHandler::write_number(const T d, const std::string varname) {
	hsize_t dim[1] = {1};
	H5::DataSpace dataspace(1, dim);
	// H5::PredType h5dt = get_datatype<T>();
	H5::DataType datatype(get_datatype<T>());
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	dataset.write(&d, datatype);
	return 0;
    }

    
    template<typename T>
    int h5FileHandler::write_array(const T *arr, const size_t len,
				   const std::string varname) {
	if(!len) {
	    std::cout << "hdf5FileHandler::write_array: "
		      << varname << " is empty."
		      << " Writing is abandoned.\n ";
	    return 1;
	}
	hsize_t dim[1]= {len};
	H5::DataSpace dataspace(1, dim);

	H5::PredType h5dt = get_datatype<T>();
	H5::DataType datatype(h5dt);
	H5::DataSet dataset;

	if(len > 1000) {
	    hsize_t chunkSize = dim[0] / 10;
	    hsize_t chunkdim[1] = {chunkSize};
	    H5::DSetCreatPropList plist;
	    plist.setChunk(1, chunkdim);
	    plist.setDeflate(6);
	    dataset = _h5file.createDataSet(varname, datatype, dataspace, plist);
	} else {
	    dataset = _h5file.createDataSet(varname, datatype, dataspace);
	}
	dataset.write(arr, datatype);
	return 0;
    }

    template<typename T>
    int h5FileHandler::write_2d_array(const std::vector< std::vector<T> > &arr,
				     const std::string varname) {
	if(arr.empty()) {
	    std::cout << "hdf5FileHandler::write_2d_array: "
		      << varname << " is empty."
		      << " Writing is abandoned.\n ";
	    return 1;
	}
	size_t dim[2];
	dim[0]= arr.size();
	dim[1] = arr[0].size();
	
	/* create the buffer data for writing */
	T *data = new T[dim[0]*dim[1]];
	for(hsize_t i = 0; i < dim[0]; ++i) {
	    for(hsize_t j = 0; j < dim[1]; ++j)
		*(data+dim[1]*i+j) = arr[i][j];
	}
	
	write_2d_array<T>(data, dim, varname);
	
	delete[] data;
	return 0;
    }
    template<typename T>
    int h5FileHandler::write_2d_array(const std::vector<T> &arr,
				      const size_t *len,
				      const std::string varname) {
	if(!len[0] || !len[1]) {
	    std::cout << "hdf5FileHandler::write_2d_array: "
		      << varname << " is empty."
		      << " Writing is abandoned.\n ";
	    return 1;
	}
	size_t dim[2];
	dim[0]= len[0];
	dim[1] = len[1];
	
	/* create the buffer data for writing */
	T *data = new T[dim[0]*dim[1]];
	for(hsize_t i = 0; i < dim[0]*dim[1]; ++i) {
		data[i] = arr[i];
	}
	
	write_2d_array<T>(data, dim, varname);
	
	delete[] data;
	return 0;
    }
    template<typename T>
    int h5FileHandler::write_2d_array(const T *arr, const size_t *len,
				     const std::string varname) {
	if(!len[0] || !len[1]) {
	    std::cout << "hdf5FileHandler::write_2d_array: "
		      << varname << " is empty."
		      << " Writing is abandoned.\n ";
	    return 1;
	}
	hsize_t dim[2];
	dim[0] = len[0];
	dim[1] = len[1];
	H5::DataSpace dataspace(2, dim);
	H5::DataType datatype(get_datatype<T>());
	H5::DataSet dataset;
	
	/* write to h5 dataset */
	if(dim[0] > 1000) {
	    hsize_t chunkSize = dim[0] / 10;
	    hsize_t chunkdim[2] = {chunkSize, dim[1]};
	    H5::DSetCreatPropList plist;
	    plist.setChunk(2, chunkdim);
	    plist.setDeflate(6);
	    dataset = _h5file.createDataSet(varname, datatype, dataspace, plist);
	} else {
	    dataset = _h5file.createDataSet(varname, datatype, dataspace);
	}
	dataset.write(arr, datatype);
	return 0;
    }

    template<typename T>
    int h5FileHandler::read_number(T *ptrd, const std::string varname) {
	T data[1];
	H5::DataSet dataset = _h5file.openDataSet(varname);
	H5::DataSpace dataspace = dataset.getSpace();
	// int rank_dnx = dnxspace.getSimpleExtentNdims();
	hsize_t dim[1] = {1};
	H5::DataSpace memspace(1, dim);
	dataset.read(data, get_datatype<T>(), memspace, dataspace);
	*ptrd = data[0];
	return 0;
    }

    template<typename T>
    int h5FileHandler::read_array(std::vector<T> &arr,
				  const std::string varname) {
	H5::DataSet dataset = _h5file.openDataSet(varname);
	H5::DataSpace dataspace = dataset.getSpace();
	int rank = dataspace.getSimpleExtentNdims();
	hsize_t dim[1];
	rank = dataspace.getSimpleExtentDims(dim);
	H5::DataSpace memspace(1, dim);

	T *data = new T[dim[0]];
	dataset.read(data, get_datatype<T>(), memspace, dataspace);
	arr.resize(dim[0]);
	for(hsize_t i = 0; i < dim[0]; ++i)
	    arr[i] = data[i];
	delete[] data;
	return 0;
    }

    template<typename T>
    int h5FileHandler::read_2d_array(std::vector<T> &arr, size_t *dims,
				     const std::string varname) {
	H5::DataSet dataset = _h5file.openDataSet(varname);
	H5::DataSpace dataspace = dataset.getSpace();
	// int rank = dataspace.getSimpleExtentNdims();
	hsize_t hdim[2];
	int ndims = dataspace.getSimpleExtentDims(hdim, NULL);
	H5::DataSpace memspace(2, hdim);
	dims[0] = (size_t)hdim[0];
	dims[1] = (size_t)hdim[1];
	
	T *data = new T[hdim[0]*hdim[1]];
	dataset.read(data, get_datatype<T>(), memspace, dataspace);
	arr.resize(hdim[0]*hdim[1]);
	arr.shrink_to_fit();
	for(hsize_t i = 0; i < hdim[0]*hdim[1]; ++i)
	    arr[i] = data[i];
	delete[] data;
	return 0;
    }


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

	    datafile = "controller_" + tokens[0] + "_w" + std::to_string(i) + ".h5";
	    rocs::h5FileHandler h5f(datafile, H5F_ACC_TRUNC);
	    h5f.write_problem_setting(sys);
	    h5f.write_ivec_array(w[i]->_goal, "G");
	    h5f.write_ivec_array(w[i]->_obs, "xobs");
	    h5f.write_sptree_controller(*(w[i]));
	}
    }
    

}//namespace rocs


#endif
