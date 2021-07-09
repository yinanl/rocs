/**
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
// #include <boost/multi_array.hpp>
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

    /**
     * \brief An hdf5 file input/output class.
     *
     * Read and write controller data into hdf5 data format.
     */
    class h5FileHandler {
    private:
	const H5std_string _filename; /**< The variable of the filename */
	const H5::H5File _h5file; /**< The H5File variable */

    public:
	h5FileHandler() = delete;
	/**
	 * \brief A constructor.
	 * 
	 * @param[in] f The file name
	 * @param[in] flag H5F_ACC_RDONLY(read), H5F_ACC_TRUNC(write)
	 */
	h5FileHandler(std::string f, unsigned int flag) :
	    _filename(f), _h5file(_filename, flag) {}

	/**
	 * \brief Create an hdf5 group.
	 *
	 * @param[in] The group name
	 */
	void create_group(const std::string groupname) {
	    H5::Group group(_h5file.createGroup(groupname));
	}

	/**
	 * \brief Write a number into an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in] x The number
	 * @param[in] varname The variable name to be written
	 */
	template<typename T>
	int write_number(const T x, const std::string varname);

	/**
	 * \brief Write an std::vector into an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in] arr The vector
	 * @param[in] varname The variable name to be written
	 */
	template<typename T>
	int write_array(const std::vector<T> &arr, const std::string varname) {
	    if(arr.empty()) {
		std::cout << "hdf5FileHandler::write_array: Input array is empty. Writing is abandoned.\n ";
		return 1;
	    } else {
		return write_array<T>(&(arr[0]), arr.size(), varname);
	    }
	}
	
	/**
	 * \brief Write an array into an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in] arr The array
	 * @param[in] varname The variable name to be written
	 */
	template<typename T>
	int write_array(const T *arr, const size_t len, const std::string varname);

	/**
	 * \brief Write a 2d std::vector into an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in] arr The vector
	 * @param[in] varname The variable name to be written
	 */
	template<typename T>
	int write_2d_array(const std::vector< std::vector<T> > &arr,
			   const std::string varname);

	/**
	 * \brief Write a 2d array into an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in] arr The array
	 * @param[in] varname The variable name to be written
	 */
	template<typename T>
	int write_2d_array(const std::vector<T> &arr, const size_t *len,
			  const std::string varname);

	/**
	 * \brief Write a 2d array into an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in] arr The array
	 * @param[in] varname The variable name to be written
	 */
	template<typename T>
	int write_2d_array(const T *arr, const size_t *len,
			   const std::string varname);

	/**
	 * \brief Write a std::vector of interval vectors into an hdf5 file.
	 *
	 * @param[in] arr The vector of interval vectors
	 * @param[in] varname The variable name to be written
	 */
	int write_ivec_array(const std::vector<ivec> &arr, const std::string varname);

	/**
	 * \brief Write system state space into an hdf5 file.
	 *
	 * @param[in] ws System state space
	 * @param[in] varname The variable name to be written
	 */
	int write_state_space(const ivec &ws, const std::string varname);

	/**
	 * \brief Write control input into an hdf5 file.
	 *
	 * Control input is a set of uniformly sampled points from a compact control set \f$U\f$.
	 *
	 * @param[in] ugrid A set of control inputs
	 * @param[in] varname The variable name to be written
	 */
	int write_input_values(const grid &ugrid, const std::string varname);

	/**
	 * \brief Write control system settings.
	 * 
	 * The related data includes:
	 * - the sampling time `ts`
	 * - the state space `X`
	 * - the control set `U`
	 * @param[in] sys A system object
	 */
	template<typename S>
	void write_problem_setting(const S &sys) {//, const CSolver &sol) {
	    write_number<double>(sys._tau, "ts");
	    write_state_space(sys._workspace, "X");
	    write_input_values(sys._ugrid, "U");
	    // write_ivec_array(sol._goal, "G");
	    // write_ivec_array(sol._obs, "xobs");
	}
	
	/**
	 * \brief Write the BSolver controller into an hdf5 file.
	 *
	 * @param[in] sol the head of the controller
	 */
	int write_discrete_controller(HEAD *sol);

	/**
	 * \brief Write the DSolver controller into an hdf5 file.
	 *
	 * @param[in] dsol A DSolver object
	 */
	int write_discrete_controller(const DSolver &dsol);

	/**
	 * \brief Write a finite transition system into an hdf5 file.
	 *
	 * @param[in] fts A finite transition system
	 */
	int write_transitions(const fts &trans);

	/**
	 * \brief Write the winning graph of a controller generated by the BSolver into an hdf5 file.
	 *
	 * @param[in] patcher A Patcher object
	 */
	int write_winning_graph(const Patcher &patcher);

	/**
	 * \brief Write a CSolver tree leaves from a given node into an hdf5 file.
	 *
	 * @param[in] ctlr A SPtree object
	 * @param[in] ptrn The given tree node
	 */
	int write_sptree_leaves(const SPtree &ctlr, SPnode *ptrn);

	/**
	 * \brief Write a CSolver controller into an hdf5 file.
	 *
	 * @param[in] sol A CSolver object
	 */
	int write_sptree_controller(const CSolver &sol) {
	    write_sptree_leaves(sol._ctlr, sol._ctlr._root);
	    return 0;
	}

	/**
	 * \brief Read a number from an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in,out] ptrd The number to be read
	 * @param[in] varname The variable name in the file to be read
	 */
	template<typename T>
	int read_number(T *ptrd, const std::string varname);

	/**
	 * \brief Read an std::vector from an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in,out] arr The vector to be read
	 * @param[in] varname The variable name in the file to be read
	 */
	template<typename T>
	int read_array(std::vector<T> &arr, const std::string varname);

	/**
	 * \brief Read a 2d array from an hdf5 file.
	 *
	 * The supported types: double, int, long long, unsigned int, etc.
	 * @param[in,out] arr The vector to be read
	 * @param[in,out] len The length of each dimension
	 * @param[in] varname The variable name in the file to be read
	 */
	template<typename T>
	int read_2d_array(std::vector<T> &arr, size_t *len,
			  const std::string varname);

	/**
	 * \brief Read a finite transition system from an hdf5 file.
	 *
	 * @param[in,out] fts A finite transition system to be read
	 */
	int read_transitions(fts &trans);

	/**
	 * \brief Read the winning graph of a controller generated by the BSolver from an hdf5 file.
	 *
	 * @param[in,out] patcher A Patcher object to be read
	 */
	int read_winning_graph(Patcher &patcher);

	/**
	 * \brief Read the DSolver controller from an hdf5 file.
	 *
	 * @param[in,out] win The winning set (marked by 1)
	 * @param[in,out] lsctlr The least restrict controller (marked by 1)
	 * @param[in,out] cdim The length of each dimension
	 * @param[in,out] optctlr The optimal controller
	 * @param[in,out] value The optimal value
	 */
	int read_discrete_controller(boost::dynamic_bitset<> &win,
				     boost::dynamic_bitset<> &lsctlr,
				     size_t *cdims,
				     std::vector<size_t> &optctlr,
				     std::vector<double> &value);

	/**
	 * \brief Read the BSolver controller from an hdf5 file.
	 *
	 * @param[in,out] w_x0 The winning set
	 * @param[in,out] encode3 A state index mapping
	 * @param[in,out] nts_ctrlr The indexing table for finding the right entry in ctrl
	 * @param[in,out] ctrl The controller pair (q,u), where q is the %DBA state and u is the control input
	 * @param[in,out] q_prime The lookup table for q'
	 */
	int read_discrete_controller(std::vector<long long> &w_x0,
				     std::vector<long long> &encode3,
				     std::vector<NODE_POST> &nts_ctrlr,
				     std::vector<CTRL> &ctrl,
				     std::vector<int> &q_prime);

	/**
	 * \brief Read the CSolver controller from an hdf5 file.
	 *
	 * @param[in,out] pavings A 2d array of intervals
	 * @param[in,out] pdims The length of each dimension of pavings
	 * @param[in,out] tag The indicator of winning set (1 if the interval belongs to the winning set)
	 * @param[in,out] cntl The control table (row: # of pavings, column: # of control inputs, 1 if the corresponding control to the cell is permissible)
	 * @param[in,out] cdim The length of each dimension of cntl
	 */
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
