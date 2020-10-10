/**
 *  hdf5io.cpp
 *
 *  The source file of input/output classes to .h5 files.
 *
 *  Created by Yinan Li on Aug 16, 2020.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#include "hdf5io.h"

namespace rocs {

    int h5FileHandler::write_real_number(const double x, const std::string varname) {
	hsize_t dim[1] = {1};
	H5::DataSpace dataspace(1, dim);
	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	dataset.write(&x, datatype);
	return 0;
    }


    int h5FileHandler::write_real_array(const std::vector<double> &x,
					const std::string varname) {
	hsize_t dim[1]= {x.size()};
	H5::DataSpace dataspace(1, dim);
	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	
	/* Setup chunk size and compression filter */
	hsize_t chunkSize = dim[0] / 10;
	hsize_t chunkdim[1] = {chunkSize};
	H5::DSetCreatPropList plist;
	plist.setChunk(1, chunkdim);
	plist.setDeflate(6);
	
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace, plist);
	dataset.write(&(x[0]), datatype);
	// try {
	//     H5::hsize_t dim[]{x.size()};
	//     H5::DataSpace dataspace(1, dim);
	//     H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	//     H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	//     dataset.write(&(x[0]), datatype);
	// }
	// // catch failure caused by the H5File operations
	// catch( H5::FileIException error ) {
	//     error.printError();
	//     return -1;
	// }

	return 0;
    }


    int h5FileHandler::write_2d_real_array(const std::vector< std::vector<double> > &x,
					   const std::string varname) {
	hsize_t dim[2];
	dim[0]= x.size();
	dim[1] = x[0].size();
	H5::DataSpace dataspace(2, dim);
	/* create the buffer data for writing */
	double *data = new double[dim[0]*dim[1]];
	for(hsize_t i = 0; i < dim[0]; ++i) {
	    for(hsize_t j = 0; j < dim[1]; ++j)
		*(data+dim[1]*i+j) = x[i][j];
	}
	
	/* write to h5 dataset */
	hsize_t chunkSize = dim[0] / 10;
	hsize_t chunkdim[2] = {chunkSize, dim[1]};
	H5::DSetCreatPropList plist;
	plist.setChunk(2, chunkdim);
	plist.setDeflate(6);
	
	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace, plist);
	dataset.write(data, datatype);
	delete[] data;
	return 0;
    }


    int h5FileHandler::write_uint_number(const size_t d, const std::string varname) {
	hsize_t dim[1] = {1};
	H5::DataSpace dataspace(1, dim);
	H5::IntType datatype(H5::PredType::NATIVE_UINT64);
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	dataset.write(&d, datatype);
	return 0;
    }


    int h5FileHandler::write_uint_array(const std::vector<size_t> &arr,
					const std::string varname) {
	hsize_t dim[1]= {arr.size()};
	H5::DataSpace dataspace(1, dim);

	hsize_t chunkSize = dim[0] / 10;
	hsize_t chunkdim[1] = {chunkSize};
	H5::DSetCreatPropList plist;
	plist.setChunk(1, chunkdim);
	plist.setDeflate(6);
	
	H5::IntType datatype(H5::PredType::NATIVE_UINT64);
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace, plist);
	dataset.write(&(arr[0]), datatype);
	return 0;
    }


    int h5FileHandler::write_state_space(const ivec &ws, const std::string varname) {
        int n = ws.getdim();
	hsize_t dim[2];
	dim[0] = ws.getdim();
	dim[1] = 2;
	H5::DataSpace dataspace(2, dim);
	/* create the buffer data for writing */
	double *data = new double[n*2];
	for(int i = 0; i < n; ++i) {
	    *(data+2*i) = ws[i].getinf();
	    *(data+2*i+1) = ws[i].getsup();
	}
	/* write to h5 dataset */
	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	dataset.write(data, datatype);

	delete[] data;
	return 0;
    }


    int h5FileHandler::write_input_values(const grid &ugrid, const std::string varname) {
	if(ugrid._nv > 1) {
	    hsize_t dim[2];
	    dim[0] = ugrid._nv;
	    dim[1] = ugrid._dim;
	    H5::DataSpace dataspace(2, dim);
	    double *data = new double[dim[0]*dim[1]];
	    if(!ugrid._data.empty()) {
		for(hsize_t r = 0; r < dim[0]; ++r)
		    for(hsize_t c = 0; c < dim[1]; ++c)
			*(data+dim[1]*r+c) = ugrid._data[r][c];
	    } else {//control inputs are: 1,...,ugrid._nv (e.g, switched modes)
		for(hsize_t r = 0; r < dim[0]; ++r)
		    *(data+r) = double(r+1);
	    }
	    H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	    H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	    dataset.write(data, datatype);
	    delete[] data;
	} else {
	    std::cout << "h5FileHandler::write_input_values: Write input values abandoned due to single input value.\n";
	}
	return 0;
    }


    int h5FileHandler::write_ivec_array(const std::vector<ivec> &arr,
					const std::string varname) {
	if(arr.empty()) {
	    std::cout << "h5FileHandler::write_ivec_array:Input interval vector array is empty.\n";
	    return 0;
	}
	hsize_t dim[3];
	dim[0] = arr[0].getdim();
	dim[1] = 2;
	dim[2] = arr.size();
	H5::DataSpace dataspace(3, dim);
	// hsize_t nData = dim[0]*dim[1]*dim[2];
	double *data = new double[dim[0]*dim[1]*dim[2]];
	for(hsize_t i = 0; i < dim[2]; ++i)
	    for(hsize_t j = 0; j < dim[0]; ++j) {// j row of the i element in arr, 2 cols
		*(data+i*dim[0]*dim[1]+j*dim[1]) = arr[i][j].getinf();
		*(data+i*dim[0]*dim[1]+j*dim[1] + 1) = arr[i][j].getsup();
	    }
	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	dataset.write(data, datatype);
	delete[] data;
	return 0;
    }


    int h5FileHandler::write_sptree_leaves(const SPtree &ctlr, SPnode *ptrn) {
	if (ptrn == NULL) {
	    std::cout << "h5FileHandler::write_sptree_leaves: Input sptree root is empty." << std::endl;
	    return 0;
	}
	int nLeaf = ctlr.leafcount(ptrn);
	int dimState = ptrn->_box.getdim();
	int nInput = ptrn->_cntl.size();
	double *pavings = new double[nLeaf*dimState*2];
	unsigned char *validu = new unsigned char[nLeaf*nInput];
	int *tag = new int[nLeaf];

	std::stack<SPnode*> stk;
	SPnode *current = ptrn;
	stk.push(current);
	size_t j = 0;
	while (!stk.empty()) {
	    current = stk.top();
	    stk.pop();
	    if (ctlr.isleaf(current)) { //write leaves
		for (int i = 0; i < dimState; ++ i) {
		    pavings[j*2*dimState+2*i] = (current->_box)[i].getinf();
		    pavings[j*2*dimState+2*i+1] = (current->_box)[i].getsup();
		}
		for (int k = 0; k < nInput; ++ k) {
		    validu[j*nInput+k] = current->_cntl[k];
		}
		tag[j] = current->_tag;
		j ++;
	    }
	    if (current->_left)
		stk.push(current->_left);
	    if (current->_right)
		stk.push(current->_right);
	}

	/* write leaf nodes of a SPtree to "pavings" */
	hsize_t dim1[2];
	dim1[0] = nLeaf;
	dim1[1] = 2*dimState;
	H5::DataSpace dataspace1(2, dim1);
	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	H5::DataSet dsetPavings = _h5file.createDataSet(std::string("pavings"), datatype, dataspace1);
	dsetPavings.write(pavings, datatype);
	delete[] pavings;

	/* write _ctlr to "ctlr" */
	hsize_t dim2[2];
	dim2[0] = nLeaf;
	dim2[1] = nInput;
	H5::DataSpace dataspace2(2, dim2);
	H5::DataSet dsetCtlr = _h5file.createDataSet(std::string("ctlr"),
						     H5::PredType::NATIVE_UCHAR, dataspace2);
	dsetCtlr.write(validu, H5::PredType::NATIVE_UCHAR);
	delete[] validu;

	/* write _tag to "tag" */
	hsize_t dim3[1];
	dim3[0] = nLeaf;
	H5::DataSpace dataspace3(1, dim3);
	H5::DataSet dsetTag = _h5file.createDataSet(std::string("tag"),
						    H5::PredType::NATIVE_INT, dataspace3);
	dsetTag.write(tag, H5::PredType::NATIVE_INT);
	delete[] tag;

	return 0;
    }


    int h5FileHandler::write_transitions(const fts &trans) {
	H5::Group group(_h5file.createGroup("/Transitions"));
	write_uint_number(trans._nx, "/Transitions/Nx");
	write_uint_number(trans._nu, "/Transitions/Nu");
	write_uint_number(trans._ntrans, "/Transitions/Ntrans");
	write_uint_array(trans._idpost, "/Transitions/postID");
	write_uint_array(trans._npost, "/Transitions/postNum");
	write_uint_array(trans._ptrpost, "/Transitions/postAddr");
	write_uint_array(trans._idpre, "/Transitions/preID");
	write_uint_array(trans._npre, "/Transitions/preNum");
	write_uint_array(trans._ptrpre, "/Transitions/preAddr");
	write_real_array(trans._cost, "/Transitions/cost");
	return 0;
    } //write_transitions


    int h5FileHandler::read_uint_number(size_t *ptrd, const std::string varname) {
	size_t data[1];
	H5::DataSet dataset = _h5file.openDataSet(varname);
	H5::DataSpace dataspace = dataset.getSpace();
	// int rank_dnx = dnxspace.getSimpleExtentNdims();
	hsize_t dim[1] = {1};
	H5::DataSpace memspace(1, dim);
	dataset.read(data, H5::PredType::NATIVE_UINT64, memspace, dataspace);
	*ptrd = data[0];
	return 0;
    }


    int h5FileHandler::read_uint_array(std::vector<size_t> &arr,
				       const std::string varname) {
	H5::DataSet dataset = _h5file.openDataSet(varname);
	H5::DataSpace dataspace = dataset.getSpace();
	int rank = dataspace.getSimpleExtentNdims();
	hsize_t dim[1];
	rank = dataspace.getSimpleExtentDims(dim);
	H5::DataSpace memspace(1, dim);

	// /* Get the create property list of the dataset */
	// H5::DSetCreatPropList plist(dataset.getCreatePlist());
	// int numfilt = plist.getNfilters();
	
	// std::cout << "dim[0]=" << dim[0] << '\n';
	size_t *data = new size_t[dim[0]];
	dataset.read(data, H5::PredType::NATIVE_UINT64, memspace, dataspace);
	arr.resize(dim[0]);
	for(hsize_t i = 0; i < dim[0]; ++i)
	    arr[i] = data[i];
	delete[] data;
	return 0;
    }


    int h5FileHandler::read_real_array(std::vector<double> &arr,
				       const std::string varname) {
	H5::DataSet dataset = _h5file.openDataSet(varname);
	H5::DataSpace dataspace = dataset.getSpace();
	int rank = dataspace.getSimpleExtentNdims();
	hsize_t dim[1];
	rank = dataspace.getSimpleExtentDims(dim);
	H5::DataSpace memspace(1, dim);
	double *data = new double[dim[0]];
	dataset.read(data, H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
	arr.resize(dim[0]);
	for(hsize_t i = 0; i < dim[0]; ++i)
	    arr[i] = data[i];
	delete[] data;
	return 0;
    }


    int h5FileHandler::read_transitions(fts &trans) {
	read_uint_number(&trans._nx, "/Transitions/Nx");
	read_uint_number(&trans._nu, "/Transitions/Nu");
	read_uint_number(&trans._ntrans, "/Transitions/Ntrans");
	read_uint_array(trans._idpost, "/Transitions/postID");
	read_uint_array(trans._npost, "/Transitions/postNum");
	read_uint_array(trans._ptrpost, "/Transitions/postAddr");
	read_uint_array(trans._idpre, "/Transitions/preID");
	read_uint_array(trans._npre, "/Transitions/preNum");
	read_uint_array(trans._ptrpre, "/Transitions/preAddr");
	read_real_array(trans._cost, "/Transitions/cost");
	return 0;
    }//read_transitions


}//namespace rocs
