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

    /* Some template specialization */
    template<>
    H5::PredType get_datatype<int>() {return H5::PredType::NATIVE_INT;}
    template<>
    H5::PredType get_datatype<long long>() {return H5::PredType::NATIVE_LLONG;}
    template<>
    H5::PredType get_datatype<size_t>() {return H5::PredType::NATIVE_UINT64;}
    template<>
    H5::PredType get_datatype<float>() {return H5::PredType::NATIVE_FLOAT;}
    template<>
    H5::PredType get_datatype<double>() {return H5::PredType::NATIVE_DOUBLE;}
    template<>
    H5::PredType get_datatype<unsigned char>() {return H5::PredType::NATIVE_UCHAR;}
    

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


    int h5FileHandler::write_discrete_controller(const DSolver &dsol) {
	if(dsol._nw == 0) {
	    std::cout << "Winning set is empty. Writing control result is abandoned.\n";
	    return 1;
	}
	write_number<size_t>(dsol._nw, "nWin");
	std::vector<size_t> winset(dsol._nw);
	int iw = 0;
    	for(size_t i = 0; i < dsol._ts->_nx; ++i) {
    	    if(dsol._win[i]) {
    		winset[iw] = i;
    		++iw;
    	    }
    	}
	write_array<size_t>(winset, "WinSet");
	write_array<size_t>(dsol._optctlr, "OptCtlr");
	write_array<double>(dsol._value, "Value");
	size_t dim[2] = {dsol._ts->_nx, dsol._ts->_nu};
	unsigned char *data = new unsigned char[dim[0]*dim[1]];
	for(size_t i = 0; i < dim[0]*dim[1]; ++i) {
		data[i] = dsol._leastctlr[i];
	}
	write_2d_array<unsigned char>(data, dim, "LeastCtlr");
	delete[] data;
	return 0;
    }


    int h5FileHandler::write_transitions(const fts &trans) {
	H5::Group group(_h5file.createGroup("/Transitions"));
	write_number<size_t>(trans._nx, "/Transitions/Nx");
	write_number<size_t>(trans._nu, "/Transitions/Nu");
	write_number<size_t>(trans._ntrans, "/Transitions/Ntrans");
	write_array<size_t>(trans._idpost, "/Transitions/postID");
	write_array<int>(trans._npost, "/Transitions/postNum");
	write_array<size_t>(trans._ptrpost, "/Transitions/postAddr");
	write_array<size_t>(trans._idpre, "/Transitions/preID");
	write_array<int>(trans._npre, "/Transitions/preNum");
	write_array<size_t>(trans._ptrpre, "/Transitions/preAddr");
	write_array<double>(trans._cost, "/Transitions/cost");
	return 0;
    } //write_transitions


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
	size_t dim1[2] = {(size_t)nLeaf, 2*(size_t)dimState};
	write_2d_array<double>(pavings, dim1, "pavings");
	// hsize_t dim1[2];
	// dim1[0] = nLeaf;
	// dim1[1] = 2*dimState;
	// H5::DataSpace dataspace1(2, dim1);
	// H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	// H5::DataSet dsetPavings = _h5file.createDataSet(std::string("pavings"), datatype, dataspace1);
	// dsetPavings.write(pavings, datatype);
	delete[] pavings;

	/* write _ctlr to "ctlr" */
	size_t dim2[2] = {(size_t)nLeaf, (size_t)nInput};
	write_2d_array<unsigned char>(validu, dim2, "ctlr");
	// hsize_t dim2[2];
	// dim2[0] = nLeaf;
	// dim2[1] = nInput;
	// H5::DataSpace dataspace2(2, dim2);
	// H5::DataSet dsetCtlr = _h5file.createDataSet(std::string("ctlr"),
	// 					     H5::PredType::NATIVE_UCHAR, dataspace2);
	// dsetCtlr.write(validu, H5::PredType::NATIVE_UCHAR);
	delete[] validu;

	/* write _tag to "tag" */
	write_array<int>(tag, nLeaf, "tag");
	// hsize_t dim3[1];
	// dim3[0] = nLeaf;
	// H5::DataSpace dataspace3(1, dim3);
	// H5::DataSet dsetTag = _h5file.createDataSet(std::string("tag"),
	// 					    H5::PredType::NATIVE_INT, dataspace3);
	// dsetTag.write(tag, H5::PredType::NATIVE_INT);
	delete[] tag;

	return 0;
    }


    int h5FileHandler::read_transitions(fts &trans) {
	read_number<size_t>(&trans._nx, "/Transitions/Nx");
	read_number<size_t>(&trans._nu, "/Transitions/Nu");
	read_number<size_t>(&trans._ntrans, "/Transitions/Ntrans");
	read_array<size_t>(trans._idpost, "/Transitions/postID");
	read_array<int>(trans._npost, "/Transitions/postNum");
	read_array<size_t>(trans._ptrpost, "/Transitions/postAddr");
	read_array<size_t>(trans._idpre, "/Transitions/preID");
	read_array<int>(trans._npre, "/Transitions/preNum");
	read_array<size_t>(trans._ptrpre, "/Transitions/preAddr");
	read_array<double>(trans._cost, "/Transitions/cost");
	return 0;
    }//read_transitions
    

    int h5FileHandler::read_sptree_controller(std::vector<double> &pavings, size_t *pdims,
					      std::vector<int> &tag,
					      boost::dynamic_bitset<> &cntl, size_t *cdims) {
	read_2d_array<double>(pavings, pdims, "pavings");
	read_array<int>(tag, "tag");
	std::vector<unsigned char> c;
	read_2d_array<unsigned char>(c, cdims, "ctlr");
	cntl.resize(c.size(), false);
	for(size_t i = 0; i < c.size(); ++i) {
	    if(c[i]>0)
		cntl[i] = true;
	}
	return 0;
    }

}//namespace rocs
