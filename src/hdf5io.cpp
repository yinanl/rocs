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
	    std::cout << "h5FileHandler::write_ivec_array:"
		      << varname << " is empty.\n";
	    return 0;
	}
	hsize_t dim[3] = {(hsize_t)arr[0].getdim(), 2, arr.size()};
	H5::DataSpace dataspace(3, dim);
	double *data = new double[dim[0]*dim[1]*dim[2]];
	for(hsize_t i = 0; i < dim[2]; ++i) {
	    for(hsize_t j = 0; j < dim[0]; ++j) {// j row of the i element in arr, 2 cols
		*(data+j*dim[2]*dim[1]+i) = arr[i][j].getinf();
		*(data+j*dim[2]*dim[1]+dim[2]+i) = arr[i][j].getsup();
	    }
	}
	// boost::multi_array<double, 3>  data(boost::extents[dim[0]][dim[1]][dim[2]]);
	// for(hsize_t i = 0; i < dim[2]; ++i) {
	//     for(hsize_t j = 0; j < dim[0]; ++j) {// j row of the i element in arr, 2 cols
	// 	data[j][0][i] = arr[i][j].getinf();
	// 	data[j][1][i] = arr[i][j].getsup();
	//     }
	// }
	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	H5::DataSet dataset = _h5file.createDataSet(varname, datatype, dataspace);
	dataset.write(data, datatype);
	delete[] data;
	return 0;
    }


    int h5FileHandler::write_discrete_controller(HEAD *sol) {
    	CTRLR ctrlr = sol->ctrlr;
    	if(!ctrlr.ctrl) {
    	    std::cout << "No controller generated. Writing control result is abandoned.\n";
    	    return 1;
    	}
    	SCOPE *scope = sol->nts_product.scope;
    	long long *decode2 = sol->decode2;
    	long long *decode1 = sol->decode1;
    	long long *encode3 = sol->encode3;
    	/* W_X0: the winning set */
    	long long *win = new long long[ctrlr.w0];
    	long long i = 0;
    	for(SCOPE *p = scope; i < ctrlr.w0; ++i) {
    	    p += p->next;
    	    win[i] = *(decode1+*(decode2+p->x)/sol->dba.n);
    	}
    	write_array<long long>(win, ctrlr.w0, "WinSet");
    	delete[] win;
    	/* encode3: n0->n0_2 */
    	long long *en3 = new long long[sol->nts_pre.n];
    	i = 0;
    	for(long long *p1=encode3, *p2=p1+sol->nts_pre.n; p1<p2; ++p1, ++i)
    	    en3[i] = *p1;
    	write_array<long long>(en3, sol->nts_pre.n, "encode3");
    	delete[] en3;
    	/* nts_ctrlr2 */
    	long long *ctrlr2 = new long long[ctrlr.n*3];
    	size_t cdims[2] = {(size_t)ctrlr.n, 3};
    	i = 0;
    	for(NODE_POST *p1=ctrlr.graph, *p2=p1+ctrlr.n; p1<p2; ++p1, ++i) {
    	    ctrlr2[3*i] = p1->num_a;
    	    ctrlr2[3*i+1] = p1->label;
    	    ctrlr2[3*i+2] = p1->pos;
    	}
    	write_2d_array<long long>(ctrlr2, cdims, "nts_ctrlr");
    	delete[] ctrlr2;
    	/* ctrlr.ctrl */
    	int *ctrldata = new int[ctrlr.wp*2];
    	size_t pdims[2] = {(size_t)ctrlr.wp, 2};
    	i = 0;
    	for(CTRL *p1=ctrlr.ctrl, *p2=p1+ctrlr.wp; p1<p2; ++p1, ++i) {
    	    ctrldata[2*i] = p1->q;
    	    ctrldata[2*i+1] = p1->u;
    	}
    	write_2d_array<int>(ctrldata, pdims, "OptCtlr");
    	delete[] ctrldata;
    	/* q_prime */
    	int qsize = sol->dba.q_prime_size;
    	int *q_prime = new int[qsize];
    	i = 0;
    	for(int *p1=sol->dba.q_prime, *p2=p1+qsize; p1<p2; ++p1, ++i)
    	    q_prime[i] = *p1;
    	write_array<int>(q_prime, qsize, "q_prime");
    	delete[] q_prime;

    	return 0;
    }

    int h5FileHandler::write_discrete_controller(const DSolver &dsol) {
	if(!dsol._nw) {
	    std::cout << "Winning set is empty. Writing control result is abandoned.\n";
	    return 1;
	}
	
	write_number<size_t>(dsol._nw, "nWin");
	unsigned char *win = new unsigned char[dsol._ts->_nx];
	for(size_t i = 0; i < dsol._ts->_nx; ++i) {
		win[i] = dsol._win[i];
	}
	write_array<unsigned char>(win, dsol._ts->_nx, "Win");
	delete[] win;
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


    int h5FileHandler::write_winning_graph(const Patcher &patcher) {
	H5::Group group(_h5file.createGroup("/WinGraph"));
	write_number<size_t>(patcher._nwin, "/WinGraph/Nx");
	write_number<size_t>(patcher._na, "/WinGraph/Nu");
	write_array<size_t>(patcher._winfts._idpost, "/WinGraph/postID");
	write_array<int>(patcher._winfts._npost, "/WinGraph/postNum");
	write_array<size_t>(patcher._winfts._ptrpost, "/WinGraph/postAddr");
	write_array<size_t>(patcher._winfts._idpre, "/WinGraph/preID");
	write_array<int>(patcher._winfts._npre, "/WinGraph/preNum");
	write_array<size_t>(patcher._winfts._ptrpre, "/WinGraph/preAddr");
	write_array<long long>(patcher._reachstep, "/reachSteps");
	write_array<long long>(patcher._idmap, "/idMap");
	write_array<long long>(patcher._encode, "/encode");
	write_array<long long>(patcher._decode, "/decode");
	return 0;
    } //end write_winning_graph
    

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
	delete[] pavings;

	/* write _ctlr to "ctlr" */
	size_t dim2[2] = {(size_t)nLeaf, (size_t)nInput};
	write_2d_array<unsigned char>(validu, dim2, "ctlr");
	delete[] validu;

	/* write _tag to "tag" */
	write_array<int>(tag, nLeaf, "tag");
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
    
    
    int h5FileHandler::read_winning_graph(Patcher &patcher) {
	read_number<size_t>(&patcher._nwin, "/WinGraph/Nx");
	read_number<size_t>(&patcher._na, "/WinGraph/Nu");
	read_array<size_t>(patcher._winfts._idpost, "/WinGraph/postID");
	read_array<int>(patcher._winfts._npost, "/WinGraph/postNum");
	read_array<size_t>(patcher._winfts._ptrpost, "/WinGraph/postAddr");
	read_array<size_t>(patcher._winfts._idpre, "/WinGraph/preID");
	read_array<int>(patcher._winfts._npre, "/WinGraph/preNum");
	read_array<size_t>(patcher._winfts._ptrpre, "/WinGraph/preAddr");
	read_array<long long>(patcher._reachstep, "/reachSteps");
	read_array<long long>(patcher._idmap, "idMap");
	read_array<long long>(patcher._encode, "/encode");
	read_array<long long>(patcher._decode, "/decode");
	return 0;
    } //end save_winning_graph_to_h5


    int h5FileHandler::read_discrete_controller(std::vector<long long> &w_x0,
    						std::vector<long long> &encode3,
    						std::vector<NODE_POST> &nts_ctrlr,
    						std::vector<CTRL> &ctrl,
    						std::vector<int> &q_prime) {
    	read_array<long long>(w_x0, "WinSet");
    	read_array<long long>(encode3, "encode3");
    	read_array<int>(q_prime, "q_prime");
    	std::vector<long long> ctrlr2;
    	size_t cdims[2];
    	read_2d_array<long long>(ctrlr2, cdims, "nts_ctrlr");
    	nts_ctrlr.resize(cdims[0]);
    	for(size_t i = 0; i < nts_ctrlr.size(); ++i) {
    	    nts_ctrlr[i].num_a = (int)ctrlr2[i*cdims[1]];
    	    nts_ctrlr[i].label = (int)ctrlr2[i*cdims[1]+1];
    	    nts_ctrlr[i].pos = ctrlr2[i*cdims[1]+2];
    	}
    	std::vector<int> ctrlpair;
    	size_t pdims[2];
	read_2d_array<int>(ctrlpair, pdims, "OptCtlr");
    	ctrl.resize(pdims[0]);
    	for(size_t i = 0; i < ctrl.size(); ++i) {
    	    ctrl[i].q = ctrlpair[i*pdims[1]];
    	    ctrl[i].u = ctrlpair[i*pdims[1]+1];
    	}

    	return 0;
    }


    int h5FileHandler::read_discrete_controller(boost::dynamic_bitset<> &win,
						boost::dynamic_bitset<> &lsctlr,
						size_t *cdims,
						std::vector<size_t> &optctlr,
						std::vector<double> &value) {
	std::vector<unsigned char> w, ls;
	read_array<unsigned char>(w, "Win");
	read_2d_array<unsigned char>(ls, cdims, "LeastCtlr");
	win.resize(w.size(), false);
	for(size_t i = 0; i < w.size(); ++i)
	    if(w[i])
		win[i] = true;
	lsctlr.resize(ls.size(), false);
	for(size_t i = 0; i < ls.size(); ++i)
	    if(ls[i])
		lsctlr[i] = true;
	read_array<size_t>(optctlr, "OptCtlr");
	read_array<double>(value, "Value");
	return 0;
    }


    int h5FileHandler::read_sptree_controller(std::vector<double> &pavings, size_t *pdims,
    					      std::vector<int> &tag,
    					      boost::dynamic_bitset<> &cntl, size_t *cdims) {
    	read_2d_array<double>(pavings, pdims, "pavings");
    	read_array<int>(tag, "tag");
    	std::vector<unsigned char> c;
    	read_2d_array<unsigned char>(c, cdims, "ctlr");
    	cntl.resize(c.size(), false);
    	for(size_t i = 0; i < c.size(); ++i) {
    	    if(c[i])
    		cntl[i] = true;
    	}
    	return 0;
    }

}//namespace rocs
