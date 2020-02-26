/**
 *  matlabio.cpp
 *
 *  The source file of matlab input/output classes.
 *
 *  Created by Yinan Li on May 10, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include "matlabio.h"

namespace rocs {

    
    void matWriter::write_real_number(const double x, const char *varname) {
	if (_pmat == NULL) {open();}
	mxArray *ptrx = mxCreateDoubleScalar(x);
	if (matPutVariable(_pmat, varname, ptrx)!=0) {
	    std::cout << "matWriter::write_real_number: Error writing the real number.\n";
	}
	mxDestroyArray(ptrx);
    }//matWriter::write_real_number


    void matWriter::write_state_space(const ivec &ws, const char *varname) {
	if (_pmat == NULL) {open();}

	int n = ws.getdim();
	mxArray *x = mxCreateDoubleMatrix(n, 2, mxREAL);
	double *ptrx = mxGetPr(x);
	for (int i = 0; i < n; ++i) {
	    ptrx[i] = ws[i].getinf();
	    ptrx[i+n] = ws[i].getsup();
	}
      
	if (matPutVariable(_pmat, varname, x)!=0) {
	    std::cout << "matWriter::write_state_space: Error writing state space.\n";
	}
	mxDestroyArray(x);
    }//matWriter::write_state_space


    void matWriter::write_input_values(const grid &ugrid, const char *varname) {
	if (ugrid._nv > 1) {
	    if (_pmat == NULL) {open();}
	    mxArray *u = mxCreateDoubleMatrix(ugrid._nv, ugrid._dim, mxREAL);
	    double *ptru = mxGetPr(u);
	    if (!ugrid._data.empty()) {
		for (int r = 0; r < ugrid._nv; ++r) 
		    for (int c = 0; c < ugrid._dim; ++c)
			ptru[r + c*ugrid._nv] = ugrid._data[r][c];
	    } else {
		for (int r = 0; r < ugrid._nv; ++r)
		    ptru[r] = double(r+1);
	    }
      
	    if (matPutVariable(_pmat, varname, u)!=0) {
		std::cout << "matWriter::write_input_values: Error writing input space.\n";
	    }
	    mxDestroyArray(u);
	} else
	    std::cout << "matWriter: Write input values abandoned due to single input value.\n";
    }//matWriter::write_input_values
    

    void matWriter::write_ivec_array(const std::vector<ivec> &arr, const char *varname) {
	if (!arr.empty()) {
	    mwSize dims[3] = {(mwSize)arr[0].getdim(), 2, arr.size()};
	    mxArray *us = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	    double *ptrus = mxGetPr(us);
	    
	    for (int r = 0; r < dims[0]; ++r) {
		for (int p = 0; p < dims[2]; ++p) {
		    ivec obstacle = arr[p];
		    ptrus[r + p*dims[0]*dims[1]] = obstacle[r].getinf();
		    ptrus[r + dims[0] + p*dims[0]*dims[1]] = obstacle[r].getsup();
		}
	    }
	    
	    if (matPutVariable(_pmat, varname, us)!=0) {
		std::cout << "matWriter::write_unsafe_area: Error writing unsafe area.\n";
	    }
	    mxDestroyArray(us);
	}
    }//matWriter::write_ivec_array


    void matWriter::write_sptree_leaves(const SPtree &c, SPnode *ptrn) {
	if (ptrn == NULL) {
	    std::cout << "matWriter: Nothing Generated." << std::endl;
	    return;
	}
	int nleaf = c.leafcount(ptrn); //count rows of output matrices
    
	/* define outputs */
	int nx = ptrn->_box.getdim();
	mxArray *matx = mxCreateDoubleMatrix(nleaf, 2 * nx, mxREAL);
	double *ptrx = mxGetPr(matx);
    
	int nu = ptrn->_cntl.size();
	mxArray *matu = mxCreateLogicalMatrix(nleaf, nu);
	mxLogical *ptru = mxGetLogicals(matu);
    
	mxArray *mattag = mxCreateDoubleMatrix(nleaf, 1, mxREAL);
	double *ptrtag = mxGetPr(mattag);
	
	/* write data to buffers using stacks */
	std::stack<SPnode*> stk;
	SPnode *current = ptrn;
	stk.push(current);
	std::vector<double> lower, upper;
	size_t j = 0;
	while (!stk.empty()) {
	    current = stk.top();
	    stk.pop();
	    
	    if (c.isleaf(current)) { //write leaves
		lower = current->_box.getinf();
		upper = current->_box.getsup();
		for (int i = 0; i < nx; ++ i) {
		    ptrx[j + 2 * i * nleaf] = lower[i];
		    ptrx[j + (2 * i + 1) * nleaf] = upper[i];
		}

		for (int k = 0; k < nu; ++ k) {
		    ptru[j + k * nleaf] = current->_cntl[k];
		}
		ptrtag[j] = current->_tag;
		j ++;
	    }

	    if (current->_left)
		stk.push(current->_left);

	    if (current->_right)
		stk.push(current->_right);
	}
	
	/* write to variables in .mat */
	if (matPutVariable(_pmat, "pavings", matx) ||
	    matPutVariable(_pmat, "ctlr", matu) ||
	    matPutVariable(_pmat, "tag", mattag) != 0) {
	    std::cout << "matWriter::write_sptree_leaves: Error outputting variables" << std::endl;
	    return;
	}
    
	mxDestroyArray(matx);
	mxDestroyArray(matu);
	mxDestroyArray(mattag);
    }//matWriter::write_sptree_leaves
    

    void matWriter::write_winset_boundary(const double eta[], const CSolver &sol, const char *varname) {
	if (_pmat == NULL) {open();}

	/* Generate a grid for the state space */
	ivec ws = sol._ctlr._root->_box;
	grid xg(sol._xdim, eta, ws);
	xg.gridding();
	
	if (xg._data.empty()) {
	    std::cout << "matWriter::write_winset_boundary: Error in grid generation.\n";
	    return;
	    
	} else {
	    std::vector<bool> win(xg._nv);
	    SPnode *q;
	    for (int i = 0; i < xg._nv; ++i) {
		q = sol._ctlr._root;
		while (q->_tag == 2 && !sol._ctlr.isleaf(q)) {
		    if (q->_left->_box.isin(xg._data[i]))
			q = q->_left;
		    else
			q = q->_right;
		}
		if (q->_tag == 1) {//the ith grid is inside the winning set 
		    win[i] = true;
		}
	    }
	    /* Extract the boundary by checking if neighbours are all in the winning set */
	    write_boundary(xg, win, varname);
	}
    
    }//matWriter::write_winset_boundary


    void matWriter::write_boundary(grid &g, const std::vector<bool> &win, const char *varname) {
	if (_pmat == NULL) {open();}

	std::vector<size_t> nbs;
	std::vector<size_t> boundary;
	for (int i = 0; i < g._nv; ++i) {
	    if (win[i] == true) {
		nbs = g.neighbours(i);
		for (int j = 0; j < nbs.size(); ++j) {
		    if (win[nbs[j]] == false) {
			boundary.push_back(nbs[j]);
			break;
		    }
		}
	    }
	}
	/* write to a mat file */
	int bdsize = boundary.size();
	mxArray *matbd = mxCreateDoubleMatrix(bdsize, g._dim, mxREAL);
	double *ptrbd = mxGetPr(matbd);
	for (size_t i = 0; i < bdsize; ++ i) {
	    for (size_t j = 0; j < g._dim; ++ j) {
		ptrbd[i + bdsize * j] = g._data[boundary[i]][j];
	    }
	}
	    
	matPutVariable(_pmat, varname, matbd);
	mxDestroyArray(matbd);
    }//matWriter::write_boundary


    void matWriter::write_controller_serialized(const CSolver &sol) {

	if (sol._ctlr.isempty()) {
	    std::cout << "matWriter: No Controller Generated.\n";
	    return;
	}

	if (_pmat == NULL) {open();}

	size_t H = sol._ctlr.height(sol._ctlr._root);
	mwSize L = sol._ctlr.leafcount(sol._ctlr._root);
	mwSize M = std::pow(2,H)-1;
    
	int nx = sol._xdim;
    
	/* define outputs */
	mxArray *mats = mxCreateDoubleMatrix(M, 2*nx+1, mxREAL);  // Btree
	double *ptrs = mxGetPr(mats);
	mxArray *matu = mxCreateLogicalMatrix(L, sol._nu);  // control table
	mxLogical *ptru = mxGetLogicals(matu);
	mxArray *mati = mxCreateDoubleMatrix(L, 2, mxREAL);  // indices & tags of control entries
	double *ptri = mxGetPr(mati);

	/* write data to buffers using stacks */
	std::queue<SPnode*> qnode;
	std::queue<size_t> qid;
	qnode.push(sol._ctlr._root);
	qid.push(0);

	std::vector<double> lower, upper;
	SPnode *current;
	size_t i;
	size_t lc = 0;
	while(!qnode.empty()) {

	    current = qnode.front();
	    qnode.pop();
	    i = qid.front();
	    qid.pop();

	    ptrs[i] = current->_split;  
	
	    lower = current->_box.getinf();  
	    upper = current->_box.getsup();
	    for (int j = 0; j < nx; ++j) {
		
		ptrs[i + (2*j+1)*M] = lower[j];
		ptrs[i + (2*j+2)*M] = upper[j];
	    }

	    if (sol._ctlr.isleaf(current)) {
	    
		for (int k = 0; k < sol._nu; ++k) { // control array

		    ptru[lc + k*L] = current->_cntl[k];
		}
		ptri[lc] = i + 1;
		ptri[lc + L] = current->_tag;

		++lc;
	    }
	    else {

		if (current->_left != NULL) {

		    qnode.push(current->_left);
		    qid.push(2*i+1);
		}

		if (current->_right != NULL) {

		    qnode.push(current->_right);
		    qid.push(2*i+2);
		}
	    } // end if
	} // end while

	assert(i < M);
	assert(lc-1 < L);

	/* write to variables in .mat */
	if (matPutVariable(_pmat, "ctree", mats) ||
	    matPutVariable(_pmat, "cindex", mati) ||
	    matPutVariable(_pmat, "cvalue", matu) != 0) {
	
	    std::cout << "matWriter::write_controller_serialized: Error outputting variables.\n";
	    return;
	}
    
	mxDestroyArray(mats);
	mxDestroyArray(matu);
	mxDestroyArray(mati);
    }//matWriter::write_controller_serialized


    void matWriter::write_transitions(const fts &trans,
				      const char* vtranspost, const char* vpost, const char* vptrpost,
				      const char* vtranspre, const char* vpre, const char* vptrpre) {
	if (_pmat == NULL) {open();}

	mxArray *mtranspost = mxCreateDoubleMatrix(trans._ntrans, 5, mxREAL);
	double *ptranspost = mxGetPr(mtranspost);
	mxArray *matnpost = mxCreateDoubleMatrix(trans._nx*trans._nu, 1, mxREAL);
	double *ptrnpost = mxGetPr(matnpost);
	mxArray *matppost = mxCreateDoubleMatrix(trans._nx*trans._nu, 1, mxREAL);
	double *ptrppost = mxGetPr(matppost);
    
	mxArray *mtranspre = mxCreateDoubleMatrix(trans._ntrans, 3, mxREAL);
	double *ptranspre = mxGetPr(mtranspre);
	mxArray *matnpre = mxCreateDoubleMatrix(trans._nx*trans._nu, 1, mxREAL);
	double *ptrnpre = mxGetPr(matnpre);
	mxArray *matppre = mxCreateDoubleMatrix(trans._nx*trans._nu, 1, mxREAL);
	double *ptrppre = mxGetPr(matppre);
    
	size_t i = 0;
	size_t j = 0;
	for (size_t row = 0; row < trans._nx; ++row) {
	    for (size_t col = 0; col < trans._nu; ++col) {

		for (int ip = 0; ip < trans._npost[row*trans._nu + col]; ++ip) {
		    ptranspost[i] = row;
		    ptranspost[i + trans._ntrans * 1] = trans._idpost[trans._ptrpost[row*trans._nu+col] + ip];
		    ptranspost[i + trans._ntrans * 2] = col;
		    ptranspost[i + trans._ntrans * 3] = 1;
		    ptranspost[i + trans._ntrans * 4] = trans._cost[trans._ptrpost[row*trans._nu+col] + ip];
		
		    i ++;
		}

		for (int jp = 0; jp < trans._npre[row*trans._nu + col]; ++jp) {
		    ptranspre[j] = row;
		    ptranspre[j + trans._ntrans * 1] = trans._idpre[trans._ptrpre[row*trans._nu+col] + jp];
		    ptranspre[j + trans._ntrans * 2] = col;	    
		    j ++;
		}

	    }
	}

	for (size_t j = 0; j < trans._nx*trans._nu; ++j) {
	    ptrnpost[j] = trans._npost[j];
	    ptrppost[j] = trans._ptrpost[j];
	    ptrnpre[j] = trans._npre[j];
	    ptrppre[j] = trans._ptrpre[j];
	}

	if(matPutVariable(_pmat, vtranspost, mtranspost) != 0) {
	    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
	}

	if(matPutVariable(_pmat, vpost, matnpost) != 0) {
	    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
	}
	if(matPutVariable(_pmat, vptrpost, matppost) != 0) {
	    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
	}
	if(matPutVariable(_pmat, vtranspre, mtranspre) != 0) {
	    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
	}
	if(matPutVariable(_pmat, vpre, matnpre) != 0) {
	    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
	}
	if(matPutVariable(_pmat, vptrpre, matppre) != 0) {
	    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
	}

	mxDestroyArray(mtranspost);
	mxDestroyArray(matnpost);
	mxDestroyArray(matppost);

	mxDestroyArray(mtranspre);
	mxDestroyArray(matnpre);
	mxDestroyArray(matppre);
    }//matWriter::write_transitions


    void matWriter::write_uniform_grids(const grid &g, const char* varname) {
	if (_pmat == NULL) {open();}
	if (g._data.empty()) {
	    std::cout << "No valid data in the grid.\n";
	    return;
	}

	mxArray *matg = mxCreateDoubleMatrix(g._nv, g._dim, mxREAL);
	double *ptrg = mxGetPr(matg);
	for (size_t row = 0; row < g._nv; ++row) {
	    for (int col = 0; col < g._dim; ++col) {
		ptrg[row + g._nv * col] = g._data[row][col];
	    }
	}

	if(matPutVariable(_pmat, varname, matg) != 0) {
	    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
	}
	mxDestroyArray(matg);
    }//matWriter::write_uniform_grids


    void matWriter::write_discrete_controller(const DSolver &dsol,
					      const char* vleastctlr, const char* voptctlr) {

	if (_pmat == NULL) {open();}
	
	mxArray *matopt = mxCreateDoubleMatrix(dsol._ts->_nx, 3, mxREAL);
	double *popt = mxGetPr(matopt);
	mxArray *matleast = mxCreateLogicalMatrix(dsol._ts->_nx, dsol._ts->_nu+1);
	mxLogical *pleast = mxGetLogicals(matleast);

	for (size_t row = 0; row < dsol._ts->_nx; ++row) {

	    /* write optimal controller [winset, optctlr, value] */
	    popt[row] = dsol._win[row];
	    popt[row + dsol._ts->_nx] = dsol._optctlr[row];
	    popt[row + dsol._ts->_nx * 2] = dsol._value[row];

	    /* write least restrictive controller */
	    pleast[row] = dsol._win[row];
	    for (size_t col = 0; col < dsol._ts->_nu; ++col) {
		pleast[row + dsol._ts->_nx * (col + 1)] = dsol._leastctlr[col + row * dsol._ts->_nu];
	    }
	}
    
	if (matPutVariable(_pmat, voptctlr, matopt) ||
	    matPutVariable(_pmat, vleastctlr, matleast) != 0) {
	    std::cout << "Error writing symbolic controllers to matlab variables.\n";
	}

	mxDestroyArray(matopt);
	mxDestroyArray(matleast);
    }//matWriter::write_discrete_controller


    void matReader::read_transitions(fts &trans)
    {
	/* get number of states and inputs */
	mxArray *nvar = matGetVariable(_pmat, "n");
	trans._nx = *mxGetPr(nvar);
	nvar = matGetVariable(_pmat, "n");
	trans._nu = *mxGetPr(nvar);
	

	/* get _idpost, _npost and _ptrpost */
	mxArray *matv = matGetVariable(_pmat, "trans_post");
	if (matv == NULL) {
	    printf("Error reading existing matrix LocalDouble\n");
	    return;
	}
	mwSize T = mxGetM(matv);
	trans._ntrans = T;
	trans._idpost.resize(T);
	double *ptrv = mxGetPr(matv);
	for (mwSize iter = 0; iter < T; ++iter) {
	    trans._idpost[iter] = static_cast<size_t>(ptrv[iter + T]); // extract the 2nd column
	}

	mxArray *matno = matGetVariable(_pmat, "post");
	mxArray *matptr = matGetVariable(_pmat, "postptr");
	if(mxGetM(matno) != mxGetM(matptr)) {
	    printf("Incorrect post numbers and post pointers.\n");
	    return;
	}
	mwSize NM = mxGetM(matno);
	trans._npost.resize(NM);
	trans._ptrpost.resize(NM);
	double *pno = mxGetPr(matno);
	double *pptr = mxGetPr(matptr);
	for (mwSize iter = 0; iter < NM; ++iter) {
	    trans._npost[iter] = static_cast<size_t>(pno[iter]);
	    trans._ptrpost[iter] = static_cast<size_t>(pptr[iter]);
	}

	/* get _idpre, _npre and _ptrpre */
	matv = matGetVariable(_pmat, "trans_pre");
	if (matv == NULL || mxGetM(matv)!= T) {
	    printf("Error reading existing matrix LocalDouble\n");
	    return;
	}
	trans._idpost.resize(T);
	ptrv = mxGetPr(matv);
	for (mwSize iter = 0; iter < T; ++iter) {
	    trans._idpre[iter] = static_cast<size_t>(ptrv[iter + T]); // extract the 2nd column
	}

	matno = matGetVariable(_pmat, "pre");
	matptr = matGetVariable(_pmat, "preptr");
	if(mxGetM(matno) != mxGetM(matptr) || mxGetM(matno) != NM) {
	    printf("Incorrect pre numbers and pre pointers.\n");
	    return;
	}
	trans._npre.resize(NM);
	trans._ptrpost.resize(NM);
	pno = mxGetPr(matno);
	pptr = mxGetPr(matptr);
	for (mwSize iter = 0; iter < NM; ++iter) {
	    trans._npre[iter] = static_cast<size_t>(pno[iter]);
	    trans._ptrpre[iter] = static_cast<size_t>(pptr[iter]);
	}	

	mxDestroyArray(matv);
	mxDestroyArray(nvar);
	mxDestroyArray(matno);
	mxDestroyArray(matptr);
    }//matReader::read_transitions
    
}//namespace rocs
