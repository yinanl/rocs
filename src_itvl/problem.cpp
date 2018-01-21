/**
 * problem.cpp
 *
 * A control problem class.
 *
 * Created by Yinan Li on Jan. 04, 2017.
 *
 * Hybrid Systems Group, University of Waterloo.
 */


#include "problem.h"


namespace rocs {

void CntlProb::write2mat_settings(const char* filename) {

    MATFile *pmat;
    pmat = matOpen(filename, "w");
    if (pmat == NULL) {
	std::cout << "Error creating file" << filename;
	return;
    }

    mxArray *mx = mxCreateDoubleMatrix(_xdim, 2, mxREAL); // X
    double *px = mxGetPr(mx);
    mxArray *mu = mxCreateDoubleMatrix(_vf->_ugrid._nv, _udim, mxREAL); // U
    double *pu = mxGetPr(mu);
    mxArray *mt = mxCreateDoubleScalar(_vf->_tau); // ts

    mxArray *mg = mxCreateDoubleMatrix(_xdim, 2, mxREAL); // G (goal)
    double *pg = mxGetPr(mg);

    for (int i = 0; i < _xdim; ++i) {

	px[i] = _workspace[i].getinf(); // X
	px[i+_xdim] = _workspace[i].getsup();

	pg[i] = _goal[i].getinf(); // G
	pg[i+_xdim] = _goal[i].getsup();
    }

    std::vector< std::vector<double> > udata;
    if (!_vf->_ugrid._data.empty()) {

	udata = _vf->_ugrid._data;
    } else
	assert(false);

    for (int r = 0; r < _vf->_ugrid._nv; ++r) 
	for (int c = 0; c < _udim; ++c)
	    pu[r + c*_vf->_ugrid._nv] = udata[r][c];

    /* write to variables in .mat */
    if (matPutVariable(pmat, "X", mx) ||
	matPutVariable(pmat, "U", mu) ||
	matPutVariable(pmat, "ts", mt)||
	matPutVariable(pmat, "G", mg)!= 0) {
	
	std::cout << "Error outputting variables.\n";
        return;
    }
    
    mxDestroyArray(mx);
    mxDestroyArray(mu);
    mxDestroyArray(mt);
    mxDestroyArray(mg);

    /* output _obs only if it is nonempty */
    if (!_obs.empty()) {
	mwSize dims[3] = {(mwSize)_xdim, 2, _obs.size()};
	mxArray *mo = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // obs
	double *po = mxGetPr(mo);
	for (int r = 0; r < dims[0]; ++r)
	    for (int p = 0; p < dims[2]; ++p) {
		ivec obstacle = _obs[p];
		po[r + p*dims[0]*dims[1]] = obstacle[r].getinf();
		po[r + dims[0] + p*dims[0]*dims[1]] = obstacle[r].getsup();
	    }

	matPutVariable(pmat, "xobs", mo);
	mxDestroyArray(mo);
    }
    
    
    if (EOF == matClose(pmat)) {
        return;
    }
}


/* I/O friend functions */
std::ostream& operator<<(std::ostream &out, const CntlProb &prob) {

    out << '<' << prob._name << '>' << ":\n";
    out << "- state dimension: " << prob._xdim << '\n';
    out << "- input dimension: " << prob._udim << '\n';
  
    out << "- state space: " << '\n';
    out << prob._workspace << '\n';

    out << "- number of controls: " << prob._vf->_ugrid._nv << '\n';

    // out << "- the controller:\n ";
    // prob._ctlr.print_leaves(prob._ctlr._root);
    
  
    return out;
}

} // namespace rocs
