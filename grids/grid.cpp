/**
 *  grid.cpp
 *
 *  Source file for grid class
 *
 *  Created by Yinan Li on Jan. 21, 2017.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include "grid.h"


void grid::init(const int n, const double eta[], const double lb[], const double ub[]) {
    _dim = n;
    _gw.assign(eta, eta + n);

    ivec cp(n);
    
    _nv = 1;
    _size.resize(n);
    _valmin.resize(n);
    for (int i = 0; i < _dim; ++i) {

	cp[i] = interval(lb[i], ub[i]);

	_size[i] = floor((ub[i] - lb[i]) / _gw[i]) + 1;

	_nv *= _size[i];

	_valmin[i] = ((ub[i]-lb[i])-floor((ub[i]-lb[i])/eta[i])*eta[i]) / 2. + lb[i];
    }

    _bds = cp;  // copy assignment (deep copy)
}

void grid::init(const int n, const double eta[], const ivec &x) {
    _dim = n;
    _gw.assign(eta, eta + n);

    ivec cp(n);
    
    _nv = 1;
    _size.resize(n);
    _valmin.resize(n);

    std::vector<double> lb = x.getinf();
    std::vector<double> ub = x.getsup();
    
    for (int i = 0; i < _dim; ++i) {
	_size[i] = floor((ub[i] - lb[i]) / _gw[i]) + 1;
	_nv *= _size[i];
	cp[i] = interval(lb[i], ub[i]);
	_valmin[i] = ((ub[i]-lb[i])-floor((ub[i]-lb[i])/eta[i])*eta[i]) / 2. + lb[i];
    }

    _bds = cp;  // copy assignment (deep copy)
}


void grid::gridding() {

    _data.resize(_nv, std::vector<double>(_dim));

    griddingHelper(_data, _valmin, _gw, _size, _nv);

    // size_t a, b, r, ik;

    // for (size_t i = 0; i < _nv; ++i) {

    // 	r = i;
    // 	a = _size[0];
    // 	b = 1;

    // 	for (int k = 0; k < _dim; ++k) {

    // 	    ik = (r % a) / b;
    // 	    _data[i][k] = _valmin[k] + ik * _gw[k];
	    
    // 	    if (k < _dim - 1) {
		
    // 		r -= ik * b;
    // 		a *= _size[k+1];
    // 		b *= _size[k];
    // 	    }
    // 	}
	
    // }
    
}

void grid::gridding(const int n, const double eta[],
		    const double lb[], const double ub[]) {
    
    init(n, eta, lb, ub);

    gridding();
}

void grid::gridding(const int n, const double eta[], const ivec &x) {

    init(n, eta, x);

    gridding();
}

std::vector< std::vector<double> >
grid::subgridding(std::vector<double> &xc, std::vector<double> &rp) {

    int nv = 1;
    std::vector<size_t> number(_dim);
    std::vector<double> subgw(_dim);
    std::vector<double> xmin(_dim);
    for (int k = 0; k < _dim; ++k) {
	number[k] = ceil(1.0 / rp[k]);
	subgw[k] = _gw[k] / number[k];
	nv *= number[k];
	// xmin[k] = xc[k] - _gw[k]/2. + rp[k]*_gw[k]/2.;
	xmin[k] = xc[k] - _gw[k]/2. + subgw[k]/2.;
    }
    
    std::vector< std::vector<double> > sub(nv, std::vector<double> (_dim));
    griddingHelper(sub, xmin, subgw, number, nv);
    
    // size_t a, b, r, ik;

    // for (size_t i = 0; i < nv; ++i) {

    // 	r = i;
    // 	a = number[0];
    // 	b = 1;

    // 	for (int k = 0; k < _dim; ++k) {

    // 	    ik = (r % a) / b;
    // 	    sub[i][k] = xmin[k] + ik * rp[k]*_gw[k];
	    
    // 	    if (k < _dim - 1) {
		
    // 		r -= ik * b;
    // 		a *= number[k+1];
    // 		b *= number[k];
    // 	    }
    // 	}
	
    // }

    return sub;
}

void grid::griddingHelper(std::vector<std::vector<double> > &data,
			  const std::vector<double> &xmin,
			  const std::vector<double> &gw,
			  const std::vector<size_t> number,
			  const size_t nv) {

    size_t a, b, r, ik;

    for (size_t i = 0; i < nv; ++i) {
	r = i;
	a = number[0];
	b = 1;

	for (int k = 0; k < _dim; ++k) {
	    ik = (r % a) / b;
	    data[i][k] = xmin[k] + ik * gw[k];
	    
	    if (k < _dim - 1) {
		r -= ik * b;
		a *= number[k+1];
		b *= number[k];
	    } //end if
	} //end for k
    } //end for i
}


size_t grid::val_to_id(std::vector<double> val) {

    assert(val.size() == _dim);

    size_t index = 0;
    size_t a = 1;
    
    for (int k = 0; k < _dim; ++k) {

	index += a * round((val[k] - _valmin[k]) / _gw[k]);

	a *= _size[k];
    }

    return index;
}


std::vector<size_t> grid::subset(ivec &box, bool boxin, bool strictin) {

    assert(_bds.getdim() == _dim);

    std::vector<size_t> ss;

    if (_bds.isout(box))  // if box is out of range, return empty
	return ss;

    ivec x;
    if (_bds.isin(box)) {  // if box is not fully inside, check boxin

	x = box;
    }
    else {

	if (boxin) {
	    
	    return ss;
	}
	else {

	    x = intersect(box, _bds);  // take the intersection
	}
    }

    
    /* box and _bds have intersections */
    std::vector<int> il(_dim);
    std::vector<int> iu(_dim);

    std::vector<double> xl = x.getinf();
    std::vector<double> xu = x.getsup();

    // /********** logging **********/
    // std::cout << "Index range of each dimension:\n";
    // /********** logging **********/
    if (strictin) {  // only collect those fully inside box area

	for (int k = 0; k < _dim; ++k) {

	    /* compute the index of the lower bound */
	    if (fabs(xl[k] - _bds.getinf()[k]) < EPSIVAL) {

		// case 1: box.inf == _bds.inf
		il[k] = 0;
	    }
	    else if (fabs(fmod(xl[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {

		// case 2: box.inf is on the boundary of a grid
		il[k] = ceil((xl[k] - _valmin[k]) / _gw[k]);
	    }
	    else {

		// case 3: rest of the cases
		il[k] = ceil((xl[k] - _valmin[k]) / _gw[k] + 0.5);
	    }
	    // /********** logging **********/
	    // std::cout << il[k] << ", ";
	    // /********** logging **********/
	    

	    /* compute the index of the upper bound */
	    if (fabs(xu[k] - _bds.getsup()[k]) < EPSIVAL) {

		iu[k] = _size[k] - 1;
	    }
	    else if (fabs(fmod(xu[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {

		iu[k] = floor((xu[k] - _valmin[k]) / _gw[k]);
	    }
	    else {

		iu[k] = floor((xu[k] - _valmin[k]) / _gw[k] - 0.5);
	    }
	    // /********** logging **********/
	    // std::cout << iu[k] << '\n';
	    // /********** logging **********/
	}
	
    }
    else {

	for (int k = 0; k < _dim; ++k) {

	    if (fabs(fmod(xl[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {
		
		il[k] = ceil((xl[k] - _valmin[k]) / _gw[k]);
	    }
	    else {
		
		il[k] = round((xl[k] - _valmin[k]) / _gw[k]);
	    }

	    if (fabs(fmod(xu[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {
		
		iu[k] = floor((xu[k] - _valmin[k]) / _gw[k]);
	    }
	    else {
		
		iu[k] = round((xu[k] - _valmin[k]) / _gw[k]);
	    }
	}
    }

    // /********** logging **********/
    // for (int k = 0; k < _dim; ++k) {
    // 	std::cout << "[" << il[k] << " " << iu[k] << "]";
    // 	if (k < _dim - 1)
    // 	    std::cout << 'x';
    // 	else
    // 	    std::cout << '\n';
    // }
    // /********** logging **********/

    size_t n = 1;  // n is the number of subset grids
    std::vector<size_t> base(_dim);
    base[0] = 1;
    std::vector<size_t> range(_dim);
    for (int k = 0; k < _dim; ++k) {

	if ( (iu[k] - il[k]) < 0) {
	    
	    return ss;
	}
	else {
	    
	    range[k] = iu[k] - il[k] + 1;

	    n *= range[k];

	    if (k > 0)
		base[k] = _size[k - 1] * base[k - 1];
	}
    }

    /* generate indices using the index range in each dimension */
    ss.resize(n);
    ss[0] = 0;
    for (int k = 0; k < _dim; ++k) {

	ss[0] += base[k] * il[k];
    }

    size_t len = 1;
    for (int k = 0; k < _dim; ++k) {

	for (size_t j = 1; j < range[k]; ++j) {

	    for (size_t h = 0; h < len; ++h) {

		ss[h + j*len] = ss[h] + base[k] * j;
	    }
	}

	len *= range[k];
    }
    
    return ss;
}


/* 
 * displaying, saving, and logging
 */
void grid::write2txt_data(const char *filename) {

    std::ofstream txtfile;
    txtfile.open(filename, std::ios::out);
    if (!txtfile.is_open()) {

	std::cout << "Unable to open file " << filename << ".";
	return;
    }


    for (size_t row = 0; row < _nv; ++row) {

	txtfile << row << " ";

	for (int col = 0; col < _dim; ++col) {

	    txtfile << _data[row][col] << ' ';
	}

	txtfile << '\n';
    }

    txtfile.close();
}

void grid::write2mat_data(const char *filename, const char *varname) {

    if (_data.empty()) {

	std::cout << "No valid data in the grid.\n";
	return;
    }
    
    MATFile *pmat;
    pmat = matOpen(filename, "u");
    if (pmat == NULL) {

	pmat = matOpen(filename, "w");
	// std::cout << "Error creating file " << filename << '\n';
	// return ;
    }

    mxArray *matg = mxCreateDoubleMatrix(_nv, _dim, mxREAL);
    double *ptrg = mxGetPr(matg);

    for (size_t row = 0; row < _nv; ++row) {

	for (int col = 0; col < _dim; ++col) {

	    ptrg[row + _nv * col] = _data[row][col];
	}
    }

    if(matPutVariable(pmat, varname, matg) != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
    }

    if (EOF == matClose(pmat)) {
        return ;
    }
    
}


/* display functions */
void grid::print_info() {

    std::cout << "Problem scale: " << _dim << '\n';

    std::cout << "Gridding area: " << _bds << '\n';
    
    std::cout << "Grid widths: [ ";
    for (int i = 0; i < _dim; ++i) {

	std::cout << _gw[i] << ' ';
    }
    std::cout << "]\n";

    std::cout << "Dimensional number of intervals: [ ";
    for (int i = 0; i < _dim; ++i) {

	std::cout << _size[i] << ' ';
    }
    std::cout << "]\n";

    std::cout << "Number of grids: " << _nv << '\n';

    std::cout << "Minimum grid center: [ ";
    for (int i = 0; i < _dim; ++i) {

	std::cout << _valmin[i] << ' ';
    }
    std::cout << "]\n";
    
}


void grid::print_data() {

    if (_data.empty()) {

	std::cout << "No valid data in the grid.\n";
	return;
    }

    for (size_t row = 0; row < _nv; ++row) {

	std::cout << row << " ";

	for (int col = 0; col < _dim; ++col) {

	    std::cout << _data[row][col] << ' ';
	}

	std::cout << '\n';
    }
    
}
