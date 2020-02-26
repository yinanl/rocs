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


namespace rocs {

    void grid::gridding() {
	_data.resize(_nv, std::vector<double>(_dim));
	griddingHelper(_data, _valmin, _gw, _size, _nv);
    }

    void grid::gridding(const double eta[], const double lb[], const double ub[]) {
	init(eta, lb, ub);
	gridding();
    }

    void grid::gridding(const double eta[], const ivec &x) {
	init(eta, x);
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
	
	return sub;
    }

    // void grid::griddingHelper(std::vector<std::vector<double> > &data,
    // 			      const std::vector<double> &xmin,
    // 			      const std::vector<double> &gw,
    // 			      const std::vector<size_t> number,
    // 			      const size_t nv) {

    // 	size_t a, b, r, ik;

    // 	for (size_t i = 0; i < nv; ++i) {
    // 	    r = i;
    // 	    a = number[0];
    // 	    b = 1;

    // 	    for (int k = 0; k < _dim; ++k) {
    // 		ik = (r % a) / b;
    // 		data[i][k] = xmin[k] + ik * gw[k];
	    
    // 		if (k < _dim - 1) {
    // 		    r -= ik * b;
    // 		    a *= number[k+1];
    // 		    b *= number[k];
    // 		} //end if
    // 	    } //end for k
    // 	} //end for i
    // }
    void grid::griddingHelper(std::vector<std::vector<double> > &data,
    			      const std::vector<double> &xmin,
    			      const std::vector<double> &gw,
    			      const std::vector<size_t> number,
    			      const size_t nv) {

    	size_t r, ik;
    	for (size_t i = 0; i < nv; ++i) {
    	    r = i;
    	    for (int k = 0; k < _dim; ++k) {
    		ik = r % number[k];
    		data[i][k] = xmin[k] + ik * gw[k];
    		r = r / number[k];
    	    } //end for k
    	} //end for i
    }


    // size_t grid::val_to_id(std::vector<double> val) {

    // 	assert(val.size() == _dim);

    // 	size_t id = 0;
    // 	size_t a = 1;
    
    // 	for (int k = 0; k < _dim; ++k) {

    // 	    id += a * round((val[k] - _valmin[k]) / _gw[k]);

    // 	    a *= _size[k];
    // 	}

    // 	return id;
    // }
    size_t grid::val_to_id(std::vector<double> val) {
    	assert(val.size() == _dim);
	if (_bds.isout(val)) {
	    throw std::runtime_error("rocs::grid:val_to_id: the input value is not in the state space.\n");
	}
    	size_t id = 0;
    	for (int k = _dim - 1; k > 0; --k) {
    	    id = _size[k-1] * (id + round((val[k]-_valmin[k])/_gw[k]));
    	}
    	id += round((val[0]-_valmin[0])/_gw[0]);
    	return id;
    }

    void grid::id_to_val(std::vector<double> &val, size_t id) const {
	if (id > _nv) {
	    throw std::runtime_error("rocs::grid:id_to_val: input id is out of range.\n");
	}
	size_t ik;
	for (int k = 0; k < _dim; ++k) {
	    ik = id % _size[k];
	    val[k] = _valmin[k] + ik * _gw[k];
	    id = id / _size[k];
	}
    }


    std::vector<size_t> grid::subset(ivec &box, bool boxin, bool strictin) {
	assert(_bds.getdim() == _dim);
	std::vector<size_t> ss;

	if (_bds.isout(box))  /* if box is out of range, return empty */
	    return ss;

	ivec x;
	if (_bds.isin(box)) {  /* if box is not fully inside, check boxin */
	    x = box;
	} else {
	    if (boxin)
		return ss;
	    else 
		x = intersect(box, _bds);  /* take the intersection */
	}
    
	/* box and _bds have intersections */
	std::vector<int> il(_dim), iu(_dim);
	std::vector<double> xl(_dim), xu(_dim), xmin(_dim), xmax(_dim);
	x.getinf(xl); x.getsup(xu);
	_bds.getinf(xmin); _bds.getsup(xmax);
	// /********** logging **********/
	// std::cout << "Index range of each dimension:\n";
	// /********** logging **********/
	if (strictin) { /* only collect those fully inside box area */
	    for (int k = 0; k < _dim; ++k) {
		/* compute the index of the lower bound */
		if (fabs(xl[k]-xmin[k]) < EPSIVAL) {
		    /* case 1: box.inf == _bds.inf */
		    il[k] = 0;
		} else if (fabs(fmod(xl[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {
		    /* case 2: box.inf is on the boundary of a grid */
		    il[k] = ceil((xl[k] - _valmin[k]) / _gw[k]);
		} else {/* case 3: rest of the cases */
		    il[k] = ceil((xl[k] - _valmin[k]) / _gw[k] + 0.5);
		}
		// /********** logging **********/
		// std::cout << il[k] << ", ";
		// /********** logging **********/

		/* compute the index of the upper bound */
		if (fabs(xu[k] - xmax[k]) < EPSIVAL) {
		    iu[k] = _size[k] - 1;
		} else if (fabs(fmod(xu[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {
		    iu[k] = floor((xu[k] - _valmin[k]) / _gw[k]);
		} else {
		    iu[k] = floor((xu[k] - _valmin[k]) / _gw[k] - 0.5);
		}
		// /********** logging **********/
		// std::cout << iu[k] << '\n';
		// /********** logging **********/
	    }
	} else { /* intersected box area */
	    for (int k = 0; k < _dim; ++k) {
		if (fabs(fmod(xl[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {
		    il[k] = ceil((xl[k] - _valmin[k]) / _gw[k]);
		} else {
		    il[k] = round((xl[k] - _valmin[k]) / _gw[k]);
		}
		if (fabs(fmod(xu[k]-_valmin[k], _gw[k]) - _gw[k]/2.) < EPSIVAL) {
		    iu[k] = floor((xu[k] - _valmin[k]) / _gw[k]);
		} else {
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
	size_t n = 1;
	std::vector<size_t> range(_dim);
	for (int k = 0; k < _dim; ++k) {
	    if ( (iu[k] - il[k]) < 0) {
		return ss;
	    } else {
		range[k] = iu[k] - il[k] + 1;
		n *= range[k];
	    }
	}

	/* generate indices using the index range in each dimension */
	ss.resize(n);
	ss[0] = 0;
	for (int k = 0; k < _dim; ++k)
	    ss[0] += _base[k] * il[k];
	size_t r, id;
	for (size_t l = 1; l < n; ++l) {
	    r = l; id = ss[0];
	    for (int k = 0; k < _dim; ++k) {
		// std::cout << "r=" << r << ", range=" << range[k] << ", base=" << _base[k] <<'\n';
		id += (r % range[k])*_base[k];
		r /= range[k];
	    }
	    ss[l] = id;
	}
	// size_t len = 1;
	// for (int k = 0; k < _dim; ++k) {
	//     for (size_t j = 1; j < range[k]; ++j) {
	// 	for (size_t h = 0; h < len; ++h) {
	// 	    ss[h + j*len] = ss[h] + _base[k] * j;
	// 	}
	//     }
	//     len *= range[k];
	// }
    
	return ss;
    }


    std::vector<size_t> grid::neighbours(size_t id) {
	std::vector<size_t> nbs;

	size_t a, b, ik, i;
	i = id;
	a = _size[0];
	b = 1;
	
	for (int k = 0; k < _dim; ++k) {
	    ik = (i % a) / b;

	    if(ik > 0)
		nbs.push_back(id-b);
	    if(ik < _size[k])
		nbs.push_back(id+b);
	    
	    if (k < _dim - 1) {
		i -= ik * b;
		a *= _size[k+1];
		b *= _size[k];
	    } //end if
	} //end for k
	
	return nbs;
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

} // namespace rocs
