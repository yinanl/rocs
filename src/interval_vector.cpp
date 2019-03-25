/**
 *  interval_vector.cpp
 *  
 *  An interval vector class with a given dimension size
 *
 *  Created by Yinan Li on Aug. 08, 2016.
 *  
 *  Hybrid Systems Group, University of Waterloo.
 */


#include "interval_vector.h"


namespace rocs {

    /**
     * A copy constructor: y = ivec(x)
     */
    ivec::ivec(const ivec &x) {

	_dim = x._dim;

	if (x._itvls == NULL) {

	    _itvls = NULL;
	}
	else {

	    _itvls = new interval[_dim];

	    for (int i = 0; i < _dim; i++)
		_itvls[i] = x._itvls[i];
	}
    }


    /**
     * A copy assignment: y = x 
     */
    ivec& ivec::operator=(const ivec &x) {

	if (this != &x) {
            if (_dim != x._dim)
                _dim = x._dim;
            
            if (x._itvls == NULL) {
                delete[] _itvls;
                _itvls = NULL;
                
            } else {
                interval* temp = new interval[_dim];
                delete[] _itvls;  // delete memory first to avoid memory leaks.
                _itvls = temp;
                for (int i = 0; i < _dim; i++)
                    _itvls[i] = x._itvls[i];
            }
        }
        
        return (*this);
    }


    /**
     * Basic properties.
     */
    bool ivec::isempty() const {

	if (_itvls == NULL || _dim == 0)
	    return true;
    
	for (int i = 0; i < _dim; ++i) {

	    if (_itvls[i].isempty()) {

		return true;
	    }
	}

	return false;
    }
    
    double ivec::maxwidth() const {

	double wid = 0;
    
	for (int i = 0; i < _dim; ++i) {

	    if (wid < _itvls[i].width()) {

		wid = _itvls[i].width();
	    }
	}

	return wid;
    }
    
    std::vector<double> ivec::width() const {

	std::vector<double> width(_dim);
    
	for (int i = 0; i < _dim; ++i) {

	    width[i] = _itvls[i].width();
	}

	return width;
    }
    
    std::vector<double> ivec::radius() const {

	std::vector<double> radius(_dim);
    
	for (int i = 0; i < _dim; ++i) {

	    radius[i] = _itvls[i].width() / 2;
	}

	return radius;
    }
    
    std::vector<double> ivec::mid() const {

	std::vector<double> mid(_dim);
    
	for (int i = 0; i < _dim; ++i) {

	    mid[i] = _itvls[i].mid();
	}

	return mid;
    }
    
    int ivec::maxdim() const {

	int maxi = 0;
	double wid = 0;
    
	for (int i = 0; i < _dim; ++i) {

	    if (wid < _itvls[i].width()) {

		wid = _itvls[i].width();
		maxi = i;
	    }
	}

	return maxi;
    }


    /**
     * Getters and setters.
     */
    std::vector<double> ivec::getinf() const {

	std::vector<double> inf(_dim);
    
	for (int i = 0; i < _dim; ++i) {

	    inf[i] = _itvls[i].getinf();
	}

	return inf;
    }
    
    std::vector<double> ivec::getsup() const {

	std::vector<double> sup(_dim);
    
	for (int i = 0; i < _dim; ++i) {

	    sup[i] = _itvls[i].getsup();
	}

	return sup;
    }


    /**
     * Intersection tests.
     */
    bool ivec::isout(const ivec &y) const {

	for (int i = 0; i < _dim; ++i) {

	    if (_itvls[i].isout(y[i])) {

		return true;
	    }
	}

	return false;
    }

    bool ivec::isout(const std::vector<double> &y) const {

	for (int i = 0; i < _dim; ++i) {

	    if (_itvls[i].isout(y[i])) {

		return true;
	    }
	}

	return false;
    }

    
    bool ivec::isin(const ivec &y) const {

	for (int i = 0; i < _dim; ++i) {

	    if (! _itvls[i].isin(y[i])) {

		return false;
	    }
	}

	return true;
    }

    bool ivec::isin(const std::vector<double> &y) const {

	for (int i = 0; i < _dim; ++i) {

	    if (! _itvls[i].isin(y[i])) {

		return false;
	    }
	}

	return true;
    }


    /**
     * Set operations. 
     */
    ivec intersect(const ivec &x, const ivec &y) {

	assert(x._dim == y._dim);

	ivec r(x._dim);
    
	for (int i = 0; i < x._dim; ++i) {

	    r.setval(i, intersect(x[i], y[i]));
	}

	return r;
    }

    ivec hull(const ivec &x, const ivec &y) {

	assert(x._dim == y._dim);

	ivec r(x._dim);
    
	for (int i = 0; i < x._dim; ++i) {

	    r.setval(i, hull(x[i], y[i]));
	}

	return r;
    }


    /**
     * Boolean operation overloads.
     */
    bool operator==(const ivec &lhs, const ivec &rhs) {

	if (lhs.getdim() == rhs.getdim()) {

	    for (int i = 0; i < lhs.getdim(); ++i) {

		if (lhs[i] != rhs[i])
		    return false;
	    }

	    return true;
	} else {

	    return false;
	}
    }

    bool operator!=(const ivec &lhs, const ivec &rhs) {

	return !(lhs == rhs);
    }

    ivec operator-(const ivec &lhs, const ivec &rhs) {

	int dim = rhs.getdim();
    
	if (lhs.getdim() == dim) {

	    ivec r(dim);
	    for (int i = 0; i < dim; ++i)
		r[i] = lhs[i] - rhs[i];

	    return r;
	
	} else {

	    ivec r;
	    return r;
	}
    }

    ivec operator-(const double val, const ivec &rhs) {
    
	ivec r(rhs.getdim());
	for (int i = 0; i < rhs.getdim(); ++i)
	    r[i] = val - rhs[i];
	
	return r;
    }

    ivec operator-(const ivec &lhs, const double val) {

	ivec r(lhs.getdim());
	for (int i = 0; i < lhs.getdim(); ++i)
	    r[i] = lhs[i] - val;
	
	return r;
    }
    ivec operator-(const std::vector<double> &val, const ivec &rhs) {
	assert(val.size() == rhs.getdim());
	
	ivec r(val.size());
	for (int i = 0; i < val.size(); ++i)
	    r[i] = val[i] - rhs[i];
	
	return r;
    }
    ivec operator-(const ivec &lhs, const std::vector<double> &val) {
	assert(val.size() == lhs.getdim());
	
	ivec r(val.size());
	for (int i = 0; i < val.size(); ++i)
	    r[i] = lhs[i] - val[i];
	
	return r;
    }
    

    ivec operator+(const ivec &lhs, const ivec &rhs) {

	int dim = rhs.getdim();
    
	if (lhs.getdim() == dim) {

	    ivec r(dim);
	    for (int i = 0; i < dim; ++i)
		r[i] = lhs[i] + rhs[i];
	
	    return r;
	
	} else {

	    ivec r;
	    return r;
	}
    }

    ivec operator+(const double val, const ivec &rhs) {
    
	ivec r(rhs.getdim());
	for (int i = 0; i < rhs.getdim(); ++i)
	    r[i] = val + rhs[i];
	
	return r;
    }

    ivec operator+(const ivec &lhs, const double val) {

	ivec r(lhs.getdim());
	for (int i = 0; i < lhs.getdim(); ++i)
	    r[i] = lhs[i] + val;
	
	return r;
    }
    ivec operator+(const std::vector<double> &val, const ivec &rhs) {
	assert(val.size() == rhs.getdim());
	
	ivec r(val.size());
	for (int i = 0; i < val.size(); ++i)
	    r[i] = val[i] + rhs[i];
	
	return r;
    }
    ivec operator+(const ivec &lhs, const std::vector<double> &val) {
	return val + lhs;
    }


    ivec operator*(const double val, const ivec &rhs) {
    
	ivec r(rhs.getdim());
	for (int i = 0; i < rhs.getdim(); ++i)
	    r[i] = val * rhs[i];
	
	return r;
    }

    ivec operator*(const ivec &lhs, const double val) {

	ivec r(lhs.getdim());
	for (int i = 0; i < lhs.getdim(); ++i)
	    r[i] = lhs[i] * val;
	
	return r;
    }


    /**
     * Bisection.
     */
    ivec lowerhalf(const ivec &self, const int axis) {

	if (self.isempty()) {

	    return self;
	}
	else {

	    ivec left = self;
	    left.setval(axis, lowerhalf(self[axis]));
	
	    return left;
	} 
    }

    ivec upperhalf(const ivec &self, const int axis) {

	if (self.isempty()) {

	    return self;
	}
	else {

	    ivec right = self;
	    right.setval(axis, upperhalf(self[axis]));
	
	    return right;
	} 
    }


    /**
     * A linear operation. 
     */
    ivec linmap(const arma::mat &A, const arma::vec &b, const ivec &x) {

	ivec y(x.getdim());

	arma::vec xc(x.mid());
	arma::vec xr(x.radius());

	arma::vec yl = A * xc - arma::abs(A) * arma::abs(xr) + b;
	arma::vec yu = A * xc + arma::abs(A) * arma::abs(xr) + b;

	for (int i = 0; i < x.getdim(); ++ i) 
	    y[i] = interval(yl[i], yu[i]);
    
	return y;
    }

    
    std::ostream& operator<<(std::ostream &out, const ivec &x) {

	for (int i = 0; i < x._dim - 1; ++i) {

	    out << x._itvls[i] << "x";
	}

	out << x._itvls[x._dim - 1];

	return out;
    }


} // namespace rocs
