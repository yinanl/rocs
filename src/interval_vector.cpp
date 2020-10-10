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
    	assert(val.size() == (size_t)rhs.getdim());
	
    	ivec r(val.size());
    	for (size_t i = 0; i < val.size(); ++i)
    	    r[i] = val[i] - rhs[i];
	
    	return r;
    }
    ivec operator-(const ivec &lhs, const std::vector<double> &val) {
    	assert(val.size() == (size_t)lhs.getdim());
	
    	ivec r(val.size());
    	for (size_t i = 0; i < val.size(); ++i)
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
    	assert(val.size() == (size_t)rhs.getdim());
	
    	ivec r(val.size());
    	for (size_t i = 0; i < val.size(); ++i)
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
	if (x.isempty())
	    out << "The interval vector is empty.\n";
	else {
	    for (int i = 0; i < x._dim - 1; ++i) 
		out << x._itvls[i] << "x";

	    out << x._itvls[x._dim - 1];
	}

	return out;
    }


} // namespace rocs
