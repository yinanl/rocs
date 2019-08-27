/**
 *  interval_vector.h
 *  
 *  An interval vector class with a given dimension size.
 *
 *  Created by Yinan Li on Aug. 08, 2016.
 *  
 *  Hybrid Systems Group, University of Waterloo.
 */




#ifndef _interval_vector_h_
#define _interval_vector_h_

#include <vector>

#include <armadillo>
#include "interval.h"


namespace rocs {

    class ivec
    {
    private:
	interval* _itvls; /**< A pointer to an interval. */
	int _dim; /**< The dimension of the interval vector. */

    public:

	/**
	 * \brief Constructors.
	 */
	ivec(): _itvls(NULL), _dim(0) {};
	ivec(int n) { _dim = n; _itvls = new interval[n];}
	/**
	 * \brief A copy constructor: y = ivec(x).
	 */
	ivec(const ivec &x) : _dim(x._dim) {
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
	 * \brief A copy assignment: y = x.
	 */
	ivec& operator=(const ivec &x) {
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
		} // end if
		
	    } // end if

	    return (*this);
	}
	/**
	 * \brief A destructor.
	 */
	~ivec() { delete[] _itvls;}

	/**
	 * \brief A [] operator accessing elements.
	 */
	interval& operator[](const int i) const { return _itvls[i];}
	
	/**
	 * \brief An empty test.
	 *
	 * x = \empty if exists an empty dimension.
	 */
	bool isempty() const;
  
	/**
	 * \brief Get the lower bound of an interval vector.
	 *
	 * (a1.inf,a2.inf,...,an.inf]), vec = [a1] x [a2] x... [an]
	 */
	std::vector<double> getinf() const;
	
	/**
	 * \brief Get the upper bound of an interval vector.
	 *
	 * (a1.sup, a2.sup,..., an.sup), vec = [a1] x [a2] x... [an]
	 */
	std::vector<double> getsup() const;
	
	/**
	 * \brief Get the dimension of an interval vector.
	 */
	int getdim() const { return _dim;}
	/**
	 * \brief Set the value to a specified dimension of the interval vector.
	 */
	void setval(int axis, const interval &x) { _itvls[axis] = x;}
	
	/**
	 * \brief Get the width of each dimension.
	 *
	 * (w[a1],w[a2],...w[an]), vec = [a1] x [a2] x... [an]
	 */
	std::vector<double> width() const;
	
	/**
	 * \brief Get the radius of each dimension.
	 *
	 * (w[a1]/2,w[a2]/2,...w[an]/2), vec = [a1] x [a2] x... [an]
	 */
	std::vector<double> radius() const;
	
	/**
	 * \brief Get the center point of an interval vector.
	 *
	 * (mid[a1],mid[a2],...mid[an]), vec = [a1] x [a2] x... [an]
	 */
	std::vector<double> mid() const;
	
	/**
	 * \brief Get the maximum width of all dimensions.
	 *
	 * max([a1],[a2],...[an]), vec = [a1] x [a2] x... [an]
	 */
	double maxwidth() const;
	
	/**
	 * \brief Determine the dimension which has the maximum width.
	 */
	int maxdim() const;

	/**
	 * \brief A test to see if an interval vector is outside the other.
	 *
	 * x.isout(y), 1 if x, y are disjoint
	 */
	bool isout(const ivec &y) const;
	bool isout(const std::vector<double> &y) const;
	
	/**
	 * \brief A test to see if an interval vector is inside the other.
	 *
	 * x.isin(y), 1 if y is contained in x
	 */
	bool isin(const ivec &y) const;
	bool isin(const std::vector<double> &y) const;
  
	/**
	 * Intersection of two interval vectors.
	 */
	friend ivec intersect(const ivec&, const ivec&);
	
	/**
	 * Interval hull of two interval vectors.
	 */
	friend ivec hull(const ivec&, const ivec&);

	/**
	 * \brief Bisect x along a given dimension to lowerhalf and upperhalf.
	 */
	friend ivec lowerhalf(const ivec&, const int);
	friend ivec upperhalf(const ivec&, const int);

	/**
	 * Operation overloads between interval vectors.
	 */
	friend bool operator==(const ivec&, const ivec&);
	friend bool operator!=(const ivec&, const ivec&);

	friend ivec operator-(const ivec&, const ivec&);
	friend ivec operator-(const double, const ivec&);
	friend ivec operator-(const ivec&, const double);
	friend ivec operator-(const std::vector<double>&, const ivec&);
	friend ivec operator-(const ivec&, const std::vector<double>&);

	friend ivec operator+(const ivec&, const ivec&);
	friend ivec operator+(const double, const ivec&);
	friend ivec operator+(const ivec&, const double);
	friend ivec operator+(const std::vector<double>&, const ivec&);
	friend ivec operator+(const ivec&, const std::vector<double>&);

	friend ivec operator*(const double, const ivec&);
	friend ivec operator*(const ivec&, const double);
  
	/**
	 * \brief A linear affine operation on an interval vector.
	 * y = A[x] + b (b=0 -> y = A[x])
	 * y = A * xc + |A|*|xr|, x = xc + [xr]
	 */
	friend ivec linmap(const arma::mat&, const arma::vec&, const ivec&);

  
	/**
	 * I/O 
	 */
	friend std::ostream& operator<<(std::ostream&, const ivec&);

    };


    /**************** Inline functions ***************/
    inline bool ivec::isempty() const {
	if (_itvls == NULL || _dim == 0)
	    return true;
    
	for (int i = 0; i < _dim; ++i) {
	    if (_itvls[i].isempty())
		return true;
	}

	return false;
    }

    inline std::vector<double> ivec::getinf() const {
	std::vector<double> inf(_dim);
	for (int i = 0; i < _dim; ++i)
	    inf[i] = _itvls[i].getinf();

	return inf;
    }
    inline std::vector<double> ivec::getsup() const {
	std::vector<double> sup(_dim);
	for (int i = 0; i < _dim; ++i) 
	    sup[i] = _itvls[i].getsup();

	return sup;
    }

    inline std::vector<double> ivec::width() const {
	std::vector<double> width(_dim);
	for (int i = 0; i < _dim; ++i)
	    width[i] = _itvls[i].width();

	return width;
    }

    inline std::vector<double> ivec::radius() const {
	std::vector<double> radius(_dim);
	for (int i = 0; i < _dim; ++i)
	    radius[i] = _itvls[i].width() / 2;

	return radius;
    }

    inline std::vector<double> ivec::mid() const {
	std::vector<double> mid(_dim);
	for (int i = 0; i < _dim; ++i)
	    mid[i] = _itvls[i].mid();

	return mid;
    }

    inline double ivec::maxwidth() const {
	double wid = 0;
	for (int i = 0; i < _dim; ++i) {
	    if (wid < _itvls[i].width()) {
		wid = _itvls[i].width();
	    }
	}

	return wid;
    }

    inline int ivec::maxdim() const {
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

    inline bool ivec::isout(const ivec &y) const {
	for (int i = 0; i < _dim; ++i) {
	    if (_itvls[i].isout(y[i]))
		return true;
	}

	return false;
    }

    inline bool ivec::isout(const std::vector<double> &y) const {
	for (int i = 0; i < _dim; ++i) {
	    if (_itvls[i].isout(y[i]))
		return true;
	}

	return false;
    }

    inline bool ivec::isin(const ivec &y) const {
	for (int i = 0; i < _dim; ++i) {
	    if (! _itvls[i].isin(y[i]))
		return false;
	}

	return true;
    }
    
    inline bool ivec::isin(const std::vector<double> &y) const {
	for (int i = 0; i < _dim; ++i) {
	    if (! _itvls[i].isin(y[i]))
		return false;
	}

	return true;
    }
    
    inline ivec intersect(const ivec &x, const ivec &y) {
    	assert(x._dim == y._dim);
    	ivec r(x._dim);
    	for (int i = 0; i < x._dim; ++i) {
    	    r.setval(i, intersect(x[i], y[i]));
    	}
    	return r;
    }

    inline ivec hull(const ivec &x, const ivec &y) {
    	assert(x._dim == y._dim);
    	ivec r(x._dim);
    	for (int i = 0; i < x._dim; ++i) {
    	    r.setval(i, hull(x[i], y[i]));
    	}
    	return r;
    }
    

    // inline bool operator==(const ivec &lhs, const ivec &rhs) {
    // 	if (lhs.getdim() == rhs.getdim()) {
    // 	    for (int i = 0; i < lhs.getdim(); ++i) {
    // 		if (lhs[i] != rhs[i])
    // 		    return false;
    // 	    }
    // 	    return true;
    // 	} else {
    // 	    return false;
    // 	}
    // }

    // inline bool operator!=(const ivec &lhs, const ivec &rhs) {
    // 	return !(lhs == rhs);
    // }

    // inline ivec operator-(const ivec &lhs, const ivec &rhs) {
    // 	int dim = rhs.getdim();
    // 	if (lhs.getdim() == dim) {
    // 	    ivec r(dim);
    // 	    for (int i = 0; i < dim; ++i)
    // 		r[i] = lhs[i] - rhs[i];
    // 	    return r;
    // 	} else {
    // 	    ivec r;
    // 	    return r;
    // 	}
    // }

    // inline ivec operator-(const double val, const ivec &rhs) {
    // 	ivec r(rhs.getdim());
    // 	for (int i = 0; i < rhs.getdim(); ++i)
    // 	    r[i] = val - rhs[i];
    // 	return r;
    // }

    // inline ivec operator-(const ivec &lhs, const double val) {
    // 	ivec r(lhs.getdim());
    // 	for (int i = 0; i < lhs.getdim(); ++i)
    // 	    r[i] = lhs[i] - val;
    // 	return r;
    // }
    
    // inline ivec operator-(const std::vector<double> &val, const ivec &rhs) {
    // 	assert(val.size() == rhs.getdim());
    // 	ivec r(val.size());
    // 	for (int i = 0; i < val.size(); ++i)
    // 	    r[i] = val[i] - rhs[i];
    // 	return r;
    // }
    
    // inline ivec operator-(const ivec &lhs, const std::vector<double> &val) {
    // 	assert(val.size() == lhs.getdim());
    // 	ivec r(val.size());
    // 	for (int i = 0; i < val.size(); ++i)
    // 	    r[i] = lhs[i] - val[i];
    // 	return r;
    // }

    // inline ivec operator+(const ivec &lhs, const ivec &rhs) {
    // 	int dim = rhs.getdim();
    // 	if (lhs.getdim() == dim) {
    // 	    ivec r(dim);
    // 	    for (int i = 0; i < dim; ++i)
    // 		r[i] = lhs[i] + rhs[i];
	
    // 	    return r;
	
    // 	} else {
    // 	    ivec r;
    // 	    return r;
    // 	}
    // }

    // inline ivec operator+(const double val, const ivec &rhs) {
    // 	ivec r(rhs.getdim());
    // 	for (int i = 0; i < rhs.getdim(); ++i)
    // 	    r[i] = val + rhs[i];
	
    // 	return r;
    // }

    // inline ivec operator+(const ivec &lhs, const double val) {
    // 	ivec r(lhs.getdim());
    // 	for (int i = 0; i < lhs.getdim(); ++i)
    // 	    r[i] = lhs[i] + val;
	
    // 	return r;
    // }
    
    // inline ivec operator+(const std::vector<double> &val, const ivec &rhs) {
    // 	assert(val.size() == rhs.getdim());
    // 	ivec r(val.size());
    // 	for (int i = 0; i < val.size(); ++i)
    // 	    r[i] = val[i] + rhs[i];
	
    // 	return r;
    // }
    
    // inline ivec operator+(const ivec &lhs, const std::vector<double> &val) {
    // 	return val + lhs;
    // }

    // inline ivec operator*(const double val, const ivec &rhs) {
    // 	ivec r(rhs.getdim());
    // 	for (int i = 0; i < rhs.getdim(); ++i)
    // 	    r[i] = val * rhs[i];
	
    // 	return r;
    // }

    // inline ivec operator*(const ivec &lhs, const double val) {
    // 	ivec r(lhs.getdim());
    // 	for (int i = 0; i < lhs.getdim(); ++i)
    // 	    r[i] = lhs[i] * val;
	
    // 	return r;
    // }
    
    
    // inline ivec lowerhalf(const ivec &self, const int axis) {
    // 	if (self.isempty()) {
    // 	    return self;
    // 	}
    // 	else {
    // 	    ivec left = self;
    // 	    left.setval(axis, lowerhalf(self[axis]));
    // 	    return left;
    // 	} 
    // }

    // inline ivec upperhalf(const ivec &self, const int axis) {
    // 	if (self.isempty()) {
    // 	    return self;
    // 	}
    // 	else {
    // 	    ivec right = self;
    // 	    right.setval(axis, upperhalf(self[axis]));
    // 	    return right;
    // 	}
    // }

    // inline ivec linmap(const arma::mat &A, const arma::vec &b, const ivec &x) {
    // 	ivec y(x.getdim());
    // 	arma::vec xc(x.mid());
    // 	arma::vec xr(x.radius());

    // 	arma::vec yl = A * xc - arma::abs(A) * arma::abs(xr) + b;
    // 	arma::vec yu = A * xc + arma::abs(A) * arma::abs(xr) + b;

    // 	for (int i = 0; i < x.getdim(); ++ i) 
    // 	    y[i] = interval(yl[i], yu[i]);
    
    // 	return y;
    // }

} // namespace rocs

#endif
