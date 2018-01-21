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

/********** member functions ***********/

/*---------------------
   copy constructor: y = ivec(x)
---------------------*/
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


/*---------------------
   assignment: y = x 
---------------------*/
ivec& ivec::operator=(const ivec &x) {

    if (this == &x) {

	return (*this);
    }
    else {

	if (_dim != x._dim) {

	    _dim = x._dim;
	    delete[] _itvls;
	}
	if (x._itvls == NULL) {

	    _itvls = NULL;
	} else {

	    _itvls = new interval[_dim];
	}

	for (int i = 0; i < _dim; i++)
	    _itvls[i] = x._itvls[i];

	return (*this);
    }
}


/*---------------------
   isempty():
   x = \empty if any dimension is empty
---------------------*/
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


/*---------------------
   maxwidth():
   max([a1],[a2],...[an]), vec = [a1] x [a2] x... [an]
---------------------*/
double ivec::maxwidth() const {

    double wid = 0;
    
    for (int i = 0; i < _dim; ++i) {

	if (wid < _itvls[i].width()) {

	    wid = _itvls[i].width();
	}
    }

    return wid;
}


/*---------------------
   width():
   (w[a1],w[a2],...w[an]), vec = [a1] x [a2] x... [an]
---------------------*/
std::vector<double> ivec::width() const {

    std::vector<double> width(_dim);
    
    for (int i = 0; i < _dim; ++i) {

	width[i] = _itvls[i].width();
    }

    return width;
}


/*---------------------
   radius():
   (w[a1]/2,w[a2]/2,...w[an]/2), vec = [a1] x [a2] x... [an]
---------------------*/
std::vector<double> ivec::radius() const {

    std::vector<double> radius(_dim);
    
    for (int i = 0; i < _dim; ++i) {

	radius[i] = _itvls[i].width() / 2;
    }

    return radius;
}


/*---------------------
   mid():
   (mid[a1],mid[a2],...mid[an]), vec = [a1] x [a2] x... [an]
---------------------*/
std::vector<double> ivec::mid() const {

    std::vector<double> mid(_dim);
    
    for (int i = 0; i < _dim; ++i) {

	mid[i] = _itvls[i].mid();
    }

    return mid;
}


/*---------------------
   maxdim():
   Return the dimension with maximal width.
---------------------*/
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


/*---------------------
   getinf():
   (a1.inf,a2.inf,...,an.inf]), vec = [a1] x [a2] x... [an]
---------------------*/
std::vector<double> ivec::getinf() const {

    std::vector<double> inf(_dim);
    
    for (int i = 0; i < _dim; ++i) {

	inf[i] = _itvls[i].getinf();
    }

    return inf;
}


/*---------------------
   getsup():
   (a1.sup, a2.sup,..., an.sup), vec = [a1] x [a2] x... [an]
---------------------*/
std::vector<double> ivec::getsup() const {

    std::vector<double> sup(_dim);
    
    for (int i = 0; i < _dim; ++i) {

	sup[i] = _itvls[i].getsup();
    }

    return sup;
}


/*---------------------
   isout():

   x.isout(y), 1 if x, y are disjoint
---------------------*/
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


/*---------------------
   isin():

   x.isin(y), 1 if y is contained in x
---------------------*/
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




/********** friend functions ***********/

/*
 * Intersection of two interval vectors
 */
ivec intersect(const ivec &x, const ivec &y) {

    assert(x._dim == y._dim);

    ivec r(x._dim);
    
    for (int i = 0; i < x._dim; ++i) {

	r.setval(i, intersect(x[i], y[i]));
    }

    return r;
}


/*
 * interval hull of two interval vectors
 */
ivec hull(const ivec &x, const ivec &y) {

    assert(x._dim == y._dim);

    ivec r(x._dim);
    
    for (int i = 0; i < x._dim; ++i) {

	r.setval(i, hull(x[i], y[i]));
    }

    return r;
}


/*---------------------
   boolean operation overload
---------------------*/
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


/*---------------------
   bisect x along the axis dimension:
   lowerhalf
   upperhalf
---------------------*/
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


/*---------------------
  Linear affine operation on [x]:
  y = A[x] + b (b=0 -> y = A[x])

  y = A * xc + |A|*|xr|, x = xc + [xr]
---------------------*/
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



/*---------------------
   I/O
---------------------*/
std::ostream& operator<<(std::ostream &out, const ivec &x) {

    for (int i = 0; i < x._dim - 1; ++i) {

	out << x._itvls[i] << "x";
    }

    out << x._itvls[x._dim - 1];

    return out;
}


} // namespace rocs
