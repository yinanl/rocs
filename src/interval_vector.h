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
#include <initializer_list>
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
		for (int i = 0; i < _dim; ++i)
		    _itvls[i] = x._itvls[i];
	    }
	}
	/**
	 * \brief Constructor by initializer_list.
	 */
	ivec(std::initializer_list<interval> l) {
	    _dim = l.size();
	    _itvls = new interval[_dim];
	    std::initializer_list<interval>::iterator iter;
	    int i = 0;
	    for (iter = l.begin(); iter != l.end(); ++iter){
		_itvls[i] = *iter;
		++i;
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
	void getinf(std::vector<double> &inf) const;
	
	/**
	 * \brief Get the upper bound of an interval vector.
	 *
	 * (a1.sup, a2.sup,..., an.sup), vec = [a1] x [a2] x... [an]
	 */
	std::vector<double> getsup() const;
	void getsup(std::vector<double> &sup) const;
	
	/**
	 * \brief Get the dimension of an interval vector.
	 */
	int getdim() const { return _dim;}
	
	/**
	 * \brief Set the interval in the axis (the first argument) dimension.
	 */
	void setval(int axis, const interval &x) { _itvls[axis] = x;}

	/**
	 * \brief Get the interval in the axis (the first argument) dimension.
	 */
	interval getval(int axis) { return _itvls[axis];}
	
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
	 * \brief Compute the volume of the interval vector.
	 *
	 */
	double volume() const;

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
	 * \brief Add an interval vector to (*this)
	 */
	ivec& operator += (const ivec &x);
	/**
	 * \brief Add a real vector to (*this)
	 */
	ivec& operator += (const std::vector<double> &x);
	/**
	 * \brief Add a real number a (broadcasted to the dim of (*this)) to (*this).
	 */
	ivec& operator += (const double a);

	/**
	 * \brief Subtract an interval vector from (*this)
	 */
	ivec& operator -= (const ivec &x);
	/**
	 * \brief Subtract a real vector from (*this)
	 */
	ivec& operator -= (const std::vector<double> &x);
	/**
	 * \brief Subtract a real number (broadcasted to the dim of (*this)) a from (*this).
	 */
	ivec& operator -= (const double a);
	
  
	/**
	 * \brief Intersection of two interval vectors.
	 */
	friend ivec intersect(const ivec&, const ivec&);
	
	/**
	 * \brief Interval hull of two interval vectors.
	 */
	friend ivec hull(const ivec&, const ivec&);

	/**
	 * \brief Bisect x along a given dimension to lowerhalf and upperhalf.
	 */
	friend ivec lowerhalf(const ivec&, const int);
	friend ivec upperhalf(const ivec&, const int);

	/**
	 * \brief Operation overloads between interval vectors.
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
    inline void ivec::getinf(std::vector<double> &inf) const {
	for (int i = 0; i < _dim; ++i)
	    inf[i] = _itvls[i].getinf();
    }
    inline std::vector<double> ivec::getsup() const {
	std::vector<double> sup(_dim);
	for (int i = 0; i < _dim; ++i) 
	    sup[i] = _itvls[i].getsup();

	return sup;
    }
    inline void ivec::getsup(std::vector<double> &sup) const {
	for (int i = 0; i < _dim; ++i) 
	    sup[i] = _itvls[i].getsup();
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

    inline double ivec::volume() const {
	double v = 1;
	for (int i = 0; i < _dim; ++i)
	    v *= _itvls[i].width();

	return v;
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

    inline ivec& ivec::operator += (const ivec &x) {
	if (!this->isempty()) {
	    assert(this->_dim == x.getdim());
	    for (int i = 0; i != this->_dim; ++i)
		_itvls[i] += x[i];
	}
	return *this;
    }

    inline ivec& ivec::operator += (const std::vector<double> &a) {
	if (!this->isempty()) {
	    assert(this->_dim == int(a.size()));
	    for (int i = 0; i != this->_dim; ++i)
		_itvls[i] += a[i];
	}
	return *this;
    }

    inline ivec& ivec::operator += (const double a) {
	if (!this->isempty()) {
	    for (int i = 0; i != this->_dim; ++i)
		_itvls[i] += a;
	}
	return *this;
    }

    inline ivec& ivec::operator -= (const ivec &x) {
	if (!this->isempty()) {
	    assert(this->_dim == x.getdim());
	    for (int i = 0; i != this->_dim; ++i)
		_itvls[i] -= x[i];
	}
	return *this;
    }

    inline ivec& ivec::operator -= (const std::vector<double> &a) {
	if (!this->isempty()) {
	    assert(this->_dim == int(a.size()));
	    for (int i = 0; i != this->_dim; ++i)
		_itvls[i] -= a[i];
	}
	return *this;
    }

    inline ivec& ivec::operator -= (const double a) {
	if (!this->isempty()) {
	    for (int i = 0; i != this->_dim; ++i)
		_itvls[i] -= a;
	}
	return *this;
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

} // namespace rocs

#endif
