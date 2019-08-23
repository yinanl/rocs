/**
 *  interval.h
 *
 *  An interval class.
 *
 *  Created by Yinan Li on May 24, 2016.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _interval_h_
#define _interval_h_

#include <iostream>
#include <cmath>
#include "config.h"
#include "validated.h"


namespace rocs {
  
    class interval
    {
    private:
	double m_inf; /**< the lower bound */
	double m_sup; /**< the upper bound */
        
    public:

	/**
	 * \brief Constructors.
	 */
	interval(): m_inf(0), m_sup(0){};
	interval(double l, double u): m_inf(l), m_sup(u){};
	interval(double r): m_inf(r), m_sup(r){};
	/**
	 * \brief A copy constructor.
	 * @param a another interval.
	 */
	interval(const interval &a) : m_inf(a.m_inf), m_sup(a.m_sup) {}

	/**
	 * \brief A copy assignment.
	 */
	interval& operator=(const interval &a) {
	    if (this != &a) {
		m_inf = a.m_inf;
		m_sup = a.m_sup;
	    }
	    return *this;
	}

	/**
	 * \brief Empty set test.
	 */
	bool isempty() const {
	    return (m_inf > m_sup) | std::isnan(m_inf) | std::isnan(m_sup);
	}
	/**
	 * \brief Return lower/upper bound of (*this).
	 */
	double getinf() const { return m_inf; }
	double getsup() const { return m_sup; }
	/**
	 * \brief Set values for lower/upper bound of (*this).
	 */
	void setinf(const double val) { m_inf= val; }
	void setsup(const double val) { m_sup= val; }
	/**
	 * \brief Return the width of (*this).
	 */
	double width() const { return roundup(m_sup - m_inf); }
	/**
	 * \brief Return the center point of (*this).
	 */
	double mid() const { return roundup(m_inf + (m_sup - m_inf)/2.0); }

	/**
	 * \brief [x].isout([a]) is true if [x],[a] are disjoint.
	 */
	bool isout(const interval &a) const {
	    return isempty() || (m_sup < a.m_inf) || (m_inf > a.m_sup);
	}
	/**
	 * \brief [x].isout(val) is true if a real value val is outside of [x].
	 */
	bool isout(const double val) const {
	    return isempty() || (m_sup < val) || (m_inf > val);
	}
	/**
	 * \brief [x].isin([a]) is true if [a] is inside [x].
	 */
	bool isin(const interval &a) const { // if a is contained in it
	    return ~isempty() && m_inf <= a.m_inf && m_sup >= a.m_sup;
	}
	/**
	 * \brief [x].isin([a]) is true if a real value val is inside [x].
	 */
	bool isin(const double val) const {
	    return ~isempty() && m_inf <= val && m_sup >= val;
	}

	/**
	 * \brief Return *this.
	 */
	interval& operator+() {return *this;}
	/**
	 * \brief Negation
	 * Return a new interval with lower and upper bounds negated.
	 */
	interval operator-();
        /**
	 * \brief Add a real value x to (*this) interval.
	 */
	interval& operator += (const double x);
	/**
	 * \brief Add an interval [x] to (*this) interval.
	 */
	interval& operator += (const interval &x);
	/**
	 * \brief Subtract a real value x to (*this) interval.
	 */
	interval& operator -= (const double x);
	/**
	 * \brief Subtract an interval [x] to (*this) interval.
	 */
	interval& operator -= (const interval &x);
	/**
	 * \brief Multiply a real value x to (*this) interval.
	 */
	interval& operator *= (const double x);
	/**
	 * \brief Multiply an interval [x] to (*this) interval.
	 */
	interval& operator *= (const interval &x);
	/**
	 * \brief Divide (*this) interval by a real value x.
	 */
	interval& operator /= (const double x);
	/**
	 * \brief Divide (*this) interval by an interval [x].
	 */
	interval& operator /= (const interval &x);

	/**
	 * \brief x == y iff they have the same lower and upper bounds. Otherwise, x != y.
	 */
	friend bool operator==(const interval&, const interval&);
	friend bool operator!=(const interval&, const interval&);

	/**
	 * \brief [x]<[y] if [x] is strictly included in [y].
	 */
	friend bool operator<(const interval&x, const interval& y);
	/**
	 * \brief [x]<=[y] if [x] is included in [y].
	 */
	friend bool operator<=(const interval&x, const interval& y);
	/**
	 * \brief [x]>[y] if [y] is strictly included in [x].
	 */
	friend bool operator>(const interval&x, const interval& y);
	/**
	 * \brief [x]>=[y] if [y] is included in [x].
	 */
	friend bool operator>=(const interval&x, const interval& y);

	/**
	 * Interval arithmetics
	 */
	friend interval operator+(const interval&, const interval&);
	friend interval operator+(const double, const interval&);
	friend interval operator+(const interval&, const double);
	friend interval operator*(const interval&, const interval&);
	friend interval operator*(const double, const interval&);
	friend interval operator*(const interval&, const double);
	friend interval operator-(const interval&, const interval&);
	friend interval operator-(const double, const interval&);
	friend interval operator-(const interval&, const double);
	friend interval operator/(const interval&, const interval&);
	friend interval operator/(const double, const interval&);
	friend interval operator/(const interval&, const double);
  
	friend bool operator==(const interval&, const interval&);
	friend bool operator!=(const interval&, const interval&);
	friend bool operator<(const interval&x, const interval& y);
	friend bool operator<=(const interval&x, const interval& y);
	friend bool operator>(const interval&x, const interval& y);
	friend bool operator>=(const interval&x, const interval& y);

	friend interval sin(const interval&); // y = sin(x)
	friend interval cos(const interval&); // y = cos(x)
	friend interval tan(const interval&); // y = tan(x)
	friend interval atan(const interval&); // y = atan(x)
	friend interval asin(const interval&); // y = asin(x)
	friend interval acos(const interval&); // y = acos(x)
	friend interval exp(const interval&); // y = e^x
	friend interval log(const interval&); // y = log_e(x)
	friend interval log2(const interval&); // y = log_2(x)
	friend interval abs(const interval&); // y=|x|
	friend interval sqrt(const interval&); // y = sqrt(x)
	friend interval sqr(const interval&); // y = x^2
	friend interval pow(const interval&, const int); // y = x^n
	friend interval pow(const interval &x, const interval &y); // r=x^y
	friend interval pow(const interval &x, const double y);
	friend interval nthroot(const interval&, int); // y = x^(1/n)

	/**
	 * \brief The intersection of two intervals.
	 */
	friend interval intersect(const interval &, const interval &);
	/**
	 * \brief The union hall of two intervals.
	 */
	friend interval hull(const interval&, const interval&);
	/**
	 * \brief Return the lower half of the input interval.
	 */
	friend interval lowerhalf(const interval&);
	/**
	 * \brief Return the upper half of the input interval.
	 */
	friend interval upperhalf(const interval&);
	
	/**
	 * \brief I/O
	 */
	friend std::ostream& operator<<(std::ostream&, const interval&);
    };


    
    /** Inline member functions */
    inline interval interval::operator-() {
	return interval(- m_sup, - m_inf);
    }
    
    inline interval& interval::operator+=(const double x) {
	if( this->isempty() ) {
	    m_inf = NAN;
	    m_sup = NAN;
	} else {
	    m_inf = add_RNDD(m_inf, x);
	    m_sup = add_RNDU(m_sup, x);
	    roundnear();
	}
	
	return *this;
    }

    inline interval& interval::operator+=(const interval &x) {
	if( this->isempty() || x.isempty()) {
	    m_inf = NAN;
	    m_sup = NAN;	    
	} else {
	    m_inf = add_RNDD(m_inf, x.m_inf);
	    m_sup = add_RNDU(m_sup, x.m_sup);
	    roundnear();
	}

	return *this;
    }

    inline interval& interval::operator-=(const double x) {
    	if( this->isempty() ) {
	    m_inf = NAN;
	    m_sup = NAN;
	} else {
	    m_inf = sub_RNDD(m_inf, x);
	    m_sup = sub_RNDU(m_sup, x);
	    roundnear();
	}
	
	return *this;
    }

    inline interval& interval::operator-=(const interval &x) {
	if( this->isempty() || x.isempty()) {
	    m_inf = NAN;
	    m_sup = NAN;
	} else {
	    m_inf = sub_RNDD(m_inf, x.m_sup);
	    m_sup = sub_RNDU(m_sup, x.m_inf);
	    roundnear();
	}

	return *this;
    }

    inline interval& interval::operator*=(const double a) {
	if (this->isempty()) {
	    m_inf = NAN;
	    m_sup = NAN;
	    
	} else {
	    double inf = m_inf;
	    double sup = m_sup;
	    if (a < 0) {
		m_inf = mul_RNDD(sup, a);
		m_sup = mul_RNDU(inf, a);
		
	    } else if (a > 0) {
		m_inf = mul_RNDD(inf, a);
		m_sup = mul_RNDU(sup, a);
		
	    } else {
		m_inf = 0.0;
		m_sup = 0.0;
	    }
	}

	return *this;
    }

    inline interval& interval::operator*=(const interval &x) {
	if ( this->isempty() || x.isempty() ) {
	    m_inf = NAN;
	    m_sup = NAN;
	    
	} else {
	    double inf = m_inf;
	    double sup = m_sup;
	    
	    if ( m_inf > 0 && x.m_inf > 0) {
		m_inf = mul_RNDD(inf, x.m_inf);
		m_sup = mul_RNDU(sup, x.m_sup);
		
	    } else if ( m_inf > 0 && x.isin(0) ) {
		if (m_sup == PINF) {
		    if (std::fabs(x.m_inf) < EPSMACHINE) {// [y]=[0,*]
			m_inf = 0.0;
			m_sup = PINF;
		    } else if (std::fabs(x.m_sup) < EPSMACHINE) {// [y]=[*,0]
			m_inf = NINF;
			m_sup = 0.0;
		    } else {
			m_inf = NINF;
			m_sup = PINF;
		    }
		} else {
		    m_inf = mul_RNDD(sup, x.m_inf);
		    m_sup = mul_RNDU(sup, x.m_sup);
		}

	    } else if (m_inf > 0 && x.m_sup <0) {
		m_inf = mul_RNDD(sup, x.m_inf);
		m_sup = mul_RNDU(inf, x.m_sup);

	    } else if (this->isin(0) && x.m_inf > 0) {
		if (x.m_sup == PINF) {
		    if (std::fabs(m_inf) < EPSMACHINE) {// [x]=[0,*]
			m_inf = 0.0;
			m_sup = PINF;
		    } else if (std::fabs(m_sup) < EPSMACHINE) {// [x]=[*,0]
			m_inf = NINF;
			m_sup = 0.0;
		    } else {
			m_inf = NINF;
			m_sup = PINF;
		    }
		} else {
		    m_inf = mul_RNDD(inf, x.m_sup);
		    m_sup = mul_RNDU(sup, x.m_sup);
		}
		
	    } else if (this->isin(0) && x.isin(0)) {
		if (x.m_inf == NINF || x.m_sup == PINF ||
		    m_inf == NINF || m_sup == PINF) {

		    if ( (x.m_inf >= 0 && m_inf >= 0) ||
			 (x.m_sup <=0 && m_sup <= 0)) {//[0,inf]*[0,inf] or [-inf,0]*[-inf,0]
			m_inf = 0.0;
			m_sup = PINF;
		    } else if ((m_inf >= 0 && x.m_sup <= 0) ||
			     (m_sup <=0 && x.m_inf >= 0)) {//[0,inf]*[-inf,0] or switch
			m_inf = NINF;
			m_sup = 0.0;
		    } else {
			m_inf = NINF;
			m_sup = PINF;
		    }
		} else {
		    m_inf = m_inf*x.m_sup < m_sup*x.m_inf ?
			mul_RNDD(inf, x.m_sup) : mul_RNDD(sup, x.m_inf);
		    m_sup = m_inf*x.m_inf > m_sup*x.m_sup ?
			mul_RNDD(inf, x.m_inf) : mul_RNDU(sup, x.m_sup);
		}
		
	    } else if (this->isin(0) && x.m_sup < 0) {
		if (x.m_inf == NINF) {
		    if (std::fabs(m_inf) < EPSMACHINE) {// [x]=[0,*]
			m_inf = NINF;
			m_sup = 0.0;
		    } else if (std::fabs(m_sup) < EPSMACHINE) {// [x]=[*,0]
			m_inf = 0.0;
			m_sup = PINF;
		    } else {
			m_inf = NINF;
			m_sup = PINF;
		    }
		} else {
		    m_inf = mul_RNDD(sup, x.m_inf);
		    m_sup = mul_RNDU(inf, x.m_inf);
		}
		
	    } else if (m_sup < 0 && x.m_inf > 0) {
		m_inf = mul_RNDD(inf, x.m_sup);
		m_sup = mul_RNDU(sup, x.m_inf);
		
	    } else if (m_sup < 0 && x.isin(0)) {
		if (m_inf == NINF) {
		    if (std::fabs(x.m_inf) < EPSMACHINE) {// [y]=[0,*]
			m_inf = NINF;
			m_sup = 0.0;
		    } else if (std::fabs(x.m_sup) < EPSMACHINE) {// [y]=[*,0]
			m_inf = 0.0;
			m_sup = PINF;
		    } else {
			m_inf = NINF;
			m_sup = PINF;
		    }
		} else {
		    m_inf = mul_RNDD(inf, x.m_sup);
		    m_sup = mul_RNDU(inf, x.m_inf);
		}
		
	    } else {
		m_inf = mul_RNDD(sup, x.m_sup);
		m_sup = mul_RNDU(inf, x.m_inf);
	    }
	}

	return *this;
    }

    inline interval& interval::operator/=(const interval &x) {
	if ( this->isempty() || x.isempty() ||
	     (std::fabs(x.m_inf) < EPSMACHINE && std::fabs(x.m_sup) < EPSMACHINE)) {
	    m_inf = NAN;
	    m_sup = NAN;
	} else {
	    double inf = m_inf;
	    double sup = m_sup;
	    if (std::fabs(x.m_inf) < EPSMACHINE) {// y.inf = 0
		if (m_sup < 0) {
		    m_inf = NINF;
		    m_sup = div_RNDU(sup, x.m_sup);
		} else if (m_inf > 0) {
		    m_inf = div_RNDD(inf, x.m_sup);
		    m_sup = PINF;
		} else {
		    m_inf = NINF;
		    m_sup = PINF;
		}
		
	    } else if (std::fabs(x.m_sup) < EPSMACHINE) {// y.sup = 0
		if (m_sup < 0) {
		    m_inf = div_RNDD(sup, x.m_inf);
		    m_sup = PINF;
		} else if (m_inf > 0) {
		    m_inf = NINF;
		    m_sup = div_RNDU(inf, x.m_inf);
		} else {
		    m_inf = NINF;
		    m_sup = PINF;
		}
		
	    } else if (x.m_inf < 0 && x.m_sup > 0) {// y.inf < 0 < y.sup
		m_inf = NINF;
		m_sup = PINF;
		
	    } else {	
		if ( m_inf > 0 && x.m_inf > 0) {
		    m_inf = div_RNDD(inf, x.m_sup);
		    m_sup = div_RNDU(sup, x.m_inf);
		} else if (m_inf > 0 && x.m_sup <0) {
		    m_inf = div_RNDD(sup, x.m_sup);
		    m_sup = div_RNDU(inf, x.m_inf);
		} else if (this->isin(0) && x.m_inf > 0) {
		    m_inf = div_RNDD(inf, x.m_inf);
		    m_sup = div_RNDU(sup, x.m_inf);
		} else if (this->isin(0) && x.m_sup < 0) {
		    m_inf = div_RNDD(sup, x.m_sup);
		    m_sup = div_RNDU(inf, x.m_sup);
		} else if (m_sup < 0 && x.m_inf > 0) {
		    m_inf = div_RNDD(inf, x.m_inf);
		    m_sup = div_RNDU(sup, x.m_sup);
		} else {
		    m_inf = div_RNDD(sup, x.m_inf);
		    m_sup = div_RNDU(inf, x.m_sup);
		}
	    }
	
	}

	return *this;
    }

    inline interval& interval::operator/=(const double a) {
	double inf = m_inf;
	double sup = m_sup;
	
	if (this->isempty() || std::fabs(a) < EPSMACHINE) {
	    m_inf = NAN;
	    m_sup = NAN;
	    
	} else if (a < 0) {
	    m_inf = div_RNDD(sup, a);
	    m_sup = div_RNDU(inf, a);
	    
	} else {
	    m_inf = div_RNDD(inf, a);
	    m_sup = div_RNDU(sup, a);
	}
	return *this;
    }

    /** Inline non-member functions. */
    inline bool operator<(const interval &x, const interval &y) {
	return !y.isempty() && (x.m_inf > y.m_inf) && (x.m_sup < y.m_sup);
    }

    inline bool operator<=(const interval &x, const interval &y) {
	return !y.isempty() && (x.m_inf >= y.m_inf) && (x.m_sup <= y.m_sup);
    }

    inline bool operator>(const interval &x, const interval &y) {
	return !x.isempty() && (x.m_inf < y.m_inf) && (x.m_sup > y.m_sup);
    }

    inline bool operator>=(const interval &x, const interval &y) {
	return !x.isempty() && (x.m_inf <= y.m_inf) && (x.m_sup >= y.m_sup);
    }

    inline bool operator==(const interval &x, const interval &y) {
	if (x.isempty() && y.isempty()) {
	    return true;	    
	} else if (x.isempty() || y.isempty()) {
	    return false;
	} else {
	    if(!std::isinf(x.m_inf) && !std::isinf(x.m_sup)) {// [a,b]
		return std::fabs(x.m_inf - y.m_inf) < EPSIVAL && std::fabs(x.m_sup - y.m_sup) < EPSIVAL;
		
	    } else if (x.m_inf == NINF && !std::isinf(x.m_sup)) {// [-oo,b]
		return x.m_inf == y.m_inf && std::fabs(x.m_sup - y.m_sup) < EPSIVAL;
		
	    } else if (!std::isinf(x.m_inf) && x.m_sup == PINF) {// [a,+oo]
		return std::fabs(x.m_inf - y.m_inf) < EPSIVAL && x.m_sup == y.m_sup;
		
	    } else {// [-oo,+oo]
		return x.m_inf == y.m_inf && x.m_sup == y.m_sup;
	    }
	}

    }

    inline bool operator!=(const interval &x, const interval &y) {
	return !(x == y);
    }

    inline interval operator+(const interval &x, const interval &y) {
	return interval(x) += y;
    }

    inline interval operator+(const interval &x, const double a) {
	return interval(x) += a;
    }

    inline interval operator+(const double a, const interval &x) {
	return interval(x) += a;
    }

    inline interval operator-(const interval &x, const interval &y) {
	return interval(x) -= y;
    }

    inline interval operator-(const interval &x, const double a) {
	return interval(x) -= a;
    }

    inline interval operator-(const double a, const interval &x) {
	if ( x.isempty() ) {
	    return interval(NAN, NAN);
	    
	} else {
	    interval z;
	    z.m_inf = sub_RNDD(a, x.m_sup);
	    z.m_sup = sub_RNDU(a, x.m_inf);
	    return z;
	}
    
    }
    
    inline interval operator*(const interval &x, const interval &y) {
	return interval(x) *= y;
    }

    inline interval operator*(const double a, const interval &x) {
	return interval(x) *= a;
    }

    inline interval operator*(const interval &x, const double a) {
	return a * x;
    }

    inline interval operator/(const interval &x, const interval &y) {
	return interval(x) /= y;
    }

    inline interval operator/(const double a, const interval &y) {
	if (y.isempty() ||
	    (std::fabs(y.m_inf) < EPSMACHINE && std::fabs(y.m_sup) < EPSMACHINE)) {
	    return interval(NAN, NAN);

	} else {
	    interval z;
	    if (std::fabs(y.m_inf) < EPSMACHINE) {// y.inf = 0
		if (a > 0) {
		    z.m_inf = div_RNDD(a, y.m_sup);
		    z.m_sup = PINF;
		} else if (a < 0) {
		    z.m_inf = NINF;
		    z.m_sup = div_RNDU(a, y.m_sup);
		} else {
		    z.m_inf = NAN;
		    z.m_sup = NAN;
		}
	    } else if (std::fabs(y.m_sup) < EPSMACHINE) {// y.sup = 0
		if (a > 0) {
		    z.m_inf = NINF;
		    z.m_sup = div_RNDU(a, y.m_inf);
		} else if (a < 0) {
		    z.m_inf = div_RNDD(a, y.m_inf);
		    z.m_sup = PINF;
		} else {
		    z.m_inf = NAN;
		    z.m_sup = NAN;
		}
	    } else if (y.m_inf < 0 && y.m_sup > 0) {// y.inf < 0 < y.sup
		if (std::fabs(a) < EPSMACHINE) {
		    z.m_inf = NAN;
		    z.m_sup = NAN;
		} else {
		    z.m_inf = NINF;
		    z.m_sup = PINF;
		}
	    } else {// [y] has no 0
		if (a > 0) {
		    z.m_inf = div_RNDD(a, y.m_sup);
		    z.m_sup = div_RNDU(a, y.m_inf);
		} else if (a < 0) {
		    z.m_inf = div_RNDD(a, y.m_inf);
		    z.m_sup = div_RNDU(a, y.m_sup);
		} else {
		    z.m_inf = 0.0;
		    z.m_sup = 0.0;
		}
	    }
	    return z;
	}
    }

    inline interval operator/(const interval &x, const double a) {
	return interval(x) /= a;
    }

    inline interval sin(const interval &x) {
	interval r;

	if (x.width() > PI2IVAL) {
	    r.m_inf = -1;
	    r.m_sup = 1;
	    
	} else { // shift x.inf to [0,2pi]

	    int q;
	    double remu, reml;
	    reml = std::remquo(x.m_inf, PI2IVAL, &q);
	    reml = rounddown(reml);
	    
	    if (reml < 0) {
		reml = add_RNDD(reml, PI2IVAL);
		remu = sub_RNDU(x.m_sup, mul_RNDD((q - 1), PI2IVAL));
		
	    } else {
		remu = sub_RNDU(x.m_sup, mul_RNDD(q, PI2IVAL));
	    }

	    if ((reml <= PIHALIVAL && remu >= PIHALIVAL) ||
		(reml <= 5*PIHALIVAL && remu >= 5*PIHALIVAL)) {
		r.m_sup = 1;
		
	    } else {
		r.m_sup = std::sin(reml) > std::sin(remu) ? std::sin(reml) : std::sin(remu);
		r.m_sup = roundup(r.m_sup);
	    }

	    if ((reml <= 3*PIHALIVAL && remu >= 3*PIHALIVAL) ||
		(reml <= 7*PIHALIVAL && remu >= 7*PIHALIVAL)) {
		r.m_inf = -1;
		
	    } else {
		r.m_inf = std::sin(reml) < std::sin(remu) ? std::sin(reml) : std::sin(remu);
		r.m_inf = rounddown(r.m_inf);
	    }
	}

	return r;
    }

    inline interval cos(const interval &x) { // y = cos(x) = sin(x + pi/2)
	return sin(x + PIHALIVAL);
    }

    inline interval tan(const interval &x) {
	interval r;
	if (x.width() >= PIIVAL) {
	    r.m_inf = NINF;
	    r.m_sup = PINF;
	    
	} else {
	    /* move x_inf into [-pi/2, pi/2),
	       x_sup into [0, pi) accordingly */
	    int k = floor(x.m_inf/PIIVAL);
	    if ((x.m_inf-k*PIIVAL) >= PIHALIVAL) {
		k += 1;
	    }

	    if (std::fabs(x.m_inf-k*PIIVAL + PIHALIVAL) <= EPSMACHINE
		|| (x.m_sup-k*PIIVAL) >= PIHALIVAL) {
		r.m_inf = NINF;
		r.m_sup = PINF;

	    } else {
		r.m_inf = std::tan(x.m_inf-k*PIIVAL);
		r.m_sup = std::tan(x.m_sup-k*PIIVAL);
		r.m_inf = rounddown(r.m_inf);
		r.m_sup = roundup(r.m_sup);
	    }
	}

	return r;
    }

    inline interval asin(const interval &x) { // std::asin returns [-pi/2, pi/2], increasing.
	assert(x.m_inf < -1.0 || x.m_sup > 1.0);
	interval z;
	z.m_inf = std::asin(x.m_inf);
	z.m_sup = std::asin(x.m_sup);
	z.m_inf = rounddown(z.m_inf);
	z.m_sup = roundup(z.m_sup);
	return z;
    }

    inline interval acos(const interval &x) { // std::asin returns [0, pi], decreasing.
	assert(x.m_inf < -1.0 || x.m_sup > 1.0);
	interval z;
	z.m_inf = std::acos(x.m_sup);
	z.m_sup = std::acos(x.m_inf);
	z.m_inf = rounddown(z.m_inf);
	z.m_sup = roundup(z.m_sup);
	return z;
    }

    inline interval atan(const interval &x) {
	interval z;
	z.m_inf = std::atan(x.m_inf);
	z.m_sup = std::atan(x.m_sup);
	z.m_inf = rounddown(z.m_inf);
	z.m_sup = roundup(z.m_sup);
	return z;
    }

    inline interval exp(const interval &x) {
	interval r;
	if (x.isempty()) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    
	} else {
	    r.m_inf = std::exp(x.m_inf);
	    r.m_sup = std::exp(x.m_sup);
	    r.m_inf = rounddown(r.m_inf);
	    r.m_sup = roundup(r.m_sup);

	    if (r.m_inf < 0.0)
		r.m_inf = 0.0;
	}
	return r;
    }

    inline interval log(const interval &x) {
	interval r;
	if (x.isempty() || x.m_inf <= 0.0) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    
	} else {
	    r.m_inf = std::log(x.m_inf);
	    r.m_sup = std::log(x.m_sup);
	    r.m_inf = rounddown(r.m_inf);
	    r.m_sup = roundup(r.m_sup);
	}
	return r;
    }
    
    inline interval log2(const interval &x) {
	interval r;
	if (x.isempty() || x.m_inf <= 0.0) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    
	} else {
	    r.m_inf = std::log2(x.m_inf);
	    r.m_sup = std::log2(x.m_sup);
	    r.m_inf = rounddown(r.m_inf);
	    r.m_sup = roundup(r.m_sup);
	}
	return r;
    }

    inline interval abs(const interval &x) {
	interval r;
	if (x.isempty()) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    
	} else {
	    if (x.m_inf > 0) {
		r.m_inf = x.m_inf;
		r.m_sup = x.m_sup;
		
	    } else if (x.m_sup < 0) {
		r.m_inf = -x.m_sup;
		r.m_sup = -x.m_inf;
		
	    } else {
		r.m_inf = 0.0;
		r.m_sup = -x.m_inf > x.m_sup ? -x.m_inf : x.m_sup;
	    }
	}
	return r;
    }

    inline interval sqr(const interval &x) {
	interval r;
	// r = abs(x);
	// r *= r;
	if (x.isempty()) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    
	} else {
	    if (x.m_inf > 0) {
		r.m_inf = mul_RNDD(x.m_inf, x.m_inf);
		r.m_sup = mul_RNDU(x.m_sup, x.m_sup);
		
	    } else if (x.m_sup < 0) {
		r.m_inf = mul_RNDD(x.m_sup, x.m_sup);
		r.m_sup = mul_RNDU(x.m_inf, x.m_inf);
		
	    } else {
		r.m_inf = 0.0;
		r.m_sup = -x.m_inf > x.m_sup ?
		    mul_RNDU(x.m_inf, x.m_inf) : mul_RNDU(x.m_sup, x.m_sup);
	    }
		
	}
	return r;
    }

    inline interval pow(const interval &x, const int n) {
	interval r;
	int p = (n < 0) ? (-n) : n;

	if (x.isempty()) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    
	} else {
	    if (p % 2 == 0) { // n is even
		if (x.m_inf > 0) {
		    rounddown();
		    r.m_inf = std::pow(x.m_inf, n);
		    roundup();
		    r.m_sup = std::pow(x.m_sup, n);
		    roundnear();
		} else if (x.m_sup < 0) {
		    rounddown();
		    r.m_inf = std::pow(x.m_sup, n);
		    roundup();
		    r.m_sup = std::pow(x.m_inf, n);
		    roundnear();
		} else {
		    r.m_inf = 0.0;
		    roundup();
		    r.m_sup = -x.m_inf > x.m_sup ? std::pow(x.m_inf, n) : std::pow(x.m_sup, n);
		    roundnear();
		}
	    
	    } else { // n is odd
		rounddown();
		r.m_inf = std::pow(x.m_inf, n);
		roundup();
		r.m_sup = std::pow(x.m_sup, n);
		roundnear();
	    }

	    if (n < 0)
		r = 1 / r;
	}
	
	return r;
    }

    inline interval pow(const interval &x, const interval &y) {
	interval r;
	if (x.isempty() || x.m_inf < 0.0) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    return r;
	    
	} else {
	    if (std::fabs(x.m_inf) <= EPSMACHINE) {
		if ( y.m_inf <= 0.0 ) {
		    r.m_inf = NAN;
		    r.m_sup = NAN;
		    return r;
		} else if (std::fabs(x.m_sup) <= EPSMACHINE) {
		    r.m_inf = 0.0;
		    r.m_sup = 0.0;
		    return r;
		} else {
		    r.m_inf = 0.0;
		    r.m_sup = std::log(x.m_sup);
		    r.m_sup = roundup(r.m_sup);
		}
		
	    } else {
		r = log(x);
	    }

	    r *= y;
	    if (std::fabs(x.m_inf) <= EPSMACHINE) {
		r.m_inf = 0.0;
		r.m_sup = std::exp(r.m_sup);
		r.m_sup = roundup(r.m_sup);
		return r;
	    } else {
		return exp(r);
	    }
	}
	
    }

    inline interval pow(const interval &x, const double y) {
	interval r;
	if (x.isempty() || x.m_inf < 0.0) { // if x contains negative value
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    return r;
	    
	} else {
	    if (std::fabs(x.m_inf) <= EPSMACHINE) { // if x_inf=0
		if ( y <= 0.0 ) {
		    r.m_inf = NAN;
		    r.m_sup = NAN;
		    return r;
		} else if (std::fabs(x.m_sup) <= EPSMACHINE) {// if y>0 & x=[0,0]
		    r.m_inf = 0.0;
		    r.m_sup = 0.0;
		    return r;
		} else {
		    r.m_inf = 0.0;
		    r.m_sup = std::log(x.m_sup);
		    r.m_sup = roundup(r.m_sup);
		}
		
	    } else {
		r = log(x);
	    }

	    r *= y;
	    if (std::fabs(x.m_inf) <= EPSMACHINE) {
		r.m_inf = 0.0;
		r.m_sup = std::exp(r.m_sup);
		r.m_sup = roundup(r.m_sup);
		return r;
	    } else {
		return exp(r);
	    }
	}
	
    }

    inline interval sqrt(const interval &x) {
	interval r;
	if (x.isempty() || x.m_sup < 0) {
	    r.m_inf = NAN;
	    r.m_sup = NAN;
	    
	} else if (x.m_inf >= 0) {
	    r.m_inf = std::sqrt(x.m_inf);
	    r.m_sup = std::sqrt(x.m_sup);
	    r.m_inf = rounddown(r.m_inf);
	    r.m_sup = roundup(r.m_sup);
	}
	else {
	    r.m_inf = 0.0;
	    r.m_sup = std::sqrt(x.m_sup);
	    r.m_sup = roundup(r.m_sup);
	}
	return r;
    }

    inline interval intersect(const interval &x, const interval &y) {
	if (x.isout(y)) 
	    return interval(NAN, NAN);
	else if (x.isin(y))
	    return y;
	else if (y.isin(x)) 
	    return x;
	else {
	    return interval(x.m_inf > y.m_inf ? x.m_inf : y.m_inf,
			    x.m_sup < y.m_sup ? x.m_sup : y.m_sup);
	}
    }

    inline interval hull(const interval &x, const interval &y) {
	if (x.isempty())
	    return y;
	else if (y.isempty())
	    return x;
	else {
	    return interval(x.m_inf < y.m_inf ? x.m_inf : y.m_inf,
			    x.m_sup > y.m_sup ? x.m_sup : y.m_sup);
	}
    }
    
    inline interval lowerhalf(const interval &self) {
	if (self.isempty())
	    return self;
	else 
	    return interval(self.m_inf, self.mid());
    }

    inline interval upperhalf(const interval &self) {
	if (self.isempty())
	    return self;
	else 
	    return interval(self.mid(), self.m_sup);
    }
    
    inline interval invhull(double x, double y) {
	double inf, sup;
	if (x<=y) {
	    inf = x;     sup = y;
	} else {
	    inf = y;     sup = x;
	}
	return interval(inf, sup);
    }
 

} // namespace rocs

#endif
