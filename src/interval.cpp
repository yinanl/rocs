/**
 *  interval.cpp
 *
 *  An interval class.
 *
 *  Created by Yinan Li on May 24, 2016.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#include "interval.h"
#include <cassert>
#include <cmath>
#include <iomanip>


namespace rocs {
    
    bool interval::isempty() const {
	return (m_inf > m_sup) | std::isnan(m_inf) | std::isnan(m_sup);
    }


    interval& interval::operator=(const interval& a) {
    
	if (this == &a) {
	    return *this;
	}
    
	m_inf = a.m_inf;
	m_sup = a.m_sup;
    
	return *this;
    }

    
    bool interval::isout(const interval &a) const {
	return isempty() || (m_sup < a.m_inf) || (m_inf > a.m_sup);
    }

    bool interval::isout(const double val) const {
	return isempty() || (m_sup < val) || (m_inf > val);
    }

    bool interval::isin(const interval &a) const {
	return ~isempty() && m_inf <= a.m_inf && m_sup >= a.m_sup;
    }
 
    bool interval::isin(const double val) const {
    
	return ~isempty() && m_inf <= val && m_sup >= val;
    }

    bool operator<(const interval &x, const interval &y) { // 1: x is strictly inside y, 0 otherwise
	return (x.m_inf > y.m_inf) && (x.m_sup < y.m_sup);
    }

    bool operator<=(const interval &x, const interval &y) { // 1: x is inside y, 0 otherwise
	return (x.m_inf >= y.m_inf) && (x.m_sup <= y.m_sup);
    }

    bool operator>(const interval &x, const interval &y) { // 1: x strictly includes y, 0 otherwise
	return (x.m_inf < y.m_inf) && (x.m_sup > y.m_sup);
    }

    bool operator>=(const interval &x, const interval &y) { // 1: x includes y, 0 otherwise
	return (x.m_inf <= y.m_inf) && (x.m_sup >= y.m_sup);
    }
    
    interval interval::operator-() const {

	return interval(- m_sup, - m_inf);
    }
    
    interval interval::operator+() const {

	return *this;
    }


    /**
     * z=x+y, z+=x
     */
    interval& interval::operator+=(const double x) {
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

    interval& interval::operator+=(const interval &x) {
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
    
    interval operator+(const interval &x, const interval &y){
	interval z = x;
	z += y;	
	return z;
    }

    interval operator+(const interval &x, const double a){
	interval z = x;
	z += a;
	return z;    
    }

    interval operator+(const double a, const interval &x){
	interval z = x;
	z += a;
	return z;
    }

    
    /**
     * z = x - y, z = a - x, z = x - a
     */
    interval& interval::operator-=(const double x) {
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

    interval& interval::operator-=(const interval &x) {
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
    
    interval operator-(const interval &x, const interval &y){
	interval z = x;
	z -= y;
	return z;
    }

    interval operator-(const interval &x, const double a){
	interval z = x;
	z -= a;
	return z;
    }

    interval operator-(const double a, const interval &x){

	if ( x.isempty() ) {
	    return interval(NAN, NAN);
	    
	} else {
	    interval z;
	    z.m_inf = sub_RNDD(a, x.m_sup);
	    z.m_sup = sub_RNDU(a, x.m_inf);
	    return z;
	}
    
    }

    /**
     * z = x * y, z *= x
     * have to discuss when 0 or inf is the bounds, e.g.:
     * for real number: 0 * inf = nan, -inf * inf = -inf...
     * for interval: [0 inf] * [1, inf] = [0, inf], different from real number
     */
    interval& interval::operator*=(const double a) {
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

    interval& interval::operator*=(const interval &x) {
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
    
    interval operator*(const interval &x, const interval &y){
	interval z = x;
	z *= y;
	return z;
    }

    interval operator*(const double a, const interval &x) {
	interval z = x;
	z *= a;
	return z;
    }

    interval operator*(const interval &x, const double a) {
	return a * x;
    }


    /**
     * z = x / y; z = a / x; z = x / a;
     * [x]/[y], a/[y], [x]/a:
     * \empty, if [y]=[0,0].
     */
    interval& interval::operator/=(const interval &x) {
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

    interval& interval::operator/=(const double a) {
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
    
    interval operator/(const interval &x, const interval &y){
	interval z = x;
	z /= y;
	return z;
	
    }

    interval operator/(const double a, const interval &y){

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

    interval operator/(const interval &x, const double a){
	interval z = x;
	return z /= a;
    }




    /**
     * equivalence test:
     * x == y; x != y
     * support [-oo,+oo], [-oo,a], [a,+oo], \empty
     */
    bool operator==(const interval &x, const interval &y){

	if (x.isempty() && y.isempty()) {// both are empty
	    return true;
	    
	} else if (x.isempty() || y.isempty()) {// only one is empty
	    return false;
	    
	} else {// none is empty
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

    bool operator!=(const interval &x, const interval &y)
    {
	return !(x == y);
    }




    /**
     * elementary functions
     */
    /* y = sin(x) */
    interval sin(const interval &x){

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


    /* y = cos(x) = sin(x + pi/2) */
    interval cos(const interval &x) {
	return sin(x + PIHALIVAL);
    }


    /* y = tan(x) */
    interval tan(const interval &x) {
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

    interval asin(const interval &x) { // std::asin returns [-pi/2, pi/2], increasing.
	if (x.m_inf < -1.0 || x.m_sup > 1.0)
	    assert("asin() argument out of range.");
	
	interval z;
	z.m_inf = std::asin(x.m_inf);
	z.m_sup = std::asin(x.m_sup);
	z.m_inf = rounddown(z.m_inf);
	z.m_sup = roundup(z.m_sup);
	return z;
    }

    interval acos(const interval &x) { // std::asin returns [0, pi], decreasing.
	if (x.m_inf < -1.0 || x.m_sup > 1.0)
	    assert("acos() argument out of range.");
	
	interval z;
	z.m_inf = std::acos(x.m_sup);
	z.m_sup = std::acos(x.m_inf);
	z.m_inf = rounddown(z.m_inf);
	z.m_sup = roundup(z.m_sup);
	return z;
    }

    /* y = atan(x), increasing in (-oo, +oo) */
    interval atan(const interval &x) {
	interval z;
	z.m_inf = std::atan(x.m_inf);
	z.m_sup = std::atan(x.m_sup);
	z.m_inf = rounddown(z.m_inf);
	z.m_sup = roundup(z.m_sup);
	return z;
    }


    /* y = exp(x) */
    interval exp(const interval &x) {
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

    /* y = log_e(x) */
    interval log(const interval &x) {
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
    /* y = log_2(x) */
    interval log2(const interval &x) {
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

    interval abs(const interval &x) {
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

    /* y = sqr(x) */
    interval sqr(const interval &x) {
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

    /* y = pow(x, n), n can be negative */
    interval pow(const interval &x, const int n) {
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

    interval pow(const interval &x, const interval &y) {
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

    interval pow(const interval &x, const double y) {
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

    /* y = sqrt(x) */
    interval sqrt(const interval &x) {
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




    /**
     * intersection, union hull 
     */
    interval intersect(const interval &x, const interval &y) {

	if (x.isout(y)) {

	    return interval(NAN, NAN);
	}
	else if (x.isin(y)) {

	    return y;
	}
	else if (y.isin(x)) {

	    return x;
	}
	else {

	    return interval(x.m_inf > y.m_inf ? x.m_inf : y.m_inf,
			    x.m_sup < y.m_sup ? x.m_sup : y.m_sup);
	}
    }

    interval hull(const interval &x, const interval &y) {

	if (x.isempty()) {

	    return y;
	}
	else if (y.isempty()) {

	    return x;
	}
	else {

	    return interval(x.m_inf < y.m_inf ? x.m_inf : y.m_inf,
			    x.m_sup > y.m_sup ? x.m_sup : y.m_sup);
	}
    }



    /**
     * bisection: split self into right & left half
     * self = left;
     * return right
     */
    interval lowerhalf(const interval &self) {

	if (self.isempty()) {

	    return self;
	}
	else {
	
	    return interval(self.m_inf, self.mid());
	}
    }

    interval upperhalf(const interval &self) {

	if (self.isempty()) {

	    return self;
	}
	else {

	    return interval(self.mid(), self.m_sup);
	}
    }



    /**
     * I/O
     */
    std::ostream& operator<< (std::ostream &out, const interval &a){
    
	out << "[" << a.m_inf << ", " << a.m_sup << "]";

	return out;
    }


} // namespace rocs
