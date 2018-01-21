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


namespace rocs {

/***** member functions *****/

bool interval::isempty() const{
    return (m_inf > m_sup) | std::isnan(m_inf) | std::isnan(m_sup);
}


interval& interval::operator=(const interval& a){
    
    if (this == &a) {
        return *this;
    }
    
    m_inf = a.m_inf;
    m_sup = a.m_sup;
    
    return *this;
}


/*---------------------
   interval overlap checks: [x],[y] intervals; a (real value)
   [x].isout([y]): [x],[y] are disjoint; 
   [x].isout(a): a is outside of [x];
   [x].isin([y]): [y] is inside [x];
   [x].isin(a): a is inside [x].
---------------------*/
bool interval::isout(const interval &a) const {
    return isempty() || (m_sup < a.getinf()) || (m_inf > a.getsup());
}

bool interval::isout(const double val) const {
    return isempty() || (m_sup < val) || (m_inf > val);
}

bool interval::isin(const interval &a) const {
    return ~isempty() && m_inf <= a.getinf() && m_sup >= a.getsup();
}
 
bool interval::isin(const double val) const {
    
    return ~isempty() && m_inf <= val && m_sup >= val;
}


/*---------------------
   unary operation: -
---------------------*/
interval interval::operator-() const {

    return interval(- m_sup, - m_inf);
}




/***** friend functions *****/

/*---------------------
   rational operations: +, -, *, /:
---------------------*/

/* z = x + y; z = a + x; z = x + a; */
interval operator+(const interval &x, const interval &y){

    if (!(x.isempty()) && !(y.isempty())) {
	
	return interval(x.m_inf + y.m_inf, x.m_sup + y.m_sup);
    }
    else {

	return interval(NAN, NAN);
    }

}

interval operator+(const interval &x, const double a){

    if (!(x.isempty())) {
	
        return interval(x.m_inf + a, x.m_sup + a);
	
    }
    else {
	
        return interval(NAN, NAN);
	
    }
    
}

interval operator+(const double a, const interval &x){
    
    if (!(x.isempty())) {
	
        return interval(x.m_inf + a, x.m_sup + a);
	
    }
    else {
	
        return interval(NAN, NAN);
	
    }
    
}

/* z = x - y; z = a - x; z = x - a; */
interval operator-(const interval &x, const interval &y){

    if (!(x.isempty()) && !(y.isempty())) {
	
        return interval(x.m_inf - y.m_sup, x.m_sup - y.m_inf);
	
    }
    else {
	
        return interval(NAN, NAN);
	
    }
    
}

interval operator-(const interval &x, const double a){
    
    if ( !(x.isempty()) ) {
	
        return interval(x.m_inf - a, x.m_sup - a);
	
    }
    else {
	
        return interval(NAN, NAN);
	
    }
}

interval operator-(const double a, const interval &x){

    if ( !(x.isempty()) ) {
	
        return interval(a - x.m_sup, a - x.m_inf);
	
    }
    else {
	
        return interval(NAN, NAN);
	
    }
    
}

/* mul:
   z = x * y; z = a * x; z = x * a;

   have to discuss when 0 or inf is the bounds, e.g.:

   for real number: 0 * inf = nan, -inf * inf = -inf...
   for interval: [0 inf] * [1, inf] = [0, inf], different from real number
*/
interval operator*(const interval &x, const interval &y){

    if ( x.isempty() || y.isempty() ){

	return interval(NAN, NAN);
	    
    }
    else {

	if ( x.m_inf > 0 && y.m_inf > 0){

	    return interval(x.m_inf * y.m_inf, x.m_sup * y.m_sup);

	}
	else if ( x.m_inf > 0 && y.isin(0) ){

	    if (x.m_sup == PINF) {

		if (std::fabs(y.m_inf) < EPSIVAL) {// [y]=[0,*]

		    return interval(0, PINF);
		}
		else if (std::fabs(y.m_sup) < EPSIVAL) {// [y]=[*,0]

		    return interval(NINF, 0);
		}
		else {

		    return interval(NINF, PINF);
		}
	    }
	    else {
		
		return interval(x.m_sup * y.m_inf, x.m_sup * y.m_sup);
	    }

	}
	else if (x.m_inf > 0 && y.m_sup <0){

	    return interval(x.m_sup * y.m_inf, x.m_inf * y.m_sup);

	}
	else if (x.isin(0) && y.m_inf > 0){

	    if (y.m_sup == PINF) {

		if (std::fabs(x.m_inf) < EPSIVAL) {// [x]=[0,*]

		    return interval(0, PINF);
		}
		else if (std::fabs(x.m_sup) < EPSIVAL) {// [x]=[*,0]

		    return interval(NINF, 0);
		}
		else {

		    return interval(NINF, PINF);
		}
	    }
	    else {
		
		return interval(x.m_inf * y.m_sup, x.m_sup * y.m_sup);
	    }
	}
	else if (x.isin(0) && y.isin(0)){

	    if (x.m_inf == NINF || x.m_sup == PINF ||
		y.m_inf == NINF || y.m_sup == PINF) {

		if ( (x.m_inf >= 0 && y.m_inf >= 0) ||
		     (x.m_sup <=0 && y.m_sup <= 0)) {//[0,inf]*[0,inf] or [-inf,0]*[-inf,0]

		    return interval(0, PINF);
		}
		else if ((x.m_inf >= 0 && y.m_sup <= 0) ||
			 (x.m_sup <=0 && y.m_inf >= 0)) {//[0,inf]*[-inf,0] or switch

		    return interval(NINF, 0);
		}
		else {

		    return interval(NINF, PINF);
		}
	    }
	    else {
		
		return interval( x.m_inf*y.m_sup < x.m_sup*y.m_inf ? x.m_inf*y.m_sup : x.m_sup*y.m_inf,
				 x.m_inf*y.m_inf > x.m_sup*y.m_sup ? x.m_inf*y.m_inf : x.m_sup*y.m_sup );
	    }
	}
	else if (x.isin(0) && y.m_sup < 0){

	    if (y.m_inf == NINF) {

		if (std::fabs(x.m_inf) < EPSIVAL) {// [x]=[0,*]

		    return interval(NINF, 0);
		}
		else if (std::fabs(x.m_sup) < EPSIVAL) {// [x]=[*,0]

		    return interval(0, PINF);
		}
		else {

		    return interval(NINF, PINF);
		}
	    }
	    else {
		
		return interval(x.m_sup * y.m_inf, x.m_inf * y.m_inf);
	    }
	}
	else if (x.m_sup < 0 && y.m_inf > 0){

	    return interval(x.m_inf * y.m_sup, x.m_sup * y.m_inf);
	}
	else if (x.m_sup < 0 && y.isin(0)){

	    if (x.m_inf == NINF) {

		if (std::fabs(y.m_inf) < EPSIVAL) {// [y]=[0,*]

		    return interval(NINF, 0);
		}
		else if (std::fabs(y.m_sup) < EPSIVAL) {// [y]=[*,0]

		    return interval(0, PINF);
		}
		else {

		    return interval(NINF, PINF);
		}
	    }
	    else {
		
		return interval(x.m_inf * y.m_sup, x.m_inf * y.m_inf);
	    }
	}
	else {

	    return interval(x.m_sup * y.m_sup, x.m_inf * y.m_inf);

	}
	
    }
    
}

interval operator*(const double a, const interval &x){

    if ( x.isempty() ) {

	return interval(NAN, NAN);

    }

    else {

	if (a < 0) {

	    return interval(a * x.m_sup, a * x.m_inf);

	}
	else if (a > 0) {

	    return interval(a * x.m_inf, a * x.m_sup);

	}
	else {

	    return interval(0.0, 0.0);

	}
    }

}

interval operator*(const interval &x, const double a){

    return a * x;
}


/* 
   z = x / y; z = a / x; z = x / a;
   [x]/[y], a/[y], [x]/a:
   \empty, if [y]=[0,0].
*/
interval operator/(const interval &x, const interval &y){

    if ( x.isempty() || y.isempty() ||
	 (std::fabs(y.m_inf) < EPSIVAL && std::fabs(y.m_sup) < EPSIVAL)) {

	return interval(NAN, NAN);
    }
    else {

	if (std::fabs(y.m_inf) < EPSIVAL) {// y.inf = 0

	    if (x.m_sup < 0) {

		return interval(NINF, x.m_sup / y.m_sup);
	    }
	    else if (x.m_inf > 0) {
		
		return interval(x.m_inf / y.m_sup, PINF);
	    }
	    else {
		
		return interval(NINF, PINF);
	    }
	}
	else if (std::fabs(y.m_sup) < EPSIVAL) {// y.sup = 0
	    
	    if (x.m_sup < 0) {
		
		return interval(x.m_sup / y.m_inf, PINF);
	    }
	    else if (x.m_inf > 0) {
		
		return interval(NINF, x.m_inf / y.m_inf);
	    }
	    else {
		
		return interval(NINF, PINF);
	    }    
	}
	else if (y.m_inf < 0 && y.m_sup > 0) {// y.inf < 0 < y.sup
	    
	    return interval(NINF, PINF);
	}
	else {
	
	    if ( x.m_inf > 0 && y.m_inf > 0){
		
		return interval(x.m_inf / y.m_sup, x.m_sup / y.m_inf);
	    }
	    else if (x.m_inf > 0 && y.m_sup <0){
		
		return interval(x.m_sup / y.m_sup, x.m_inf / y.m_inf);
	    }
	    else if (x.isin(0) && y.m_inf > 0){
		
		return interval(x.m_inf / y.m_inf, x.m_sup / y.m_inf);
	    }
	    else if (x.isin(0) && y.m_sup < 0){
		
		return interval(x.m_sup / y.m_sup, x.m_inf / y.m_sup);
	    }
	    else if (x.m_sup < 0 && y.m_inf > 0){
		
		return interval(x.m_inf / y.m_inf, x.m_sup / y.m_sup);
	    }
	    else {
		
		return interval(x.m_sup / y.m_inf, x.m_inf / y.m_sup);
	    }
	 
	}
	
    }
    
}

interval operator/(const double a, const interval &y){

    if (y.isempty() ||
	(std::fabs(y.m_inf) < EPSIVAL && std::fabs(y.m_sup) < EPSIVAL)) {

	return interval(NAN, NAN);

    }
    else {
	
	if (std::fabs(y.m_inf) < EPSIVAL) {// y.inf = 0

	    if (a > 0) {

		return interval(a / y.m_sup, PINF);
	    }
	    else if (a < 0) {

		return interval(NINF, a / y.m_sup);
	    }
	    else {

		return interval(NAN, NAN);
	    }
	}
	else if (std::fabs(y.m_sup) < EPSIVAL) {// y.sup = 0

	    if (a > 0) {

		return interval(NINF, a / y.m_inf);
	    }
	    else if (a < 0) {

		return interval(a / y.m_inf, PINF);
	    }
	    else {

		return interval(NAN, NAN);
	    }

	}
	else if (y.m_inf < 0 && y.m_sup > 0) {// y.inf < 0 < y.sup

	    if (std::fabs(a) < EPSIVAL) {

		return interval(NAN, NAN);
	    }
	    else {

		return interval(NINF, PINF);
	    }
	}
	else {// [y] has no 0

	    if (a > 0) {

		return interval(a / y.m_sup, a / y.m_inf);
	    }
	    else if (a < 0) {

		return interval(a / y.m_inf, a / y.m_sup);
	    }
	    else {

		return interval(0, 0);
	    }
	}
    }
    
}

interval operator/(const interval &x, const double a){

    if (x.isempty() || std::fabs(a) < EPSIVAL) {

	return interval(NAN, NAN);
    }
    else if (a < 0) {

	return interval(x.m_sup / a, x.m_inf / a);
    }
    else {

	return interval(x.m_inf / a, x.m_sup / a);
    }
}




/*---------------------
   equivalence test:
   x == y; x != y
   support [-oo,+oo], [-oo,a], [a,+oo], \empty
---------------------*/
bool operator==(const interval &x, const interval &y){

    if (x.isempty() && y.isempty()) {// both are empty

	return true;
    }
    else if (x.isempty() || y.isempty()) {// only one is empty

	return false;
    }
    else {// none is empty
	
	if(!std::isinf(x.m_inf) && !std::isinf(x.m_sup)) {// [a,b]

	    return std::fabs(x.m_inf - y.m_inf) < EPSIVAL && std::fabs(x.m_sup - y.m_sup) < EPSIVAL;
	}
	else if (x.m_inf == NINF && !std::isinf(x.m_sup)) {// [-oo,b]

	    return x.m_inf == y.m_inf && std::fabs(x.m_sup - y.m_sup) < EPSIVAL;
	}
	else if (!std::isinf(x.m_inf) && x.m_sup == PINF) {// [a,+oo]

	    return std::fabs(x.m_inf - y.m_inf) < EPSIVAL && x.m_inf == y.m_inf;
	}
	else {// [-oo,+oo]

	    return x.m_inf == y.m_inf && x.m_sup == y.m_sup;
	}
	
    }

}

bool operator!=(const interval &x, const interval &y)
{
	return !(x == y);
}




/*---------------------
   elementary functions
---------------------*/
/* y = sin(x) */
interval sin(const interval &x){

    interval r;

    if (x.width() > PI2IVAL) {

	r.m_inf = -1;
	r.m_sup = 1;
    }
    else {// shift x.inf to [0,2pi]

	int quot;
	double remu, reml = remquo(x.m_inf, PI2IVAL, &quot);
	
	if (reml < 0) {

	    reml = reml + PI2IVAL;
	    remu = x.m_sup - (quot - 1) * PI2IVAL;
	}
	else {

	    remu = x.m_sup - quot * PI2IVAL;
	}

	if ((reml <= PIHALIVAL && remu >= PIHALIVAL) ||
	    (reml <= 5*PIHALIVAL && remu >= 5*PIHALIVAL)) {

	    r.m_sup = 1;
	}
	else {

	    r.m_sup = std::sin(reml) > std::sin(remu) ? std::sin(reml) : std::sin(remu);
	}

	if ((reml <= 3*PIHALIVAL && remu >= 3*PIHALIVAL) ||
	    (reml <= 7*PIHALIVAL && remu >= 7*PIHALIVAL)) {

	    r.m_inf = -1;
	}
	else {

	    r.m_inf = std::sin(reml) < std::sin(remu) ? std::sin(reml) : std::sin(remu);
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

    if (x.width() >= PIIVAL) {

	return interval(NINF, PINF);
    }
    else {

	/* move x_inf into [-pi/2, pi/2),
	 x_sup into [0, pi) accordingly */
	int k = floor(x.m_inf/PIIVAL);
	if ((x.m_inf-k*PIIVAL) >= PIHALIVAL)
	    k += 1;

	if (std::fabs(x.m_inf-k*PIIVAL + PIHALIVAL) <= EPSIVAL
	    || (x.m_sup-k*PIIVAL) >= PIHALIVAL)
	    return interval(NINF, PINF);
	else
	    return interval(std::tan(x.m_inf-k*PIIVAL), std::tan(x.m_sup-k*PIIVAL));
    }
}


/* y = atan(x), increasing in (-oo, +oo) */
interval atan(const interval &x) {

    return interval(std::atan(x.m_inf), std::atan(x.m_sup));
}


/* y = exp(x) */
interval exp(const interval &x) {

    if (x.isempty()) {

	return interval(NAN, NAN);
    }
    else {

	return interval(std::exp(x.m_inf), std::exp(x.m_sup));
    }
}

/* y = log_e(x) */
interval log(const interval &x) {

    if (x.isempty() || x.m_inf <= 0) {

	return interval(NAN, NAN);
    }
    else {

	return interval(std::log(x.m_inf), std::log(x.m_sup));
    }
}
/* y = log_2(x) */
interval log2(const interval &x) {

    if (x.isempty() || x.m_inf <= 0) {

	return interval(NAN, NAN);
    }
    else {

	return interval(std::log2(x.m_inf), std::log2(x.m_sup));
    }
}


/* y = sqr(x) */
interval sqr(const interval &x) {

    if (x.isempty()) {

	return interval(NAN, NAN);
    }
    else {

	if (x.m_inf > 0) {

	    return interval(x.m_inf * x.m_inf, x.m_sup * x.m_sup);
	}
	else if (x.m_sup < 0) {
	    
	    return interval(x.m_sup * x.m_sup, x.m_inf * x.m_inf);
	}
	else {

	    return interval(0, x.m_sup * x.m_sup > x.m_inf * x.m_inf
			    ? x.m_sup * x.m_sup : x.m_inf * x.m_inf);
	}
    }
    
}

/* y = power(x, n) */
interval power(const interval &x, int n) {

    if (x.isempty()) {

	return interval(NAN, NAN);
    }
    else {

	if (n % 2 == 0) { // n is even

	    if (x.m_inf > 0) {

		return interval(pow(x.m_inf, n), std::pow(x.m_sup, n));
	    }
	    else if (x.m_sup < 0) {
	    
		return interval(pow(x.m_sup, n), std::pow(x.m_inf, n));
	    }
	    else {
		double powl = std::pow(x.m_inf, n);
		double powu = std::pow(x.m_sup, n);
		return interval(0, powu > powl ? powu : powl);
	    }
	}
	else { // n is odd

	    return interval(std::pow(x.m_inf, n), std::pow(x.m_sup, n));
	}
	
    }

}

/* y = sqrt(x) */
interval sqrt(const interval &x) {

    if (x.isempty() || x.m_sup < 0) {

	return interval(NAN, NAN);
    }
    else if (x.m_inf >= 0) {

	return interval(std::sqrt(x.m_inf), std::sqrt(x.m_sup));
    }
    else {

	return interval(0, std::sqrt(x.m_sup));
    }
}




/*--------------------
  intersection, union hull 
--------------------*/
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



/*---------------------
   bisection: split self into right & left half
   self = left;
   return right
---------------------*/
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



/*--------------------
  I/O
--------------------*/
std::ostream& operator<< (std::ostream &out, const interval &a){
    
    out << "[" << a.m_inf << ", " << a.m_sup << "]";

    return out;
}


} // namespace rocs
