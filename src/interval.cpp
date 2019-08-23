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
    /**
     * I/O
     */
    std::ostream& operator<< (std::ostream &out, const interval &a){
    
	out << "[" << a.m_inf << ", " << a.m_sup << "]";

	return out;
    }


} // namespace rocs
