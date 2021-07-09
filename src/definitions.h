/**
 *  Typedefs.
 *  
 *  Created by Yinan Li on April. 29, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _definitions_h
#define _definitions_h

#include <vector>
#include "interval_vector.h"


namespace rocs {
    
    typedef uint16_t UintSmall; /**< A small unsigned integer type */
    typedef std::vector<double> Rn; /**< An n-d real vector */
    typedef std::vector< std::vector<double> > vecRn; /**< An array of n-d vectors */
    typedef std::vector<ivec> vecIv; /**< An array of interval vectors */

    typedef ivec (*fcst)(const ivec&); /**< A function pointer */
  
}  // namespace rocs


#endif
