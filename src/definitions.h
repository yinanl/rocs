/**
 *  definitions.h (included in system.hpp)
 *
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
  typedef uint16_t UintSmall;
  typedef std::vector<double> Rn;
  typedef std::vector< std::vector<double> > vecRn;
  typedef std::vector<ivec> vecIv;

  typedef ivec (*fcst)(const ivec&);
  
}  // namespace rocs


#endif
