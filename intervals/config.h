#ifndef __config_h_
#define __config_h_

#include <cfloat>
#include <cfenv>
#include <cassert>



//#pragma STDC FENV_ACCESS ON

namespace rocs {

const double EPSIVAL = 1e-6; //DBL_EPSILON; // <= 1e-9

const double PINF = 1.0 / 0.0;
const double NINF = -1.0 / 0.0;
const double NANIVAL = 0.0 / 0.0;

const double PIIVAL = 3.14159265358979323846;
const double PI2IVAL = 2.0 * PIIVAL;
const double PIHALIVAL = PIIVAL / 2.0;


/* set roundup style */
inline void roundup(){ fesetround(FE_UPWARD); }
inline void rounddown(){ fesetround(FE_DOWNWARD); }
inline void roundnear(){ fesetround(FE_TONEAREST); }
inline void roundzero(){ fesetround(FE_TOWARDZERO); }

} // namespace rocs

#endif
