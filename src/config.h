#ifndef __config_h_
#define __config_h_

#include <cfloat>
#include <cfenv>
#include <cassert>
#include <limits>


//#pragma STDC FENV_ACCESS ON

namespace rocs {

    const double EPSIVAL = 1e-6; /* 1e-9 //DBL_EPSILON */

    const double PINF = std::numeric_limits<double>::infinity(); /* 1.0 / 0.0; */
    const double NINF = -PINF; /* -1.0 / 0.0; */
  
    const double PIIVAL = 3.14159265358979323846;
    const double PI2IVAL = 6.28318530717958647693;
    const double PIHALIVAL = 1.57079632679489661923;

    const int FUNCINVIDBIT = 0;
    const double EPSMACHINE = std::numeric_limits<double>::epsilon();
    const double ETA = std::numeric_limits<double>::min();
    const double ONEMINUSEPS = 1.0 - (1 << (FUNCINVIDBIT + 1)) * EPSMACHINE;
    const double ONEPLUSEPS = 1.0 + (1 << (FUNCINVIDBIT + 1)) * EPSMACHINE;


    /* set rounding style */
    inline void roundup(){ std::fesetround(FE_UPWARD); }
    inline void rounddown(){ std::fesetround(FE_DOWNWARD); }
    inline void roundnear(){ std::fesetround(FE_TONEAREST); }
    inline void roundzero(){ std::fesetround(FE_TOWARDZERO); }
    inline double roundup(double x) {
	return ((x < 0.0) ? ONEMINUSEPS : ONEPLUSEPS) * (x) + ETA;
    }
    inline double rounddown(double x) {
	return ((x < 0.0) ? ONEPLUSEPS : ONEMINUSEPS) * (x) - ETA;
    }

    typedef uint16_t UintSmall;

} // namespace rocs

#endif
