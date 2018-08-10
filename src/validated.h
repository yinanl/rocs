/**
 *  validated.h
 *
 *  Functions for validated numerics
 *
 *  Created by Yinan Li on Feb. 23, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */

#ifndef _validated_h_
#define _validated_h_

#include <iostream>
#include "config.h"


namespace rocs {
  
  inline double add_RNDD(const double x, const double y) {
    rounddown();
    double r = x + y;
    roundnear();
    return r;
  }

  inline double add_RNDU(const double x, const double y) {
    roundup();
    double r = x + y;
    roundnear();
    return r;
  }

  inline double add_RNDN(const double x, const double y) {
    roundnear();
    return x + y;
  }
  
  inline double add_RNDZ(const double x, const double y) {
    roundzero();
    double r = x + y;
    roundnear();
    return r;
  }
  

  inline double sub_RNDD(const double x, const double y) {
    rounddown();
    double r = x - y;
    roundnear();
    return r;
  }
  
  inline double sub_RNDU(const double x, const double y) {
    roundup();
    double r = x - y;
    roundnear();
    return r;
  }
  
  inline double sub_RNDN(const double x, const double y) {
    roundnear();
    return x - y;
  }
  
  inline double sub_RNDZ(const double x, const double y) {
    roundzero();
    double r = x - y;
    roundnear();
    return r;
  }

  
  inline double mul_RNDD(const double x, const double y) {
    rounddown();
    double r = x * y;
    roundnear();
    return r;
  }
  
  inline double mul_RNDU(const double x, const double y) {
    roundup();
    double r = x * y;
    roundnear();
    return r;
  }
  
  inline double mul_RNDN(const double x, const double y) {
    roundnear();
    return x * y;
  }
  
  inline double mul_RNDZ(const double x, const double y) {
    roundzero();
    double r = x * y;
    roundnear();
    return r;
  }

  
  inline double div_RNDD(const double x, const double y) {
    rounddown();
    double r = x / y;
    roundnear();
    return r;
  }
  
  inline double div_RNDU(const double x, const double y) {
    roundup();
    double r = x / y;
    roundnear();
    return r;
  }
  
  inline double div_RNDN(const double x, const double y) {
    roundnear();
    return x / y;
  }
  
  inline double div_RNDZ(const double x, const double y) {
    roundzero();
    double r = x / y;
    roundnear();
    return r;
  }

}  //namespace rocs


#endif
