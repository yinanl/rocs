/**
 * walkstep.hpp
 *
 * Dynamics of different modes for humanoid robots:
 * - Prismatic Inverted Pendulum Model (PIPM)
 * - Prismatic Pendulum Model (PPM)
 * - Stop-Launch Model (SLM)
 * - Multi-contact Model (MCM)
 * - Hopping Model (HM)
 * - Sliding Model (SM)
 *
 * Created by Yinan Li on July 24, 2018.
 *
 * Hybrid Systems Group, University of Waterloo.
 */


#ifndef _walkstep_h
#define _walkstep_h

#include <vector>
#include <cmath>

#include "src/system.hpp"



/* The bounds of disturbances */
// const double D[] = {0.0, 0.0};
// const double D[] = {0.01, 0.02};
const double D[] = {0.05, 0.1};
// const double D[] = {0.07, 0.12};
// const double D[] = {0.09, 0.16};
// const double D[] = {0.1, 0.2};
// const double D[] = {0.15, 0.3};
// const double D[] = {0.2, 0.4};


/* Define 6 different modes */
enum Mode {
    PIPM,
    PPM,
    SLM,
    MCM,
    HM,
    SM,
};


/**
 * A class for one walking step.
 */
class WalkStep : public rocs::System {
public:
    Mode _md; /**< mode */
    double _xloc; /**< contact location */

    WalkStep(const char *name, const double t, const int xd, const int ud,
	     Mode md, const double xc) : System(name, t, xd, ud),
					 _md(md), _xloc(xc) {}

    void get_reach_set(std::vector<rocs::ivec> &y, const rocs::ivec &x) {    
	// std::vector<rocs::ivec> y(_ugrid._nv, _xdim);
	double d1, d2;
	rocs::Rn w;
	
	switch (_md) {
	case PIPM:
	    y[0][0] = _ugrid._data[0][0]*(x[0]-_xloc) -x[1]; // asymptotes, use min w
	    y[0][1] = -_ugrid._data[0][0]*(x[0]-_xloc) -x[1];
	    if (y[0][0].getsup()<=0 && y[0][1].getsup()<=0) { // above the asymptotes
		for (size_t i = 0; i < _ugrid._nv; ++i) {
		    w = _ugrid._data[i];
		
		    d1 = (0.5/w[0])*(exp(w[0]*_tau)-exp(-w[0]*_tau))*D[0]
			+(0.5/(w[0]*w[0]))*fabs(exp(w[0]*_tau)+exp(-w[0]*_tau)-2)*D[1];
		    d2 = 0.5*fabs(exp(w[0]*_tau)+exp(-w[0]*_tau)-2)*D[0]
			+(0.5/w[0])*(exp(w[0]*_tau)-exp(-w[0]*_tau))*D[1];

		    y[i][0] = _xloc + 0.5*(exp(w[0]*_tau)+exp(-w[0]*_tau))*(x[0]-_xloc)
			+ (0.5/w[0])*(exp(w[0]*_tau)-exp(-w[0]*_tau))*x[1]
			+ rocs::interval(-d1,d1);
		    y[i][1] = 0.5*w[0]*(exp(w[0]*_tau)-exp(-w[0]*_tau))*(x[0]-_xloc)
			+ 0.5*(exp(w[0]*_tau)+exp(-w[0]*_tau))*x[1]
			+ rocs::interval(-d2,d2);
		}
	    } else {
		y.clear();
	    }
	    break;
	    
	case PPM: // an ellipse, no asymptotes
	    for (size_t i = 0; i < _ugrid._nv; ++i) {
		w = _ugrid._data[i];
		
		d1 = fabs(sin(w[0]*_tau)/w[0])*D[0] + (1-cos(w[0]*_tau))/(w[0]*w[0])*D[1];
		d2 = (1-cos(w[0]*_tau))*D[0] + fabs(sin(w[0]*_tau)/w[0])*D[1];
		
		y[i][0] = _xloc + cos(w[0]*_tau)*(x[0]-_xloc) + sin(w[0]*_tau)/w[0]*x[1]
		    + rocs::interval(-d1,d1);
		y[i][1] = -w[0]*sin(w[0]*_tau)*(x[0]-_xloc) + cos(w[0]*_tau)*x[1]
		    + rocs::interval(-d2,d2);
	    }
	    break;

	case SLM:
	    break;

	case MCM: // a simplified version
	    for (size_t i = 0; i < _ugrid._nv; ++i) {
		w = _ugrid._data[i];
		y[i][0] = x[0] + (x[1]+rocs::interval(-D[0],D[0]))*_tau
		    + 0.5*_tau*_tau*(w[0]+rocs::interval(-D[1],D[1]));
		y[i][1] = x[1] + (rocs::interval(-D[1],D[1])+w[0])*_tau;
	    }
	    break;

	case HM: // same as MCM for a zero acceleration
	    for (size_t i = 0; i < _ugrid._nv; ++i) {
		w = _ugrid._data[i];
		y[i][0] = x[0] + (x[1]+rocs::interval(-D[0],D[0]))*_tau
		    + 0.5*_tau*_tau*rocs::interval(-D[1],D[1]);
		y[i][1] = x[1] + rocs::interval(-D[1],D[1])*_tau;
	    }
	    break;

	case SM:
	    break;
	}

	// return y;
    }// get_reach_set()
    
};


/**
 * A robust margin set for PIPM and PPM.
 */
class RobustSet {
public:
    double _xc; /**< keyframe state (_xc,_yc) */
    double _yc;
    double _wsqr; /**< w^2 for PIPM, -w^2 for PPM */
    
    double _z0; /**< for zeta computation */
    double _dx;
    
    double _ds[2]; /**< robust margins */
    double _dz[2];

    double _s; 
    double _vx0;

    RobustSet(const double x, const double vx, const double w2,
	      const double z0, const double dx,
	      const double ds[], const double dz[]) :
	_xc(x), _yc(vx), _wsqr(w2), _z0(z0), _dx(dx) {
	_ds[0] = ds[0];
	_ds[1] = ds[1];
	_dz[0] = dz[0];
	_dz[1] = dz[1];
	_s = _yc*_yc/_wsqr;
	_vx0 = sqrt(_yc*_yc+_wsqr*_dx*_dx);
    }

    /**
     * Determine if a box [x] is inside the robust set.
     * It is an inner-approximation.
     * @param x the box to be tested.
     * @return a boolean value.
     */
    bool in_robust_set(const rocs::ivec &x) {
	rocs::ivec y(2);
	y[0] = _s * (x[1]*x[1] - _yc*_yc - _wsqr*(x[0]-_xc)*(x[0]-_xc));
	y[1] = _z0 + _z0 * pow(x[1]/_vx0, _wsqr) * (x[0]-_xc)/_dx;

	return y[0].getinf()>=_ds[0] && y[0].getsup()<=_ds[1] &&
	    y[1].getinf()>=_dz[0] && y[1].getsup()<=_dz[1];
    }

    /**
     * Determine if a box [x] is inside the robust tube determined by sigma bounds.
     * It is an inner-approximation.
     * @param x the box to be tested.
     * @return a boolean value.
     */
    bool in_robust_tube(const rocs::ivec &x) {
	rocs::interval y;
	y = _s * (x[1]*x[1] - _yc*_yc - _wsqr*(x[0]-_xc)*(x[0]-_xc));

	return y.getinf()>=_ds[0] && y.getsup()<=_ds[1];
    }
};


/**
 * A robust margin set for MCM.
 */
class RobustSetMCM {
public:
    double _xc; /**< keyframe state (_xc,_vxc) */
    double _yc; 
    double _u; /**< control input */

    double _z0;
    
    double _ds[2]; /**< robust margins */
    double _dz[2];

    RobustSetMCM(const double x, const double vx, const double u,
		 const double z0, const double ds[], const double dz[]) :
	_xc(x), _yc(vx), _u(u), _z0(z0) {

	_ds[0] = ds[0];
	_ds[1] = ds[1];
	_dz[0] = dz[0];
	_dz[1] = dz[1];
    }

    /**
     * Determine if a box [x] is inside the robust set.
     * It is an inner-approximation.
     * @param x the box to be tested.
     * @return a boolean value.
     */
    bool in_robust_set(const rocs::ivec &x) {
	rocs::ivec y(2);
	y[0] = 2*_u*(x[0]-_xc) - (x[1]-_yc)*(x[1]+_yc); // sigma
	y[1] = _z0 + _z0 * (-_u*log(x[1]/_yc) - (x[0]-_xc)); // zeta

	return y[0].getinf()>=_ds[0] && y[0].getsup()<=_ds[1] &&
	    y[1].getinf()>=_dz[0] && y[1].getsup()<=_dz[1];
    }

    /**
     * Determine if a box [x] is inside the robust tube determined by sigma bounds.
     * It is an inner-approximation.
     * @param x the box to be tested.
     * @return a boolean value.
     */
    bool in_robust_tube(const rocs::ivec &x) {
	rocs::interval y;
	y = 2*_u*(x[0]-_xc) - (x[1]-_yc)*(x[1]+_yc); // sigma

	return y.getinf()>=_ds[0] && y.getsup()<=_ds[1];
    }
};


/**
 * An intermediate set of two PIPM robust tubes.
 */
class InterSetPIPMs {
public:
    RobustSet *_pipm1;
    RobustSet *_pipm2;

    InterSetPIPMs(RobustSet *r1, RobustSet *r2): _pipm1(r1), _pipm2(r2){}
    
    bool in_interset(const rocs::ivec &x) {
	rocs::ivec y(2);
	y[0] = _pipm1->_s * (x[1]*x[1] - _pipm1->_yc*_pipm1->_yc - _pipm1->_wsqr*(x[0]-_pipm1->_xc)*(x[0]-_pipm1->_xc));
	y[1] = _pipm2->_s * (x[1]*x[1] - _pipm2->_yc*_pipm2->_yc - _pipm2->_wsqr*(x[0]-_pipm2->_xc)*(x[0]-_pipm2->_xc));
	return y[0].getinf()>=_pipm1->_ds[0] && y[0].getsup()<=_pipm1->_ds[1] &&
	    y[1].getinf()>=_pipm2->_ds[0] && y[1].getsup()<=_pipm2->_ds[1];

    }
};


/**
 * An intermediate set of PIPM and MCM robust tubes.
 */
class InterSetPIPMMCM {
public:
    RobustSet *_pipm;
    RobustSetMCM *_mcm;

    InterSetPIPMMCM(RobustSet *r1,RobustSetMCM *r2):_pipm(r1), _mcm(r2) {};
    
    bool in_interset(const rocs::ivec &x) {
	rocs::ivec y(2);
	y[0] = _pipm->_s*(x[1]*x[1]-_pipm->_yc*_pipm->_yc-_pipm->_wsqr*(x[0]-_pipm->_xc)*(x[0]-_pipm->_xc));
	y[1] = 2*_mcm->_u*(x[0]-_mcm->_xc) - 3*(x[1]-2/3.0*_mcm->_yc)*(x[1]-2/3.0*_mcm->_yc) + 1.0/3.0*_mcm->_yc*_mcm->_yc;
	return y[0].getinf()>=_pipm->_ds[0] && y[0].getsup()<=_pipm->_ds[1] &&
	    y[1].getinf()>=_mcm->_ds[0] && y[1].getsup()<=_mcm->_ds[1];
    }
};


/**
 * An intermediate set of PIPM and PPM robust tubes.
 */
class InterSetPIPMPPM {
public:
    static const int N = 2;
    
    RobustSet *_pipm;
    RobustSet *_ppm;

    double _zeta_pipm[N];
    double _zeta_ppm[N];

    InterSetPIPMPPM(RobustSet *r1,RobustSet *r2, const double z1[], const double z2[]):
	_pipm(r1), _ppm(r2) {
	for (int i = 0; i < N; ++i) {
	    _zeta_pipm[i] = z1[i];
	    _zeta_ppm[i] = z2[i];
	}
    }
    
    bool in_interset_sigma(const rocs::ivec &x) {
	rocs::ivec y(2);
	y[0] = _pipm->_s*(x[1]*x[1]-_pipm->_yc*_pipm->_yc-_pipm->_wsqr*(x[0]-_pipm->_xc)*(x[0]-_pipm->_xc));
	y[1] = _ppm->_s*(x[1]*x[1]-_ppm->_yc*_ppm->_yc-_ppm->_wsqr*(x[0]-_ppm->_xc)*(x[0]-_ppm->_xc));
	return y[0].getinf()>=_pipm->_ds[0] && y[0].getsup()<=_pipm->_ds[1] &&
	    y[1].getinf()>=_ppm->_ds[0] && y[1].getsup()<=_ppm->_ds[1];
    }

    bool in_interset_zeta(const rocs::Rn x) {
	double y[2];
	y[0] = _pipm->_z0+_pipm->_z0*pow(x[1]/_pipm->_vx0,_pipm->_wsqr)*(x[0]-_pipm->_xc)/_pipm->_dx;
	y[1] = _ppm->_z0+_ppm->_z0*pow(x[1]/_ppm->_vx0,_ppm->_wsqr)*(x[0]-_ppm->_xc)/_ppm->_dx;
	return y[0]>=_zeta_pipm[0] && y[0]<=_zeta_pipm[1] && y[1]>=_zeta_ppm[0] && y[1]<=_zeta_ppm[1];
    }
};


#endif
