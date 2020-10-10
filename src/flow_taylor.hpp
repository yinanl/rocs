/**
 *  flowTaylor.hpp
 *
 *  A flowpipe class based on Taylor models.
 *
 *  Created by Yinan Li on Mar. 23, 2018.
 *
 *  Hybrid Systems Group, University of Waterloo.
 */


#ifndef _flowtaylor_h
#define _flowtaylor_h

#include <vector>
#include <iostream>

#include "interval_vector.h"
#include "definitions.h"

#include "adutils.h"
#include "FADBAD++/tadiff.h"



namespace rocs {

    /**
     * \brief A parameter class.
     *
     * It controls the precision of reachable set computation for the flowTaylor class.
     */
    class params {
    public:
	int kmax;  /**< maximum order */
	double tol;  /**< tolerance to determine the order k */
	double alpha;  /**< the error distribution factor */
	double beta; /**< the bloating factor */

	double eps;  /**< the precision control parameter */

	/**
	 * Constructors.
	 */
	params():kmax(10),tol(0.01),alpha(0.5),beta(2),eps(0.01) {}
    
	params(const int maxorder, const double tolerance,
	       const double a, const double b)
	    : kmax(maxorder), tol(tolerance), alpha(a), beta(b), eps(0.01) {}
    
    };  // class params
    

    /**
     * A Taylor model class for reachable set computation.
     */
    template<typename F, typename S, typename P>
    class flowTaylor {
    public:
	/**
	 * Constructors (with or without controls).
	 */
	flowTaylor(const P u, params *pdata,
		   const double T=0.01, const double delta=0,
		   const double K=1.0):
	    _tau(T),_delta(delta),_K(K),_kbar(1),_wbar(0),
	    _parameters(pdata),_p(F::n),_u(F::n),_xenc(F::n) {
	    
	    F(frems, yrems, u);
	    F(fterms, yterms, u);
	    construct_helper();
	    
// #ifdef LOGGING
// 	    std::cout << "LOGGING: Taylor model.\n"
// 		      << "K=" << _K
// 		      << ", log((1-alpha)*delta/K)=" << _logdel1
// 		      << ", log(tau)=" << _logdel2
// 		      << '\n';
// #endif
	}
	flowTaylor(params *pdata,
		   const double T=0.01, const double delta=0,
		   const double K=1.0):
	    _tau(T),_delta(delta),_K(K),_kbar(1),_wbar(0),
	    _parameters(pdata),_p(F::n),_u(F::n),_xenc(F::n) {
	    
	    F(frems, yrems);
	    F(fterms, yterms);
	    construct_helper();
// #ifdef LOGGING
// 	    std::cout << "LOGGING: Taylor model.\n"
// 		      << "K=" << _K
// 		      << ", log((1-alpha)*delta/K)=" << _logdel1
// 		      << ", log(tau)=" << _logdel2
// 		      << '\n';
// #endif
	}
	/**
	 * Pre-computation of some coefficients:
	 * _parameters->eps: (al*t)*del/(Ke^t),
	 * _logdel1: log((1-al)*del/K),
	 * _logdel2: log(t).
	 *
	 * If K is known a priori, should be set in constructor (see @flowTaylor).
	 */
	void construct_helper() {
	    _parameters->eps = _parameters->alpha*_delta*_tau/(_K*std::exp(_tau));
	    _logdel1 = std::log((1-_parameters->alpha) * _delta/_K);
	    _logdel2 = std::log(_tau);
	}

	/**
	 * Operations to Taylor coefficients.
	 */
	void reset_taylor_coeffs() {
	    for (int j = 0; j < F::n; ++j)
		fterms[j].reset();
	}
    
	void init_taylor_coeffs(S *x) {
	    for (int j = 0; j < F::n; ++j) {
		fterms[j].reset();
		yterms[j][0] = x[j];
	    }
	}

	void eval_taylor_coeffs(int order = 5) {
	    for (int i = 0; i < order; ++i) {
		for (int j = 0; j < F::n; ++j) {
		    fterms[j].eval(i);
		    yterms[j][i+1] = fterms[j][i]/double(i+1);
		}
	    }
	}

	void print_taylor_coeffs() {
	    for(int j = 0; j < F::n; ++j) {
		std::cout << j << " :";
		for(int i = 0; i < yterms[j].length(); ++i)
		    std::cout << " " << yterms[j][i];
		std::cout << std::endl;
	    }
	}


	/**
	 * Evaluate the bound of Taylor terms for intervals.
	 * @param x a given interval.
	 * @param k the order of the Taylor model (=10 by default).
	 */
	double eval_taylorterm_bound(const ivec &x, int k = 10) {
	    double K = 0;
	    
	    /* set the interval */
	    for (int j = 0; j < F::n; ++j) {
		frems[j].reset();
		yrems[j][0] = x[j];
	    }
	    double w0 = x.maxwidth();
	    
	    double w1 = 0;
	    for (int i = 1; i <= k; ++i) {
		w1 = 0;
		for (int j = 0; j < F::n; ++j) {
		    frems[j].eval(i-1);
		    yrems[j][i] = frems[j][i-1];
		    w1 = yrems[j][i].width() > w1 ? yrems[j][i].width() : w1;
		}
#ifdef LOGGING
		std::cout << "K for " << i << "th derivative= " << w1/w0 << '\n';
#endif
		
		K = w1/w0 > K ? w1/w0 : K;
	    }
#ifdef LOGGING
	    std::cout << "K= " << K << '\n';
#endif
	    return K;
	}

	/**
	 * Evaluate the Taylor model approximation of order k for intervals.
	 * @param[inout] y an interval storing the evaluation results.
	 * @param x a given interval.
	 * @param k the order of the Taylor model (=10 by default).
	 */
	void eval_taylor_terms(ivec &y, const ivec &x, int k = 10) {
	    /* Taylor coefficients are saved in yterms */
	    for (int j = 0; j < F::n; ++j) {
		fterms[j].reset();
		yterms[j][0] = x[j];
	    }
	    for (int i = 1; i <= k; ++i)
		eval_taylor_kthterm(y, i);
	    
	}
	void eval_taylor_kthterm(ivec &y, int k) {
	    
	    for (int j = 0; j < F::n; ++j) {
		fterms[j].eval(k-1);
		yterms[j][k] = _tau * yterms[j][k-1]/double(k);
		y[j] += yterms[j][k];  // tight enclosure: f^[i]([x])t^i
		_p[j] += yterms[j][k] * I;  // 1st part of apriori enclosure: f^[i]([x]) t^i [0,1]
	    }
	    
	}

	/**
	 * Evaluate the apriori enclosure of an interval ODE solution.
	 * @param x a given interval.
	 * @param k the order of the Taylor model (=10 by default).
	 */
	void eval_taylor_apriori(const ivec &x, int k = 10) {
	    /* results are saved in yrems */
	    for (int j = 0; j < F::n; ++j) {
		frems[j].reset();
		yrems[j][0] = x[j];
	    }
	    for (int i = 1; i <= k; ++i) {
		for (int j = 0; j < F::n; ++j) {
		    frems[j].eval(i-1);
		    yrems[j][i] = _tau * frems[j][i-1]/double(i);
		}
	    }
	}

	
	/**
	 * Compute a valid reachable set.
	 * @param[inout] y the reachable set in terms of interval.
	 * @param x a given interval.
	 * @return 1 if success and 0 otherwise.
	 */
	bool compute_reachset_valid(ivec &y, const ivec &x) {
#ifdef LOGGING
	    std::cout << "LOGGING: compute_reachset_valid.\n";
#endif
	    /* compute valid reach set first */
	    bool accept = false;
	    for (int j = 0; j < F::n; ++j) {
		y[j] = x[j];
		_p[j] = x[j];
	    }

	    /*** initialize k ***/
	    int k = 1;
	    /* compute i=0~k-1 terms of f^[i]([x])t^i */
	    for (int j = 0; j < F::n; ++j) {
	    	fterms[j].reset();
	    	yterms[j][0] = x[j];
	    }
	    double a, ratio = PINF;
	    while (k <= _parameters->kmax && ratio>_parameters->tol) {
		ratio = 0;
		for (int j = 0; j < F::n; ++j) {
		    fterms[j].eval(k-1);
// #ifdef LOGGING
// 		    std::cout << fterms[j][k-1] << ',';
// #endif
		    yterms[j][k] = _tau * fterms[j][k-1]/double(k);

		    y[j] += yterms[j][k];  // tight enclosure: f^[i]([x])t^i

		    _u[j] = yterms[j][k] * I;
		    _p[j] += _u[j];  // 1st part of apriori enclosure: f^[i]([x]) t^i [0,1]

		    a = _u[j].width()/_p[j].width();
		    ratio = ratio >= a ? ratio : a;
		}
#ifdef LOGGING
		std::cout << "\ny[" << k << "] = ";
		for (int j = 0; j < F::n; ++j) {
		    std::cout << y[j];
		    if (j < F::n-1)
			std::cout << 'x';
		    else
			std::cout << '\n';
		}
		std::cout << "\np[" << k << "] = ";
		for (int j = 0; j < F::n; ++j) {
		    std::cout << _p[j];
		    if (j < F::n-1)
			std::cout << 'x';
		    else
			std::cout << '\n';
		}
#endif
		++k;
	    } /* end while */
	    
#ifdef LOGGING
	    std::cout << "Before validation: "
		      << "k=" << k-1 << ", kbar=" << _kbar << '\n';
#endif

	    /*** verify k ***/
#ifdef LOGGING
	    ivec enc(_xenc);
#endif
	    do {/* the current k is the actual k+1 because of ++k */
		eval_taylor_apriori(_p, k);
		for (int j = 0; j < F::n; ++j) { // compute u: f^[k](p)t^k[0,1]
		    _u[j] = yrems[j][k] * I;
		    _xenc[j] = _u[j].mid() + _p[j];
		    _xenc[j] += _parameters->beta*(_u[j].width()/2)*II;
#ifdef LOGGING
		    enc[j] = _u[j].mid() + _p[j];
		    enc[j] = enc[j] +  _parameters->beta*(_u[j].width()/2)*II;
#endif
		}
#ifdef LOGGING
		std::cout << "\nxenc[" << k << "] = ";
		for (int j = 0; j < F::n; ++j) {
		    std::cout << _xenc[j];
		    if (j < F::n-1)
			std::cout << 'x';
		    else
			std::cout << '\n';
		}
		std::cout << "\nenc[" << k << "] = ";
		for (int j = 0; j < F::n; ++j) {
		    std::cout << enc[j];
		    if (j < F::n-1)
			std::cout << 'x';
		    else
			std::cout << '\n';
		}
#endif
		eval_taylor_apriori(_xenc, k);
		for (int j = 0; j < F::n; ++j) { // update u using xenc
		    _u[j] = yrems[j][k] * I;
		}
		/* verify if k is valid */
#ifdef LOGGING
		std::cout << "enclosure= " << _xenc << '\n';
		std::cout << "reachset= " << _p + _u << '\n';
#endif
		if (_xenc.isin(_p+_u)) {
		    _wbar = _xenc.maxwidth();
		    _kbar = k - 1;
		    accept = true;		    
		    break;
		} else {
		    if (k == _parameters->kmax)
			break;
		}

		/* increase k by 1: update p and y */
		eval_taylor_kthterm(y, k);		
		++k;
		
	    } while (k <= _parameters->kmax);
#ifdef LOGGING
	    std::cout << "After validation: "
		      << "k=" << k << ", kbar=" << _kbar
		      << ", accept= " << accept << '\n';
#endif
	    
	    /*** output a tight enclosure ***/
	    for (int j = 0; j < F::n; ++j)  // the remainder term (k term)
		y[j] += yrems[j][k];

	    
	    return accept;
	}

	
	/**
	 * Compute a reachable set that satisfies robustly complete condition.
	 * @param[inout] y the reachable set in terms of interval.
	 * @param x a given interval.
	 */
	void compute_reachset_robust(ivec &y) {
	    int kfac = 1;
	    for (int i = _kbar+1; i > 0; --i)
		kfac *= i;
	    
	    int k = std::ceil((_logdel1-std::log(_wbar)+std::log(kfac))/_logdel2);
	    if (_kbar < k) {
		
		for (int i = _kbar+1; i <= k; ++i) {
		    for (int j = 0; j < F::n; ++j) {
			fterms[j].eval(i-1);
			yterms[j][i] = _tau * yterms[j][i-1]/double(i);
			y[j] += yterms[j][i];
		    }
		}

		for (int i = _kbar+2; i <= k+1; ++i) {
		    for (int j = 0; j < F::n; ++j) {
			frems[j].eval(i-1);
			yrems[j][i] = _tau * frems[j][i-1]/double(i);
		    }
		}
		for (int j = 0; j < F::n; ++j) {
		    y[j] += yrems[j][k+1];
		}
	    }
	    
	}

	/**
	 * Consider robustness when width([x]) < epsilon.
	 */
	bool reachset_robust(ivec &y, const ivec &x, double eps) {
	    
	    bool accept = compute_reachset_valid(y, x);
	    
	    if (accept && x.maxwidth() < eps) 
		compute_reachset_robust(y);

	    return accept;
	}
	
    // protected:
	static const interval I;
	static const interval II;
    
	double _tau;  /**< The continuous time horizon. */
	double _delta;  /**< The bound of the perturbation. */
	double _K;  /**< The maximum operator norm of jacobian (df^[i]/dx). */

	int _kbar;  /**< The valid order (locally updated). */
	double _wbar;  /**< The width of the enclosure (locally updated). */
	
	T<S> yterms[F::n];  /**< Taylor terms of the solution (independent variables) */
	T<S> fterms[F::n];  /**< The dependent variables of yterms. */
	T<S> yrems[F::n];  /**< The remainder (evaluated by apriori enclosure). */
	T<S> frems[F::n];  /**< The dependent variables of yrems. */

	params *_parameters;  /**< The pointer to a parameters class. */

	double _logdel1, _logdel2;  /**< The variables to store intermediate info. */
	ivec _p;  /**< sum of yterms[i]*[0,1] */
	ivec _u;  /**< yrem*/
	ivec _xenc;

	/**
	 * Evaluate the local epsilon using the local a prior enclosure.
	 */
	friend double eval_epsilon(flowTaylor& f);
    };
    
    
    template<typename F, typename S, typename P>
    const interval flowTaylor<F,S,P>::I = interval(0,1);
    
    template<typename F, typename S, typename P>
    const interval flowTaylor<F,S,P>::II = interval(-1,1);

    
    template<typename F, typename S, typename P>
    double eval_epsilon(flowTaylor<F,S,P> &f) {
	
	double Kloc = f.eval_taylorterm_bound(f._xenc, f._kbar+1);
	return f._parameters->eps/Kloc;
    }

} // namespace rocs


#endif
